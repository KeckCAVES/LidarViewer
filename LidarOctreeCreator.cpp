/***********************************************************************
LidarOctreeCreator - Class to create LiDAR octrees from point clouds
using an out-of-core algorithm.
Copyright (c) 2007-2008 Oliver Kreylos

This file is part of the LiDAR processing and analysis package.

The LiDAR processing and analysis package is free software; you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

The LiDAR processing and analysis package is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with the LiDAR processing and analysis package; if not, write to the
Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307 USA
***********************************************************************/

#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <Misc/ThrowStdErr.h>
#include <Misc/PriorityHeap.h>
#include <Math/Math.h>
#include <Math/Constants.h>
#include <Geometry/ArrayKdTree.h>

#include "TempOctree.h"
#include "SplitPoints.h"

#include "LidarOctreeCreator.h"

namespace {

/************************************************************
Helper class to find nearest neighbors of points in kd-trees:
************************************************************/

class NearestNeighborFinder // Functor class to find nearest neighbors in a kd-tree
	{
	/* Elements: */
	private:
	const LidarPoint* point; // The point whose nearest neighbor we're looking for
	const LidarPoint* nnc; // The current nearest neighbor candidate
	Scalar nncDist2; // Its squared distance from the source point
	
	/* Constructors and destructors: */
	public:
	NearestNeighborFinder(const LidarPoint& sPoint) // Must be called with an actual element of the kd-tree's point array
		:point(&sPoint),
		 nnc(0),nncDist2(Math::Constants<Scalar>::max)
		{
		};
	
	/* Methods: */
	const LidarPoint& getQueryPosition(void) const
		{
		return *point;
		};
	bool operator()(const LidarPoint& node,int splitDimension)
		{
		if(&node!=point)
			{
			/* Check if this point is closer than the previous candidate: */
			Scalar dist2=Geometry::sqrDist(*point,node);
			if(nncDist2>dist2)
				{
				nnc=&node;
				nncDist2=dist2;
				}
			}
		
		/* Stop traversing the tree if the splitting plane is farther away than the current nearest neighbor: */
		return Math::sqr(node[splitDimension]-(*point)[splitDimension])<=nncDist2;
		};
	const LidarPoint* getNeighbor(void) const // Returns the nearest neighbor
		{
		return nnc;
		};
	Scalar getDistance(void) const // Returns the distance to the closest neighbor
		{
		return Math::sqrt(nncDist2);
		};
	};

class PointRemover // Functor class to remove points from a kd-tree
	{
	/* Elements: */
	private:
	Point point; // Center point of removal sphere
	Scalar radius,radius2; // Radius and squared radius of removal sphere
	const LidarPoint* nodeBase; // Base pointer to calculate node indices
	bool* removeFlags; // Array of point removal flags
	unsigned int numRemovedPoints; // Number of (newly) removed points
	
	/* Constructors and destructors: */
	public:
	PointRemover(const Point& sPoint,Scalar sRadius,const LidarPoint* sNodeBase,bool* sRemoveFlags)
		:point(sPoint),
		 radius(sRadius),radius2(Math::sqr(radius)),
		 nodeBase(sNodeBase),
		 removeFlags(sRemoveFlags),
		 numRemovedPoints(0)
		{
		};
	
	/* Methods: */
	const Point& getQueryPosition(void) const
		{
		return point;
		};
	bool operator()(const LidarPoint& node,int splitDimension)
		{
		if(Geometry::sqrDist(node,point)<=radius2)
			{
			/* Mark the point for removal: */
			unsigned int nodeIndex=&node-nodeBase;
			if(!removeFlags[nodeIndex])
				++numRemovedPoints;
			removeFlags[nodeIndex]=true;
			}
		
		/* Stop traversing the tree if the splitting plane is farther away than the removal radius: */
		return Math::abs(node[splitDimension]-point[splitDimension])<radius;
		};
	unsigned int getNumRemovedPoints(void) const
		{
		return numRemovedPoints;
		};
	};

/************************************************************************
Helper structure to subsample point sets by collapsing nearest neighbors:
************************************************************************/

struct NeighborPair
	{
	/* Elements: */
	public:
	unsigned int point;
	unsigned int neighbor;
	Scalar distance;
	
	/* Methods: */
	static bool lessEqual(const NeighborPair& n1,const NeighborPair& n2) // Comparison function needed by the priority heap
		{
		return n1.distance<=n2.distance;
		};
	};

}

/***********************************
Methods of class PointOctreeCreator:
***********************************/

void LidarOctreeCreator::writeNodePoints(LidarOctreeCreator::Node& node,unsigned int level)
	{
	/* Check if the current level is bigger than the previous maximum level in the tree: */
	if(maxLevel<level)
		{
		/* Add new temporary point file structures to the vector: */
		for(unsigned int i=maxLevel+1;i<=level;++i)
			tempPointFiles.push_back(TempPointFile());
		
		/* Remember the maximum tree level: */
		maxLevel=level;
		}
	
	/* Check if the temporary point file for this level needs to be created: */
	if(tempPointFiles[level].file==0)
		{
		/* Create a temporary point file name: */
		char fnt[1024];
		strcpy(fnt,tempPointFileNameTemplate.c_str());
		int pointFileFd=mkstemp(fnt);
		if(pointFileFd<0)
			Misc::throwStdErr("LidarOctreeCreator::writeNodePoints: Unable to open temporary point file %s",fnt);
		
		/* Create the temporary point file: */
		tempPointFiles[level].file=new File(pointFileFd,"w+b",File::DontCare);
		tempPointFiles[level].fileName=fnt;
		}
	
	/* Sort the node's points into a temporary kd-tree: */
	typedef Geometry::ArrayKdTree<LidarPoint> KdTree;
	KdTree pointTree(node.numPoints);
	LidarPoint* tpPtr=pointTree.accessPoints();
	for(unsigned int i=0;i<node.numPoints;++i)
		tpPtr[i]=node.points[i];
	pointTree.releasePoints();
	
	/* Write the node's points to the appropriate point file: */
	node.pointsOffset=tempPointFiles[level].file->tell();
	tempPointFiles[level].file->write(pointTree.accessPoints(),node.numPoints);
	}

void LidarOctreeCreator::subsample(LidarOctreeCreator::Node& node,bool deleteAllPoints)
	{
	typedef Geometry::ArrayKdTree<LidarPoint> KdTree;
	
	/* Count the total number of points in all children's point sets: */
	unsigned int totalNumPoints=0;
	for(int childIndex=0;childIndex<8;++childIndex)
		totalNumPoints+=node.children[childIndex].numPoints;
	
	/* Create a kd-tree containing all the children's points: */
	KdTree* pointTree=new KdTree(totalNumPoints);
	LidarPoint* tpPtr=pointTree->accessPoints();
	Scalar largestCollapsedDetail=Scalar(0);
	for(int childIndex=0;childIndex<8;++childIndex)
		{
		Node& child=node.children[childIndex];
		if(largestCollapsedDetail<child.detailSize)
			largestCollapsedDetail=child.detailSize;
		for(unsigned int i=0;i<child.numPoints;++i,++tpPtr)
			*tpPtr=child.points[i];
		
		/* Delete the child's point set if it has private points: */
		if(deleteAllPoints||child.children!=0)
			delete[] child.points;
		child.points=0;
		}
	pointTree->releasePoints();
	
	bool* removeFlags=0;
	unsigned int numPointsLeft;
	while(true)
		{
		/* Create a priority queue of nearest-neighbor pairs: */
		Misc::PriorityHeap<NeighborPair,NeighborPair> neighborPairs(totalNumPoints);
		
		/* Find each point's closest neighbor: */
		const LidarPoint* treePoints=pointTree->accessPoints();
		for(unsigned int i=0;i<totalNumPoints;++i)
			{
			NearestNeighborFinder nnf(pointTree->getNode(i));
			pointTree->traverseTreeDirected(nnf);
			NeighborPair np;
			np.point=i;
			np.neighbor=nnf.getNeighbor()-treePoints;
			np.distance=nnf.getDistance();
			neighborPairs.insert(np);
			}
		
		/* Remove points from the current set: */
		removeFlags=new bool[totalNumPoints];
		for(unsigned int i=0;i<totalNumPoints;++i)
			removeFlags[i]=false;
		numPointsLeft=totalNumPoints;
		while(numPointsLeft>maxNumPointsPerNode&&!neighborPairs.isEmpty())
			{
			/* Get the pair of closest neighbors: */
			const NeighborPair& cnp=neighborPairs.getSmallest();
			
			/* Only remove around points that are not already marked for removal themselves: */
			if(!removeFlags[cnp.point])
				{
				/* Mark the point temporarily so it won't be counted by the removal process: */
				removeFlags[cnp.point]=true;
				
				/* Remove points in the point's neighborhood: */
				Scalar collapseSize=cnp.distance*Scalar(1.9);
				if(largestCollapsedDetail<collapseSize)
					largestCollapsedDetail=collapseSize;
				PointRemover pr(pointTree->getNode(cnp.point),collapseSize,treePoints,removeFlags);
				pointTree->traverseTreeDirected(pr);
				
				/* Reset the point's removal flag and update the number of remaining points: */
				removeFlags[cnp.point]=false;
				numPointsLeft-=pr.getNumRemovedPoints();
				}
			
			neighborPairs.removeSmallest();
			}
		
		unsigned int actualNumPointsLeft=0;
		for(unsigned int i=0;i<totalNumPoints;++i)
			if(!removeFlags[i])
				++actualNumPointsLeft;
		if(numPointsLeft!=actualNumPointsLeft)
			std::cerr<<"Mismatch in number of points left after subsampling; "<<actualNumPointsLeft<<" vs "<<numPointsLeft<<std::endl;
		
		/* Stop subsampling if number of points is small enough (should almost always be the case): */
		if(numPointsLeft<=maxNumPointsPerNode)
			break;
		
		/* Otherwise copy the leftover points into another kd-tree and start over: */
		KdTree* newPointTree=new KdTree(numPointsLeft);
		LidarPoint* tpPtr=newPointTree->accessPoints();
		for(unsigned int i=0;i<totalNumPoints;++i)
			{
			if(!removeFlags[i])
				{
				*tpPtr=pointTree->getNode(i);
				++tpPtr;
				}
			}
		newPointTree->releasePoints();
		delete pointTree;
		pointTree=newPointTree;
		delete[] removeFlags;
		totalNumPoints=numPointsLeft;
		}
	
	/* Store the leftover points in the node's point array: */
	node.detailSize=largestCollapsedDetail;
	node.numPoints=numPointsLeft;
	node.points=new LidarPoint[numPointsLeft];
	LidarPoint* npPtr=node.points;
	for(unsigned int i=0;i<totalNumPoints;++i)
		if(!removeFlags[i])
			{
			*npPtr=pointTree->getNode(i);
			++npPtr;
			}
	delete pointTree;
	delete[] removeFlags;
	
	if(maxNumPointsPerInteriorNode<node.numPoints)
		maxNumPointsPerInteriorNode=node.numPoints;
	}

void LidarOctreeCreator::createSubTree(LidarOctreeCreator::Node& node,const Cube& nodeDomain,unsigned int level)
	{
	/* Get an upper bound on the number of points contained in this node's domain: */
	size_t numPointsBound=0;
	for(TempOctreeList::const_iterator toIt=tempOctrees.begin();toIt!=tempOctrees.end();++toIt)
		numPointsBound+=(*toIt)->boundNumPointsInCube(nodeDomain);
	
	/* Compare the estimated number of points against the allowed maximum: */
	if(numPointsBound>maxNumCachablePoints)
		{
		/* There are too many points in this domain; split the node and delegate to its children: */
		node.children=new Node[8];
		totalNumNodes+=8;
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			Node& child=node.children[childIndex];
			createSubTree(child,Cube(nodeDomain,childIndex),level+1);
			}
		
		/* Subsample the children's point sets: */
		subsample(node,true);
		
		/* Save the node's points to the temporary point file for the node's level: */
		writeNodePoints(node,level);
		}
	else if(numPointsBound>0)
		{
		/* Get the actual points contained in this node's domain: */
		//std::cout<<"Loading points from temporary octrees..."<<std::flush;
		TempOctree::LidarPointList points;
		points.reserve(maxNumCachablePoints);
		for(TempOctreeList::const_iterator toIt=tempOctrees.begin();toIt!=tempOctrees.end();++toIt)
			(*toIt)->getPointsInCube(nodeDomain,points);
		totalNumReadPoints+=points.size();
		//std::cout<<" done reading "<<points.size()<<" points"<<std::endl;
		std::cout<<"\rCreating octree... "<<Math::floor(double(totalNumReadPoints)*100.0/double(totalNumPoints)+0.5)<<"% done"<<std::flush;
		
		/* Call the subtree creation method for this node again, this time with the actual points: */
		createSubTreeWithPoints(node,nodeDomain,&points[0],points.size(),level);
		
		/* If this node is a leaf, we need to copy its points so they don't get deleted with the vector: */
		if(node.children==0)
			{
			LidarPoint* newPoints=new LidarPoint[node.numPoints];
			for(unsigned int i=0;i<node.numPoints;++i)
				newPoints[i]=node.points[i];
			node.points=newPoints;
			}
		}
	else
		{
		/* Create an empty leaf node: */
		node.detailSize=0.0f;
		node.numPoints=0;
		node.points=0;
		}
	}

void LidarOctreeCreator::createSubTreeWithPoints(LidarOctreeCreator::Node& node,const Cube& nodeDomain,LidarPoint* points,size_t numPoints,unsigned int level)
	{
	/* Check if the number of points is smaller than the maximum: */
	if(numPoints<=size_t(maxNumPointsPerNode))
		{
		/* Create a leaf node: */
		node.detailSize=0.0f; // Detail size of leaf nodes is always 0.0; they can't be subdivided anyways
		node.numPoints=(unsigned int)(numPoints);
		node.points=points;
		if(maxNumPointsPerInteriorNode<node.numPoints)
			maxNumPointsPerInteriorNode=node.numPoints;
		}
	else
		{
		/* Split the point array between the node's children: */
		LidarPoint* childPoints[8];
		size_t childNumPoints[8];
		childPoints[0]=points;
		childNumPoints[0]=numPoints;
		
		/* Split the point set along the three dimensions, according to the node's center: */
		int numSplits=1;
		int splitSize=4;
		for(int i=2;i>=0;--i,numSplits<<=1,splitSize>>=1)
			{
			int leftIndex=0;
			for(int j=0;j<numSplits;++j,leftIndex+=splitSize*2)
				{
				size_t leftNumPoints=splitPoints(childPoints[leftIndex],childNumPoints[leftIndex],i,nodeDomain.getCenter(i));
				childPoints[leftIndex+splitSize]=childPoints[leftIndex]+leftNumPoints;
				childNumPoints[leftIndex+splitSize]=childNumPoints[leftIndex]-leftNumPoints;
				childNumPoints[leftIndex]=leftNumPoints;
				}
			}
		
		/* Initialize the child nodes and create their subtrees: */
		node.children=new Node[8];
		totalNumNodes+=8;
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			Node& child=node.children[childIndex];
			createSubTreeWithPoints(node.children[childIndex],Cube(nodeDomain,childIndex),childPoints[childIndex],childNumPoints[childIndex],level+1);
			}
		
		/* Subsample the children's point sets: */
		subsample(node,false);
		}
	
	/* Save the node's points to the temporary point file: */
	writeNodePoints(node,level);
	}

void LidarOctreeCreator::calcFileOffsets(LidarOctreeCreator::Node& node,unsigned int level,LidarFile::Offset& octreeFilePos,LidarFile::Offset& dataFilePos)
	{
	if(level==0)
		{
		/* Calculate the node's offset: */
		node.octreeNodeOffset=octreeFilePos;
		octreeFilePos+=LidarFile::Offset(LidarOctreeFileNode::getFileSize());
		
		/* Calculate the node's points' offsets: */
		node.octreeDataOffset=dataFilePos;
		dataFilePos+=LidarFile::Offset(node.numPoints);
		}
	else if(node.children!=0)
		{
		/* Recurse into the node's children: */
		for(int childIndex=0;childIndex<8;++childIndex)
			calcFileOffsets(node.children[childIndex],level-1,octreeFilePos,dataFilePos);
		}
	}

void LidarOctreeCreator::writeSubtree(const LidarOctreeCreator::Node& node,unsigned int level,File& tempPointFile,LidarFile& octreeFile,LidarFile& pointsFile,LidarPoint* pointBuffer)
	{
	if(level==0)
		{
		LidarOctreeFileNode ofn;
		
		/* Find the node's children's offset: */
		ofn.childrenOffset=LidarFile::Offset(0);
		if(node.children!=0)
			{
			/* Store the offset of the node's first child: */
			ofn.childrenOffset=node.children[0].octreeNodeOffset;
			
			/* Check if the node's children have consecutive offsets (extra paranoia): */
			for(int childIndex=1;childIndex<8;++childIndex)
				if(node.children[childIndex].octreeNodeOffset!=ofn.childrenOffset+LidarFile::Offset(LidarOctreeFileNode::getFileSize()*childIndex))
					Misc::throwStdErr("LidarOctreeCreator::writeSubtree: Node offset error in node %u",node.octreeNodeOffset);
			}
		
		/* Write the node's structure: */
		ofn.detailSize=node.detailSize;
		ofn.numPoints=node.numPoints;
		ofn.dataOffset=node.octreeDataOffset;
		ofn.write(octreeFile);
		
		/* Copy the node's points from the temporary point file: */
		tempPointFile.seekSet(node.pointsOffset);
		tempPointFile.read(pointBuffer,node.numPoints);
		pointsFile.write(pointBuffer,node.numPoints);
		++numWrittenNodes;
		if(numWrittenNodes>nextNumWrittenNodesUpdate)
			{
			int percent=int(Math::floor(double(numWrittenNodes)*100.0/double(totalNumNodes)+0.5));
			std::cout<<"\b\b\b\b"<<std::setw(3)<<percent<<"%"<<std::flush;
			nextNumWrittenNodesUpdate=((2U*(percent+1U)+1U)*totalNumNodes+199U)/200U;
			}
		}
	else if(node.children!=0)
		{
		/* Recurse into the node's children: */
		for(int childIndex=0;childIndex<8;++childIndex)
			writeSubtree(node.children[childIndex],level-1,tempPointFile,octreeFile,pointsFile,pointBuffer);
		}
	}

LidarOctreeCreator::LidarOctreeCreator(size_t sMaxNumCachablePoints,unsigned int sMaxNumPointsPerNode,const LidarOctreeCreator::TempOctreeList& sTempOctrees,std::string sTempPointFileNameTemplate)
	:maxNumCachablePoints(sMaxNumCachablePoints),
	 maxNumPointsPerNode(sMaxNumPointsPerNode),
	 tempOctrees(sTempOctrees),
	 domainBox(Box::empty),
	 tempPointFileNameTemplate(sTempPointFileNameTemplate),
	 totalNumPoints(0),totalNumReadPoints(0),
	 totalNumNodes(1),
	 maxLevel(0),
	 maxNumPointsPerInteriorNode(0)
	{
	/* Calculate the total number of points and the union of all temporary octrees' bounding boxes: */
	for(TempOctreeList::const_iterator toIt=tempOctrees.begin();toIt!=tempOctrees.end();++toIt)
		{
		domainBox.addBox((*toIt)->getPointBbox());
		totalNumPoints+=(*toIt)->getTotalNumPoints();
		}
	
	/* Calculate the root node's domain: */
	rootDomain=Cube(domainBox);
	
	/* Initialize the temporary point file vector: */
	tempPointFiles.push_back(TempPointFile());
	
	/* Create the root's subtree: */
	std::cout<<"Creating octree for "<<totalNumPoints<<" points"<<std::endl;
	std::cout<<"Creating octree... 0% done"<<std::flush;
	createSubTree(root,rootDomain,0);
	std::cout<<std::endl;
	delete[] root.points;
	root.points=0;
	if(totalNumReadPoints!=totalNumPoints)
		std::cout<<"Read "<<totalNumReadPoints<<" from temporary octree files instead of "<<totalNumPoints<<std::endl;
	std::cout<<"Octree contains "<<totalNumNodes<<" nodes with up to "<<maxNumPointsPerInteriorNode<<" points per node in "<<maxLevel+1<<" resolution levels"<<std::endl;
	
	/* Calculate the octree nodes' file offsets: */
	LidarFile::Offset octreeFilePos=LidarFile::Offset(LidarOctreeFileHeader::getFileSize());
	LidarFile::Offset dataFilePos=LidarFile::Offset(0);
	for(unsigned int level=0;level<=maxLevel;++level)
		{
		std::cout<<"Processing octree level "<<level<<std::endl;
		calcFileOffsets(root,level,octreeFilePos,dataFilePos);
		}
	std::cout<<"Octree file sizes are "<<octreeFilePos<<" bytes and "<<LidarFile::Offset(LidarDataFileHeader::getFileSize())+dataFilePos*LidarFile::Offset(sizeof(LidarPoint))<<" bytes"<<std::endl;
	}

LidarOctreeCreator::~LidarOctreeCreator(void)
	{
	/* Delete the temporary point files: */
	for(TempPointFileList::iterator tpfIt=tempPointFiles.begin();tpfIt!=tempPointFiles.end();++tpfIt)
		{
		if(tpfIt->file!=0)
			{
			delete tpfIt->file;
			unlink(tpfIt->fileName.c_str());
			}
		}
	}

void LidarOctreeCreator::write(const char* lidarFileName)
	{
	/*********************************************************************
	Try creating the new LiDAR file base directory (oh, this can fail in
	so many ways...):
	*********************************************************************/
	
	/* Check if a file or directory of the given name already exists: */
	struct stat statBuffer;
	if(stat(lidarFileName,&statBuffer)==0)
		{
		/* Check if it's a directory: */
		if(S_ISDIR(statBuffer.st_mode))
			{
			/* Create a list of all files or subdirectories in the directory: */
			std::vector<std::string> files;
			DIR* dir=opendir(lidarFileName);
			if(dir!=0)
				{
				struct dirent* entry;
				while((entry=readdir(dir))!=0)
					if(strcmp(entry->d_name,".")!=0&&strcmp(entry->d_name,"..")!=0)
						{
						std::string file=lidarFileName;
						file.push_back('/');
						file.append(entry->d_name);
						files.push_back(file);
						}
				closedir(dir);
				
				/* Try deleting all files: */
				for(std::vector<std::string>::const_iterator fIt=files.begin();fIt!=files.end();++fIt)
					{
					if(unlink(fIt->c_str())<0)
						Misc::throwStdErr("LidarOctreeCreator::write: Directory %s already exists, and file %s cannot be deleted",lidarFileName,fIt->c_str());
					}
				}
			else
				Misc::throwStdErr("LidarOctreeCreator::write: Directory %s already exists but could not be opened",lidarFileName);
			}
		else
			{
			if(unlink(lidarFileName)==0)
				{
				if(mkdir(lidarFileName,0777)<0)
					Misc::throwStdErr("LidarOctreeCreator::write: Could not create LiDAR file %s",lidarFileName);
				}
			else
				Misc::throwStdErr("LidarOctreeCreator::write: File %s already exists and could not be removed",lidarFileName);
			}
		}
	else
		{
		if(mkdir(lidarFileName,0777)<0)
			Misc::throwStdErr("LidarOctreeCreator::write: Could not create LiDAR file %s",lidarFileName);
		}
	
	/* Create the octree file: */
	std::string octreeFileName=lidarFileName;
	octreeFileName.push_back('/');
	octreeFileName.append("Index");
	LidarFile octreeFile(octreeFileName.c_str(),"wb",LidarFile::LittleEndian);
	
	/* Write the octree file header: */
	LidarOctreeFileHeader ofh(rootDomain,maxNumPointsPerInteriorNode);
	ofh.write(octreeFile);
	
	/* Create the point data file: */
	std::string pointFileName=lidarFileName;
	pointFileName.push_back('/');
	pointFileName.append("Points");
	LidarFile pointFile(pointFileName.c_str(),"wb",LidarFile::LittleEndian);
	
	/* Write the point data file header: */
	LidarDataFileHeader dfh((unsigned int)(sizeof(LidarPoint)));
	dfh.write(pointFile);
	
	/* Write all tree levels to the files: */
	std::cout<<"Writing final octree files...   0%"<<std::flush;
	numWrittenNodes=0;
	nextNumWrittenNodesUpdate=(totalNumNodes+199U)/200U;
	LidarPoint* pointBuffer=new LidarPoint[maxNumPointsPerInteriorNode];
	for(int level=0;level<=maxLevel;++level)
		{
		/* Write the level's nodes: */
		writeSubtree(root,level,*tempPointFiles[level].file,octreeFile,pointFile,pointBuffer);
		
		/* Delete the level's temporary point file: */
		delete tempPointFiles[level].file;
		unlink(tempPointFiles[level].fileName.c_str());
		}
	delete[] pointBuffer;
	tempPointFiles.clear();
	std::cout<<"\b\b\b\bdone"<<std::endl;
	}
