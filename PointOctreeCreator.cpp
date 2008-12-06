/***********************************************************************
PointOctreeCreator - Class to create point octrees from point clouds
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
#include <unistd.h>
#include <iostream>
#include <Misc/PriorityHeap.h>
#include <Math/Math.h>
#include <Math/Constants.h>
#include <Geometry/ArrayKdTree.h>

#include "SplitPoints.h"
#include "TempPointOctree.h"

#include "LidarOctreeFile.h"

#include "PointOctreeCreator.h"

namespace {

/************************************************************
Helper class to find nearest neighbors of points in kd-trees:
************************************************************/

class NearestNeighborFinder // Functor class to find nearest neighbors in a kd-tree
	{
	/* Elements: */
	private:
	const OctreePoint* point; // The point whose nearest neighbor we're looking for
	const OctreePoint* nnc; // The current nearest neighbor candidate
	float nncDist2; // Its squared distance from the source point
	
	/* Constructors and destructors: */
	public:
	NearestNeighborFinder(const OctreePoint& sPoint) // Must be called with an actual element of the kd-tree's point array
		:point(&sPoint),
		 nnc(0),nncDist2(Math::Constants<float>::max)
		{
		};
	
	/* Methods: */
	const OctreePoint& getQueryPosition(void) const
		{
		return *point;
		};
	bool operator()(const OctreePoint& node,int splitDimension)
		{
		if(&node!=point)
			{
			/* Check if this point is closer than the previous candidate: */
			float dist2=Geometry::sqrDist(*point,node);
			if(nncDist2>dist2)
				{
				nnc=&node;
				nncDist2=dist2;
				}
			}
		
		/* Stop traversing the tree if the splitting plane is farther away than the current nearest neighbor: */
		return Math::sqr(node[splitDimension]-(*point)[splitDimension])<=nncDist2;
		};
	const OctreePoint* getNeighbor(void) const // Returns the nearest neighbor
		{
		return nnc;
		};
	float getDistance(void) const // Returns the distance to the closest neighbor
		{
		return Math::sqrt(nncDist2);
		};
	};

class PointRemover // Functor class to remove points from a kd-tree
	{
	/* Elements: */
	private:
	Point point; // Center point of removal sphere
	float radius,radius2; // Radius and squared radius of removal sphere
	const OctreePoint* nodeBase; // Base pointer to calculate node indices
	bool* removeFlags; // Array of point removal flags
	unsigned int numRemovedPoints; // Number of (newly) removed points
	
	/* Constructors and destructors: */
	public:
	PointRemover(const Point& sPoint,float sRadius,const OctreePoint* sNodeBase,bool* sRemoveFlags)
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
	bool operator()(const OctreePoint& node,int splitDimension)
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
	float distance;
	
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

void PointOctreeCreator::calcDetailSize(PointOctreeCreator::Node& node)
	{
	/* Sort the node's points into a kd-tree for fast nearest-neighbor lookup: */
	Geometry::ArrayKdTree<OctreePoint> pointTree(node.numPoints,node.points); // Copy the points for now; what does it matter
	
	/* For each point, find its nearest neighbor: */
	float minDetailSize=Math::Constants<float>::max;
	const OctreePoint* treePoints=pointTree.accessPoints();
	float minDistance=Math::Constants<float>::max;
	for(unsigned int i=0;i<node.numPoints;++i)
		{
		NearestNeighborFinder nnf(treePoints[i]);
		pointTree.traverseTreeDirected(nnf);
		float distance=nnf.getDistance();
		if(minDistance>distance)
			minDistance=distance;
		}
	
	/* Set the node's detail size to the minimum distance between any two neighbors: */
	node.detailSize=minDistance;
	}

void PointOctreeCreator::writeNodePoints(PointOctreeCreator::Node& node,unsigned int level)
	{
	/* Check if the current level is bigger than the previous maximum level in the tree: */
	if(maxLevel<level)
		{
		/* Add new temporary point file structures to the vector: */
		for(unsigned int i=maxLevel+1;i<=level;++i)
			pointFiles.push_back(PointFile());
		
		/* Remember the maximum tree level: */
		maxLevel=level;
		}
	
	/* Check if the temporary point file for this level needs to be created: */
	if(pointFiles[level].file==0)
		{
		char fnt[1024];
		strcpy(fnt,fileNameTemplate.c_str());
		pointFiles[level].file=new PointOctreeFile(mkstemp(fnt),"w+b",PointOctreeFile::DontCare);
		pointFiles[level].fileName=fnt;
		}
	
	/* Write the node's points to the appropriate point file: */
	node.pointsOffset=pointFiles[level].file->tell();
	pointFiles[level].file->write(node.points,node.numPoints);
	}

#if 1

void PointOctreeCreator::subsample(PointOctreeCreator::Node& node,bool deleteAllPoints)
	{
	typedef Geometry::ArrayKdTree<OctreePoint> KdTree;
	
	/* Count the total number of points in all children's point sets: */
	unsigned int totalNumPoints=0;
	for(int childIndex=0;childIndex<8;++childIndex)
		totalNumPoints+=node.children[childIndex].numPoints;
	
	/* Create a kd-tree containing all the children's points: */
	KdTree* pointTree=new KdTree(totalNumPoints);
	OctreePoint* tpPtr=pointTree->accessPoints();
	float largestCollapsedDetail=0.0f;
	for(int childIndex=0;childIndex<8;++childIndex)
		{
		Node& child=node.children[childIndex];
		if(largestCollapsedDetail<child.detailSize)
			largestCollapsedDetail=child.detailSize;
		for(unsigned int i=0;i<child.numPoints;++i,++tpPtr)
			*tpPtr=child.points[i];
		
		/* Delete the child's point set if it has private points: */
		if(deleteAllPoints||child.children!=0)
			{
			delete[] child.points;
			child.points=0;
			}
		}
	pointTree->releasePoints();
	
	bool* removeFlags=0;
	unsigned int numPointsLeft;
	while(true)
		{
		/* Create a priority queue of nearest-neighbor pairs: */
		Misc::PriorityHeap<NeighborPair,NeighborPair> neighborPairs(totalNumPoints);
		
		/* Find each point's closest neighbor: */
		const OctreePoint* treePoints=pointTree->accessPoints();
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
				if(largestCollapsedDetail<cnp.distance*1.9f)
					largestCollapsedDetail=cnp.distance*1.9f;
				PointRemover pr(pointTree->getNode(cnp.point),cnp.distance*1.9f,treePoints,removeFlags);
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
		OctreePoint* tpPtr=newPointTree->accessPoints();
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
	node.points=new OctreePoint[numPointsLeft];
	OctreePoint* npPtr=node.points;
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

#else

void PointOctreeCreator::subsample(PointOctreeCreator::Node& node,bool deleteAllPoints)
	{
	/* Count the total number of points in all children's point sets: */
	unsigned int totalNumPoints=0;
	for(int childIndex=0;childIndex<8;++childIndex)
		totalNumPoints+=node.children[childIndex].numPoints;
	
	/* Calculate the appropriate sampling factor and this node's number of points: */
	unsigned int sampleFactor;
	unsigned int numNodePoints;
	#if 0
	sampleFactor=subsamplingFactor;
	numNodePoints=(totalNumPoints+sampleFactor-1)/sampleFactor;
	#if 0
	if(numNodePoints>maxNumPointsPerNode)
		std::cout<<"Forced to use "<<numNodePoints<<" points in interior node"<<std::endl;
	#endif
	#else
	for(sampleFactor=subsamplingFactor;(numNodePoints=(totalNumPoints+sampleFactor-1)/sampleFactor)>maxNumPointsPerNode;++sampleFactor)
		;
	#if 0
	if(sampleFactor>subsamplingFactor)
		std::cout<<"Forced to use sample factor of "<<sampleFactor<<std::endl;
	#endif
	#endif
	
	/* Subsample the children's point sets: */
	node.numPoints=numNodePoints;
	node.points=new OctreePoint[node.numPoints];
	OctreePoint* dPtr=node.points;
	unsigned int pointIndex=0;
	for(int childIndex=0;childIndex<8;++childIndex)
		{
		/* Take every sampleFactor-th point from the current child: */
		Node& child=node.children[childIndex];
		for(;pointIndex<child.numPoints;pointIndex+=sampleFactor,++dPtr)
			*dPtr=child.points[pointIndex];
		
		/* Wrap the point index such that the total number of points matches the prediction: */
		pointIndex-=child.numPoints;
		
		/* Delete the child's point set if it has private points: */
		if(deleteAllPoints||child.children!=0)
			{
			delete[] child.points;
			child.points=0;
			}
		}
	
	if(maxNumPointsPerInteriorNode<node.numPoints)
		maxNumPointsPerInteriorNode=node.numPoints;
	}

#endif

void PointOctreeCreator::createSubTree(PointOctreeCreator::Node& node,const Cube& nodeDomain,unsigned int level)
	{
	/* Get an upper bound on the number of points contained in this node's domain: */
	unsigned int numPointsBound=0;
	for(OctreeList::const_iterator oIt=octrees.begin();oIt!=octrees.end();++oIt)
		numPointsBound+=(*oIt)->boundNumPointsInCube(nodeDomain);
	
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
		
		/* Save the node's points to the temporary octree file: */
		writeNodePoints(node,level);
		}
	else if(numPointsBound>0)
		{
		/* Get the actual points contained in this node's domain: */
		//std::cout<<"Loading points from temporary octrees..."<<std::flush;
		std::vector<OctreePoint> points;
		points.reserve(maxNumCachablePoints);
		for(OctreeList::const_iterator oIt=octrees.begin();oIt!=octrees.end();++oIt)
			(*oIt)->getPointsInCube(nodeDomain,points);
		totalNumReadPoints+=points.size();
		//std::cout<<" done reading "<<points.size()<<" points"<<std::endl;
		std::cout<<"\rCreating octree... "<<Math::floor(double(totalNumReadPoints)*100.0/double(totalNumPoints)+0.5)<<"% done"<<std::flush;
		
		/* Call the subtree creation method for this node again, this time with the actual points: */
		createSubTree(node,nodeDomain,&points[0],points.size(),level);
		
		/* If this node is a leaf, we need to copy its points so they don't get deleted with the vector: */
		if(node.children==0)
			{
			OctreePoint* newPoints=new OctreePoint[node.numPoints];
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

void PointOctreeCreator::createSubTree(PointOctreeCreator::Node& node,const Cube& nodeDomain,OctreePoint* points,unsigned int numPoints,unsigned int level)
	{
	/* Check if the number of points is smaller than the maximum: */
	if(numPoints<=maxNumPointsPerNode)
		{
		/* Create a leaf node: */
		node.detailSize=0.0f; // Detail size of leaf nodes is always 0.0; they can't be subdivided anyways
		node.numPoints=numPoints;
		node.points=points;
		if(maxNumPointsPerInteriorNode<numPoints)
			maxNumPointsPerInteriorNode=numPoints;
		}
	else
		{
		/* Split the point array between the node's children: */
		OctreePoint* childPoints[8];
		unsigned int childNumPoints[8];
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
				unsigned int leftNumPoints=splitPoints(childPoints[leftIndex],childNumPoints[leftIndex],i,nodeDomain.getCenter(i));
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
			createSubTree(node.children[childIndex],Cube(nodeDomain,childIndex),childPoints[childIndex],childNumPoints[childIndex],level+1);
			}
		
		/* Subsample the children's point sets: */
		subsample(node,false);
		}
	
	/* Save the node's points to the temporary octree file: */
	writeNodePoints(node,level);
	}

void PointOctreeCreator::calcFileOffsets(PointOctreeCreator::Node& node,unsigned int level,FileOffset& octreeFilePos,FileOffset& pointsFilePos)
	{
	if(level==0)
		{
		/* Calculate the node's offset: */
		node.octreeNodeOffset=octreeFilePos;
		octreeFilePos+=FileOffset(LidarOctreeFileNode::getSize());
		
		/* Calculate the node's points' offsets: */
		node.octreePointsOffset=pointsFilePos;
		pointsFilePos+=FileOffset(node.numPoints)*FileOffset(sizeof(OctreePoint));
		}
	else if(node.children!=0)
		{
		/* Recurse into the node's children: */
		for(int childIndex=0;childIndex<8;++childIndex)
			calcFileOffsets(node.children[childIndex],level-1,octreeFilePos,pointsFilePos);
		}
	}

void PointOctreeCreator::writeSubtree(const PointOctreeCreator::Node& node,unsigned int level,PointOctreeFile& tempPointFile,PointOctreeFile& octreeFile,PointOctreeFile& pointsFile,OctreePoint* pointBuffer)
	{
	if(level==0)
		{
		LidarOctreeFileNode ofn;
		
		/* Find the node's children's offset: */
		ofn.childrenOffset=FileOffset(0);
		if(node.children!=0)
			{
			FileOffset firstChildOffset=node.children[0].octreeNodeOffset;
			for(int childIndex=0;childIndex<8;++childIndex)
				if(node.children[childIndex].octreeNodeOffset!=firstChildOffset+FileOffset(LidarOctreeFileNode::getSize()*childIndex))
					std::cout<<"Node offset error in node "<<node.octreeNodeOffset<<std::endl;
			ofn.childrenOffset=firstChildOffset;
			}
		
		/* Write the node's structure: */
		ofn.detailSize=node.detailSize;
		ofn.pointsOffset=node.octreePointsOffset;
		ofn.numPoints=node.numPoints;
		ofn.write(octreeFile);
		
		/* Copy the node's points from the temporary point file: */
		tempPointFile.seekSet(node.pointsOffset);
		tempPointFile.read(pointBuffer,node.numPoints);
		pointsFile.write(pointBuffer,node.numPoints);
		++numWrittenNodes;
		if(numWrittenNodes>nextNumWrittenNodesUpdate)
			{
			int percent=int(Math::floor(double(numWrittenNodes)*100.0/double(totalNumNodes)+0.5));
			std::cout<<"\rWriting final octree files... "<<percent<<"% done"<<std::flush;
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

PointOctreeCreator::PointOctreeCreator(const char* sFileNameTemplate,const PointOctreeCreator::OctreeList& sOctrees,unsigned int sMaxNumCachablePoints,unsigned int sMaxNumPointsPerNode,unsigned int sSubsamplingFactor)
	:octrees(sOctrees),
	 maxNumCachablePoints(sMaxNumCachablePoints),
	 maxNumPointsPerNode(sMaxNumPointsPerNode),subsamplingFactor(sSubsamplingFactor),
	 fileNameTemplate(sFileNameTemplate),
	 totalNumReadPoints(0),
	 totalNumNodes(1),maxLevel(0),maxNumPointsPerInteriorNode(0)
	{
	/* Calculate the union of all octrees' bounding boxes: */
	totalNumPoints=0;
	domainBox=Box::empty;
	for(std::vector<TempPointOctree*>::const_iterator oIt=octrees.begin();oIt!=octrees.end();++oIt)
		{
		totalNumPoints+=(*oIt)->getTotalNumPoints();
		domainBox.addBox((*oIt)->getPointBbox());
		}
	
	/* Calculate the root node's domain: */
	rootDomain=Cube(domainBox);
	
	/* Initialize the temporary point file vector: */
	pointFiles.push_back(PointFile());
	
	/* Create the root's subtree: */
	std::cout<<"Creating octree for "<<totalNumPoints<<" points"<<std::endl;
	std::cout<<"Creating octree... 0% done"<<std::flush;
	createSubTree(root,rootDomain,0);
	delete[] root.points;
	std::cout<<std::endl;
	if(totalNumReadPoints!=totalNumPoints)
		std::cout<<"Read "<<totalNumReadPoints<<" from temporary octree files instead of "<<totalNumPoints<<std::endl;
	std::cout<<"Octree contains "<<totalNumNodes<<" nodes with up to "<<maxNumPointsPerInteriorNode<<" points per node in "<<maxLevel+1<<" resolution levels"<<std::endl;
	
	/* Calculate the octree nodes' file offsets: */
	FileOffset octreeFilePos=FileOffset(LidarOctreeFileHeader::getSize());
	FileOffset pointsFilePos=FileOffset(0);
	for(unsigned int level=0;level<=maxLevel;++level)
		{
		std::cout<<"Processing octree level "<<level<<std::endl;
		calcFileOffsets(root,level,octreeFilePos,pointsFilePos);
		}
	std::cout<<"Octree file sizes are "<<octreeFilePos<<" bytes and "<<pointsFilePos<<" bytes"<<std::endl;
	}

PointOctreeCreator::~PointOctreeCreator(void)
	{
	/* Delete the temporary point files: */
	for(std::vector<PointFile>::iterator pfIt=pointFiles.begin();pfIt!=pointFiles.end();++pfIt)
		{
		if(pfIt->file!=0)
			{
			delete pfIt->file;
			unlink(pfIt->fileName.c_str());
			}
		}
	}

void PointOctreeCreator::write(const char* octreeFileName,const char* pointsFileName)
	{
	/* Open the output files: */
	PointOctreeFile octreeFile(octreeFileName,"wb",PointOctreeFile::LittleEndian);
	PointOctreeFile pointsFile(pointsFileName,"wb",PointOctreeFile::LittleEndian);
	
	/* Write the octree file header: */
	PointOctreeFileHeader ofh;
	ofh.center=Geometry::mid(rootDomain.getMin(),rootDomain.getMax());
	ofh.size=0.0f;
	for(int i=0;i<3;++i)
		if(ofh.size<rootDomain.getMax()[i]-rootDomain.getMin()[i])
			ofh.size=rootDomain.getMax()[i]-rootDomain.getMin()[i];
	ofh.size=Math::div2(ofh.size);
	ofh.maxNumPointsPerNode=maxNumPointsPerInteriorNode;
	ofh.write(octreeFile);
	
	/* Write all tree levels to the files: */
	std::cout<<"\rWriting final octree files... 0% done"<<std::flush;
	numWrittenNodes=0;
	nextNumWrittenNodesUpdate=(totalNumNodes+199U)/200U;
	OctreePoint* pointBuffer=new OctreePoint[maxNumPointsPerInteriorNode];
	for(int level=0;level<=maxLevel;++level)
		{
		/* Write the level's nodes: */
		writeSubtree(root,level,*pointFiles[level].file,octreeFile,pointsFile,pointBuffer);
		
		/* Delete the level's temporary point file: */
		delete pointFiles[level].file;
		unlink(pointFiles[level].fileName.c_str());
		}
	delete[] pointBuffer;
	pointFiles.clear();
	}
