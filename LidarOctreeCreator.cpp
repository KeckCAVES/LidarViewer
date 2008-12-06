/***********************************************************************
LidarOctreeCreator - Class to represent Octrees of LiDAR data points to
generate multiresolution point sets.
Copyright (c) 2005-2008 Oliver Kreylos

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

#include <iostream>
#include <Math/Math.h>
#include <Geometry/Box.h>

#include "LidarOctreeCreator.h"

/*****************************************
Methods of class LidarOctreeCreator::Node:
*****************************************/

LidarOctreeCreator::Node::~Node(void)
	{
	if(children!=0)
		{
		/* Delete private point array: */
		delete[] points;
		
		/* Delete child nodes: */
		delete[] children;
		}
	}

/***********************************
Methods of class LidarOctreeCreator:
***********************************/

unsigned int LidarOctreeCreator::splitPoints(LidarPoint* points,unsigned int numPoints,int dimension,float split)
	{
	unsigned int l=0;
	unsigned int r=numPoints;
	while(l<r)
		{
		/* All points <l are <split: */
		while(l<numPoints&&points[l][dimension]<split)
			++l;
		
		/* All points >=r are >=split: */
		while(r>0&&points[r-1][dimension]>=split)
			--r;
		
		/* Swap if necessary: */
		if(l+1<r)
			{
			LidarPoint temp=points[l];
			points[l]=points[r-1];
			points[r-1]=temp;
			++l;
			--r;
			}
		}
	
	/* Return the number of points <split: */
	return l;
	}

unsigned int LidarOctreeCreator::createSubtree(const LidarOctreeCreator::Traversal& t,LidarPoint* points,unsigned int numPoints,unsigned int treeLevel)
	{
	/* Check if there are multiple points: */
	if(treeLevel>40)
		{
		if(numPoints>0)
			{
			/* Collapse all points into a single point: */
			LidarPoint average;
			for(int i=0;i<3;++i)
				average[i]=t.center[i];
			double avRgb[3]={0.0,0.0,0.0};
			for(unsigned int i=0;i<numPoints;++i)
				for(int j=0;j<3;++j)
					avRgb[j]+=double(points[i].value[j]);
			for(int j=0;j<3;++j)
				average.value[j]=GLubyte(Math::floor(avRgb[j]/double(numPoints)+0.5));
			points[0]=average;
			
			std::cout<<"Collapsing "<<numPoints<<" identical points into ("<<average[0]<<", "<<average[1]<<", "<<average[2]<<")"<<std::endl;
			
			return numPoints-1;
			}
		else
			return 0;
		}
	
	/* Check if the number of points is smaller than the maximum: */
	if(numPoints<=maxPointsPerNode)
		{
		/* Create a leaf node: */
		t.node->numPoints=numPoints;
		t.node->points=points;
		if(maxNumPointsPerInteriorNode<numPoints)
			maxNumPointsPerInteriorNode=numPoints;
		}
	else
		{
		/* Create an interior node: */
		t.node->children=new Node[8];
		numNodes+=8;
		
		/* Shuffle the point array into eight subarrays for the child nodes: */
		LidarPoint* childPoints[8];
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
				unsigned int leftNumPoints=splitPoints(childPoints[leftIndex],childNumPoints[leftIndex],i,t.center[i]);
				childPoints[leftIndex+splitSize]=childPoints[leftIndex]+leftNumPoints;
				childNumPoints[leftIndex+splitSize]=childNumPoints[leftIndex]-leftNumPoints;
				childNumPoints[leftIndex]=leftNumPoints;
				}
			}
		
		/* Create the children's subtrees: */
		float childRadius=Math::div2(t.radius);
		unsigned int childNumRemovedPoints[8];
		unsigned int totalNumRemovedPoints=0;
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			/* Create a traversal structure for the child: */
			Traversal childT;
			childT.node=&t.node->children[childIndex];
			childT.center=t.center;
			for(int i=0;i<3;++i)
				if(childIndex&(1<<i))
					childT.center[i]+=childRadius;
				else
					childT.center[i]-=childRadius;
			childT.radius=childRadius;
			
			/* Create the child's subtree: */
			childNumRemovedPoints[childIndex]=createSubtree(childT,childPoints[childIndex],childNumPoints[childIndex],treeLevel+1);
			totalNumRemovedPoints+=childNumRemovedPoints[childIndex];
			}
		
		if(totalNumRemovedPoints>0)
			{
			/* Move all removed points to the end of the point array: */
			LidarPoint* dPtr=points;
			int childIndex;
			for(childIndex=0;childNumRemovedPoints[childIndex]==0;++childIndex)
				dPtr+=childNumPoints[childIndex];
			dPtr+=childNumPoints[childIndex]-childNumRemovedPoints[childIndex];
			for(++childIndex;childIndex<8;++childIndex)
				{
				const LidarPoint* sPtr=childPoints[childIndex];
				for(unsigned int i=0;i<childNumPoints[childIndex]-childNumRemovedPoints[childIndex];++i,++sPtr,++dPtr)
					*dPtr=*sPtr;
				}
			
			/* Delete the children: */
			delete[] t.node->children;
			numNodes-=8;
			t.node->children=0;
			
			/* Pass the error condition upwards: */
			return totalNumRemovedPoints;
			}
		
		/* Create this node's subsampled point set by "randomly" picking points from the children's sets: */
		unsigned int totalNumPoints=0;
		for(int childIndex=0;childIndex<8;++childIndex)
			totalNumPoints+=t.node->children[childIndex].numPoints;
		unsigned int sampleFactor;
		unsigned int numNodePoints;
		#if 1
		sampleFactor=subsamplingFactor;
		numNodePoints=(totalNumPoints+sampleFactor-1)/sampleFactor;
		if(numNodePoints>maxPointsPerNode)
			std::cout<<"Forced to use "<<numNodePoints<<" points in interior node"<<std::endl;
		#else
		for(sampleFactor=subsamplingFactor;(numNodePoints=(totalNumPoints+sampleFactor-1)/sampleFactor)>maxPointsPerNode;++sampleFactor)
			;
		if(sampleFactor>subsamplingFactor)
			std::cout<<"Forced to use sample factor of "<<sampleFactor<<std::endl;
		#endif
		t.node->numPoints=numNodePoints;
		t.node->points=new LidarPoint[t.node->numPoints];
		LidarPoint* dPtr=t.node->points;
		unsigned int pointIndex=0;
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			/* Take every sampleFactor-th point from the current child: */
			Node* child=&t.node->children[childIndex];
			for(;pointIndex<child->numPoints;pointIndex+=sampleFactor,++dPtr)
				*dPtr=child->points[pointIndex];
			
			/* Wrap the point index such that the total number of points matches the prediction: */
			pointIndex-=child->numPoints;
			}
		if(maxNumPointsPerInteriorNode<numNodePoints)
			maxNumPointsPerInteriorNode=numNodePoints;
		}
	
	return 0;
	}

bool LidarOctreeCreator::calcFileOffsets(LidarOctreeCreator::Node* node,int level,Misc::LargeFile::Offset& octreeFilePos,Misc::LargeFile::Offset& pointsFilePos)
	{
	if(level==0)
		{
		/* Calculate the node's offset: */
		node->nodeOffset=octreeFilePos;
		octreeFilePos+=Misc::LargeFile::Offset(LidarOctreeFileNode::getSize());
		
		#if 0
		/* Calculate the node's children's offsets: */
		if(node->children!=0)
			{
			node->childrenOffset=octreeFilePos;
			octreeFilePos+=Misc::LargeFile::Offset(sizeof(LidarOctreeFileNode)*8);
			}
		else
			node->childrenOffset=Misc::LargeFile::Offset(0);
		#endif
		
		/* Calculate the node's points' offsets: */
		node->pointsOffset=pointsFilePos;
		pointsFilePos+=Misc::LargeFile::Offset(node->numPoints)*Misc::LargeFile::Offset(sizeof(LidarPoint));
		
		/* Return if the node has unprocessed children: */
		return node->children!=0;
		}
	else if(node->children!=0)
		{
		/* Recurse into the node's children: */
		bool result=false;
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			bool childResult=calcFileOffsets(&node->children[childIndex],level-1,octreeFilePos,pointsFilePos);
			result=result||childResult;
			}
		
		/* Return if any of the child nodes have non-processed children: */
		return result;
		}
	else
		return false;
	}

void LidarOctreeCreator::writeSubtree(const LidarOctreeCreator::Node* node,int level,Misc::LargeFile& octreeFile,Misc::LargeFile& pointsFile) const
	{
	if(level==0)
		{
		LidarOctreeFileNode ofn;
		
		/* Find the node's children's offset: */
		ofn.childrenOffset=Misc::LargeFile::Offset(0);
		if(node->children!=0)
			{
			Misc::LargeFile::Offset firstChildOffset=node->children[0].nodeOffset;
			for(int childIndex=0;childIndex<8;++childIndex)
				if(node->children[childIndex].nodeOffset!=firstChildOffset+Misc::LargeFile::Offset(LidarOctreeFileNode::getSize()*childIndex))
					std::cout<<"Node offset error in node "<<node->nodeOffset<<std::endl;
			ofn.childrenOffset=firstChildOffset;
			}
		
		/* Write the node's structure: */
		ofn.pointsOffset=node->pointsOffset;
		ofn.numPoints=node->numPoints;
		ofn.write(octreeFile);
		
		/* Write the node's points: */
		pointsFile.write<LidarPoint>(node->points,node->numPoints);
		}
	else if(node->children!=0)
		{
		/* Recurse into the node's children: */
		for(int childIndex=0;childIndex<8;++childIndex)
			writeSubtree(&node->children[childIndex],level-1,octreeFile,pointsFile);
		}
	}

LidarOctreeCreator::LidarOctreeCreator(LidarPoint* points,unsigned int numPoints,unsigned int sMaxPointsPerNode,unsigned int sSubsamplingFactor)
	:root(0),maxPointsPerNode(sMaxPointsPerNode),subsamplingFactor(sSubsamplingFactor),
	 numNodes(0),maxNumPointsPerInteriorNode(0),numLevels(0)
	{
	/* Calculate the point set's bounding box: */
	typedef Geometry::Box<float,3> Box;
	Box bb=Box::empty;
	for(unsigned int i=0;i<numPoints;++i)
		bb.addPoint(points[i]);
	
	/* Calculate the box' center and side length: */
	Point bbCenter;
	float bbSize=0.0f;
	for(int i=0;i<3;++i)
		{
		bbCenter[i]=Math::mid(bb.getMin(i),bb.getMax(i));
		if(bbSize<bb.getSize(i))
			bbSize=bb.getSize(i);
		}
	
	/* Create the root node and the root node traversal: */
	root=new Node;
	rootTraversal.node=root;
	rootTraversal.center=bbCenter;
	rootTraversal.radius=Math::div2(bbSize);
	
	/* Create the octree: */
	while(true)
		{
		std::cout<<"Creating octree for "<<numPoints<<" points"<<std::endl;
		numNodes=1;
		unsigned int numRemovedPoints=createSubtree(rootTraversal,points,numPoints,0);
		if(numRemovedPoints==0)
			break;
		numPoints-=numRemovedPoints;
		}
	std::cout<<std::endl;
	std::cout<<"Tree contains "<<numNodes<<" nodes with up to "<<maxNumPointsPerInteriorNode<<" points per node"<<std::endl;
	
	/* Calculate the octree nodes' file offsets: */
	Misc::LargeFile::Offset octreeFilePos=Misc::LargeFile::Offset(LidarOctreeFileHeader::getSize());
	Misc::LargeFile::Offset pointsFilePos=Misc::LargeFile::Offset(0);
	numLevels=0;
	bool oneMoreLevel;
	do
		{
		std::cout<<"Processing level "<<numLevels<<std::endl;
		oneMoreLevel=calcFileOffsets(root,numLevels,octreeFilePos,pointsFilePos);
		++numLevels;
		}
	while(oneMoreLevel);
	std::cout<<"Tree contains "<<numLevels<<" levels, file sizes are "<<octreeFilePos<<" bytes and "<<pointsFilePos<<" bytes"<<std::endl;
	}

LidarOctreeCreator::~LidarOctreeCreator(void)
	{
	/* Delete the octree: */
	delete root;
	}

void LidarOctreeCreator::write(const char* octreeFileName,const char* pointsFileName) const
	{
	/* Open the output files: */
	Misc::LargeFile octreeFile(octreeFileName,"wb",Misc::LargeFile::LittleEndian);
	Misc::LargeFile pointsFile(pointsFileName,"wb",Misc::LargeFile::LittleEndian);
	
	/* Write the octree file header: */
	LidarOctreeFileHeader ofh;
	ofh.center=rootTraversal.center;
	ofh.radius=rootTraversal.radius;
	ofh.maxNumPointsPerNode=maxNumPointsPerInteriorNode;
	ofh.write(octreeFile);
	
	/* Write all tree levels to the files: */
	for(int level=0;level<numLevels;++level)
		writeSubtree(root,level,octreeFile,pointsFile);
	}
