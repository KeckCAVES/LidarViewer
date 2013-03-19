/***********************************************************************
TempOctree - Class to store points in a temporary octree for out-of-core
preprocessing of large point clouds.
Copyright (c) 2007-2013 Oliver Kreylos

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

#include "TempOctree.h"

#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string>
#include <iostream>
#include <Misc/Utility.h>
#include <Misc/ThrowStdErr.h>

#include "SplitPoints.h"

/*********************************
Methods of class TempOctree::Node:
*********************************/

size_t TempOctree::Node::estimateNumPointsInCube(const Cube& cube) const
	{
	/* Compare the cube against this node's domain: */
	int stat=domain.compareCube(cube);
	
	if(stat==Cube::SEPARATE)
		{
		/* If the cube doesn't overlap this node, this node doesn't contribute points: */
		return 0;
		}
	else if(stat==Cube::CONTAINS)
		{
		/* If the cube contains this node, this node contributes all its points: */
		return numPoints;
		}
	else if(children!=0)
		{
		/* Delegate the question to this node's children: */
		size_t result=0;
		for(int childIndex=0;childIndex<8;++childIndex)
			result+=children[childIndex].estimateNumPointsInCube(cube);
		return result;
		}
	else
		{
		/* We would have to check each individual point, so return a conservative 0: */
		return 0;
		}
	}

size_t TempOctree::Node::boundNumPointsInCube(const Cube& cube) const
	{
	/* Compare the cube against this node's domain: */
	int stat=domain.compareCube(cube);
	
	if(stat==Cube::SEPARATE)
		{
		/* If the cube doesn't overlap this node, this node doesn't contribute points: */
		return 0;
		}
	else if(stat==Cube::CONTAINS)
		{
		/* If the cube contains this node, this node contributes all its points: */
		return numPoints;
		}
	else if(children!=0)
		{
		/* Delegate the question to this node's children: */
		size_t result=0;
		for(int childIndex=0;childIndex<8;++childIndex)
			result+=children[childIndex].boundNumPointsInCube(cube);
		return result;
		}
	else
		{
		/* We would have to check each individual point, so return a conservative upper bound: */
		return numPoints;
		}
	}

/***************************
Methods of class TempOctree:
***************************/

void TempOctree::readLidarSubtree(TempOctree::Node& node,LidarFile::Offset childrenOffset,LidarFile& indexFile)
	{
	/* Check if the node is an interior node: */
	if(childrenOffset!=0)
		{
		/* Load the node's children: */
		LidarOctreeFileNode childNodes[8];
		indexFile.setReadPosAbs(childrenOffset);
		for(int childIndex=0;childIndex<8;++childIndex)
			childNodes[childIndex].read(indexFile);
		
		/* Recursively process the node's children: */
		node.children=new Node[8];
		node.numPoints=0;
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			/* Initialize the child node: */
			node.children[childIndex].domain=Cube(node.domain,childIndex);
			node.children[childIndex].numPoints=childNodes[childIndex].numPoints;
			node.children[childIndex].pointsOffset=sizeof(LidarDataFileHeader)+childNodes[childIndex].dataOffset*sizeof(LidarPoint);
			node.children[childIndex].children=0;
			
			/* Process the child node's subtree: */
			readLidarSubtree(node.children[childIndex],childNodes[childIndex].childrenOffset,indexFile);
			
			/* Accumulate the child's number of points: */
			node.numPoints+=node.children[childIndex].numPoints;
			}
		}
	}

void* TempOctree::writerThreadMethod(void)
	{
	while(true)
		{
		/* Get the next node from the write queue: */
		Node* writeNode=0;
		{
		Threads::MutexCond::Lock writeQueueLock(writeQueueCond);
		while(writerThreadRun&&writeQueue.empty())
			writeQueueCond.wait(writeQueueLock);
		if(writeQueue.empty())
			break;
		writeNode=writeQueue.front();
		writeQueue.pop_front();
		}
		
		/* Write the node to the octree file: */
		LidarPoint* points=writeNode->points;
		writeNode->pointsOffset=file.getWritePos();
		file.write(points,writeNode->numPoints);
		}
	
	return 0;
	}

void TempOctree::createSubTree(TempOctree::Node& node)
	{
	/* Check if the number of points is smaller than the maximum: */
	if(node.numPoints<=size_t(maxNumPointsPerNode))
		{
		/* Queue up this node to be written to the octree file: */
		Threads::MutexCond::Lock writeQueueLock(writeQueueCond);
		writeQueue.push_back(&node);
		writeQueueCond.signal();
		}
	else
		{
		/* Make this node an interior node: */
		node.children=new Node[8];
		
		/* Split the point array between the node's children: */
		node.children[0].numPoints=node.numPoints;
		node.children[0].points=node.points;
		
		/* Split the point set along the three dimensions, according to the node's center: */
		int numSplits=1;
		int splitSize=4;
		for(int i=2;i>=0;--i,numSplits<<=1,splitSize>>=1)
			{
			int leftIndex=0;
			for(int j=0;j<numSplits;++j,leftIndex+=splitSize*2)
				{
				size_t leftNumPoints=splitPoints(node.children[leftIndex].points,node.children[leftIndex].numPoints,i,node.domain.getCenter(i));
				node.children[leftIndex+splitSize].points=node.children[leftIndex].points+leftNumPoints;
				node.children[leftIndex+splitSize].numPoints=node.children[leftIndex].numPoints-leftNumPoints;
				node.children[leftIndex].numPoints=leftNumPoints;
				}
			}
		
		/* Create the node's children's subtrees: */
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			Node& child=node.children[childIndex];
			child.domain=Cube(node.domain,childIndex);
			createSubTree(node.children[childIndex]);
			}
		
		node.points=0;
		}
	}

LidarPoint* TempOctree::getPointsInCube(TempOctree::Node& node,const Cube& cube,LidarPoint* points)
	{
	/* Compare the cube against this node's domain: */
	int stat=node.domain.compareCube(cube);
	
	if(stat==Cube::SEPARATE)
		{
		/* Cube doesn't overlap this node, so don't bother: */
		}
	else if(node.children!=0)
		{
		/* This is an interior node, so delegate to the child nodes: */
		for(int childIndex=0;childIndex<8;++childIndex)
			points=getPointsInCube(node.children[childIndex],cube,points);
		}
	else if(stat==Cube::CONTAINS)
		{
		/* Add all points in this node to the list: */
		file.setReadPosAbs(node.pointsOffset);
		for(unsigned int i=0;i<node.numPoints;++i)
			*(points++)=file.read<LidarPoint>();
		}
	else
		{
		/* Add only those points in this node to the list that are inside the cube: */
		file.setReadPosAbs(node.pointsOffset);
		for(unsigned int i=0;i<node.numPoints;++i)
			{
			LidarPoint lp=file.read<LidarPoint>();
			if(cube.contains(lp))
				*(points++)=lp;
			}
		}
	
	return points;
	}

namespace {

/****************
Helper functions:
****************/

int createTempFile(char* fileNameTemplate)
	{
	/* Create the temporary file: */
	int result=mkstemp(fileNameTemplate);
	
	/* Check for errors: */
	if(result<0)
		{
		int error=errno;
		Misc::throwStdErr("TempOctree::TempOctree: Error %d while creating temporary octree file %s",error,fileNameTemplate);
		}
	
	return result;
	}

}

TempOctree::TempOctree(char* fileNameTemplate,unsigned int sMaxNumPointsPerNode,LidarPoint* points,size_t numPoints)
	:tempFileName(new char[strlen(fileNameTemplate)+1]),
	 file(createTempFile(fileNameTemplate),File::ReadWrite),
	 maxNumPointsPerNode(sMaxNumPointsPerNode),
	 pointBbox(Box::empty),
	 writerThreadRun(true),writerThread(this,&TempOctree::writerThreadMethod)
	{
	/* Save the temporary file name: */
	strcpy(tempFileName,fileNameTemplate);
	
	/* Immediately unlink the temporary file; it will stay alive until the file handle is closed: */
	unlink(tempFileName);
	
	/* Calculate the point set's bounding box: */
	for(unsigned int i=0;i<numPoints;++i)
		pointBbox.addPoint(points[i]);
	
	/* Extend the bounding box by a small delta to include all points in a half-open box: */
	Point newMax=pointBbox.max;
	for(int i=0;i<3;++i)
		{
		Scalar delta=Scalar(1);
		if(newMax[i]+delta!=newMax[i])
			{
			while(newMax[i]+(delta*Scalar(0.5))!=newMax[i])
				delta*=Scalar(0.5);
			}
		else
			{
			while(newMax[i]+delta==newMax[i])
				delta*=Scalar(2);
			}
		newMax[i]+=delta;
		}
	pointBbox=Box(pointBbox.min,newMax);
	
	/* Set the root's domain to contain all points: */
	root.domain=Cube(pointBbox);
	
	/* Create the root node's subtree: */
	root.numPoints=numPoints;
	root.points=points;
	createSubTree(root);
	
	/* Wait until the octree file is finished: */
	writerThreadRun=false;
	writeQueueCond.signal();
	writerThread.join();
	file.flush();
	}

namespace {

/***********************************************************
Helper functions to load LiDAR files given a base file name:
***********************************************************/

std::string getLidarPartFileName(const char* lidarFileName,const char* partFileName)
	{
	std::string result=lidarFileName;
	result.push_back('/');
	result.append(partFileName);
	return result;
	}

}

TempOctree::TempOctree(const char* lidarFileName)
	:tempFileName(0),
	 file(getLidarPartFileName(lidarFileName,"Points").c_str(),IO::File::ReadOnly)
	{
	file.setEndianness(Misc::LittleEndian);
	
	/* Open the LiDAR file's index file: */
	LidarFile indexFile(getLidarPartFileName(lidarFileName,"Index").c_str(),IO::File::ReadOnly);
	indexFile.setEndianness(Misc::LittleEndian);
	
	/* Read the octree file header: */
	LidarOctreeFileHeader ofh(indexFile);
	
	/* Initialize the tree structure: */
	maxNumPointsPerNode=ofh.maxNumPointsPerNode;
	pointBbox=Box(ofh.domain.getMin(),ofh.domain.getMax());
	
	/* Read the root node's structure: */
	LidarOctreeFileNode rootfn;
	rootfn.read(indexFile);
	root.domain=ofh.domain;
	root.numPoints=rootfn.numPoints;
	root.pointsOffset=sizeof(LidarDataFileHeader)+rootfn.dataOffset*sizeof(LidarPoint);
	root.children=0;
	
	/* Read the entire octree structure: */
	readLidarSubtree(root,rootfn.childrenOffset,indexFile);
	}

TempOctree::~TempOctree(void)
	{
	if(tempFileName!=0)
		{
		/* Delete the temporary octree file: */
		// unlink(tempFileName);
		delete[] tempFileName;
		}
	}
