/***********************************************************************
TempOctree - Class to store points in a temporary octree for out-of-core
preprocessing of large point clouds.
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
#include <Misc/Utility.h>

#include "SplitPoints.h"

#include "TempOctree.h"

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

void TempOctree::createSubTree(TempOctree::Node& node,LidarPoint* points,size_t numPoints)
	{
	/* Check if the number of points is smaller than the maximum: */
	if(numPoints<=size_t(maxNumPointsPerNode))
		{
		/* Make this node a leaf and write the points to the temporary file: */
		node.numPoints=numPoints;
		node.pointsOffset=file.tell();
		file.write(points,numPoints);
		}
	else
		{
		/* Make this node an interior node: */
		node.numPoints=numPoints;
		node.pointsOffset=0;
		
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
				size_t leftNumPoints=splitPoints(childPoints[leftIndex],childNumPoints[leftIndex],i,node.domain.getCenter(i));
				childPoints[leftIndex+splitSize]=childPoints[leftIndex]+leftNumPoints;
				childNumPoints[leftIndex+splitSize]=childNumPoints[leftIndex]-leftNumPoints;
				childNumPoints[leftIndex]=leftNumPoints;
				}
			}
		
		/* Initialize the child nodes and create their subtrees: */
		node.children=new Node[8];
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			Node& child=node.children[childIndex];
			child.domain=Cube(node.domain,childIndex);
			createSubTree(node.children[childIndex],childPoints[childIndex],childNumPoints[childIndex]);
			}
		}
	}

void TempOctree::getPointsInCube(TempOctree::Node& node,const Cube& cube,TempOctree::LidarPointList& points)
	{
	/* Compare the cube against this node's domain: */
	int stat=node.domain.compareCube(cube);
	
	if(stat==Cube::SEPARATE)
		{
		/* Cube doesn't overlap this node, so don't bother: */
		return;
		}
	else if(node.children!=0)
		{
		/* This is an interior node, so delegate to the child nodes: */
		for(int childIndex=0;childIndex<8;++childIndex)
			getPointsInCube(node.children[childIndex],cube,points);
		}
	else if(stat==Cube::CONTAINS)
		{
		/* Add all points in this node to the list: */
		file.seekSet(node.pointsOffset);
		for(unsigned int i=0;i<node.numPoints;++i)
			{
			LidarPoint lp=file.read<LidarPoint>();
			points.push_back(lp);
			}
		}
	else
		{
		/* Add only those points in this node to the list that are inside the cube: */
		file.seekSet(node.pointsOffset);
		for(unsigned int i=0;i<node.numPoints;++i)
			{
			LidarPoint lp=file.read<LidarPoint>();
			if(cube.contains(lp))
				points.push_back(lp);
			}
		}
	}

TempOctree::TempOctree(char* fileNameTemplate,unsigned int sMaxNumPointsPerNode,LidarPoint* points,size_t numPoints)
	:tempFileName(new char[strlen(fileNameTemplate)+1]),
	 file(mkstemp(fileNameTemplate),"w+b",File::DontCare),
	 maxNumPointsPerNode(sMaxNumPointsPerNode),
	 pointBbox(Box::empty)
	{
	/* Save the temporary file name: */
	strcpy(tempFileName,fileNameTemplate);
	
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
	createSubTree(root,points,numPoints);
	}

TempOctree::~TempOctree(void)
	{
	if(tempFileName!=0)
		{
		/* Delete the temporary octree file: */
		unlink(tempFileName);
		delete[] tempFileName;
		}
	}
