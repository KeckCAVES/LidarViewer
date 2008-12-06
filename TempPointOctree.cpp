/***********************************************************************
TempPointOctree - Class to store points in a temporary octree for
out-of-core preprocessing of large point clouds.
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

#include "LidarOctreeFile.h"
#include "SplitPoints.h"

#include "TempPointOctree.h"

/**************************************
Methods of class TempPointOctree::Node:
**************************************/

unsigned int TempPointOctree::Node::estimateNumPointsInCube(const Cube& cube) const
	{
	/* Compare the cube against this node's domain: */
	std::pair<bool,bool> stat=domain.compareCube(cube);
	
	if(!stat.first)
		{
		/* If the cube doesn't overlap this node, this node doesn't contribute points: */
		return 0;
		}
	else if(stat.second)
		{
		/* Return the number of points in this node: */
		return numPoints;
		}
	else if(children!=0)
		{
		/* Delegate the question to this node's children: */
		unsigned int result=0;
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

unsigned int TempPointOctree::Node::boundNumPointsInCube(const Cube& cube) const
	{
	/* Compare the cube against this node's domain: */
	std::pair<bool,bool> stat=domain.compareCube(cube);
	
	if(!stat.first)
		{
		/* If the cube doesn't overlap this node, this node doesn't contribute points: */
		return 0;
		}
	else if(stat.second)
		{
		/* Return the number of points in this node: */
		return numPoints;
		}
	else if(children!=0)
		{
		/* Delegate the question to this node's children: */
		unsigned int result=0;
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

/********************************
Methods of class TempPointOctree:
********************************/

void TempPointOctree::createSubTree(TempPointOctree::Node& node,OctreePoint* points,unsigned int numPoints)
	{
	/* Check if the number of points is smaller than the maximum: */
	if(numPoints<=maxNumPointsPerNode)
		{
		/* Make this node a leaf and write the points to the temporary file: */
		node.numPoints=numPoints;
		node.pointsOffset=file.tell();
		file.write(points,numPoints);
		node.children=0;
		}
	else
		{
		/* Make this node an interior node: */
		node.numPoints=numPoints;
		node.pointsOffset=0;
		
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
				unsigned int leftNumPoints=splitPoints(childPoints[leftIndex],childNumPoints[leftIndex],i,node.domain.getCenter(i));
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

void TempPointOctree::loadSubTree(TempPointOctree::Node& node,PointOctreeFile& octFile,FileOffset nodeOffset)
	{
	/* Read the node's state from the octree file: */
	LidarOctreeFileNode ofn;
	octFile.seekSet(nodeOffset);
	ofn.read(octFile);
	
	if(ofn.childrenOffset!=FileOffset(0))
		{
		/* Create an interior node: */
		node.numPoints=0;
		node.pointsOffset=FileOffset(0);
		node.children=new Node[8];
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			node.children[childIndex].domain=Cube(node.domain,childIndex);
			loadSubTree(node.children[childIndex],octFile,ofn.childrenOffset);
			node.numPoints+=node.children[childIndex].numPoints;
			ofn.childrenOffset+=FileOffset(LidarOctreeFileNode::getSize());
			}
		}
	else
		{
		/* Create a leaf node: */
		node.numPoints=ofn.numPoints;
		node.pointsOffset=ofn.pointsOffset;
		}
	}

void TempPointOctree::getPointsInCube(TempPointOctree::Node& node,const Cube& cube,std::vector<OctreePoint>& points)
	{
	/* Compare the cube against this node's domain: */
	std::pair<bool,bool> stat=node.domain.compareCube(cube);
	
	if(!stat.first)
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
	else if(stat.second)
		{
		/* Add all points in this node to the list: */
		file.seekSet(node.pointsOffset);
		for(unsigned int i=0;i<node.numPoints;++i)
			{
			OctreePoint op=file.read<OctreePoint>();
			if(!cube.contains(op))
				std::cout<<"Read out-of-box point at ("<<op[0]<<", "<<op[1]<<", "<<op[2]<<")"<<std::endl;
			points.push_back(op);
			}
		}
	else
		{
		/* Add only those points in this node to the list that are inside the cube: */
		file.seekSet(node.pointsOffset);
		for(unsigned int i=0;i<node.numPoints;++i)
			{
			OctreePoint op=file.read<OctreePoint>();
			if(cube.contains(op))
				points.push_back(op);
			}
		}
	}

const char* TempPointOctree::getObinFileName(const char* octreeFileNameStem)
	{
	static char obinFileName[1024];
	snprintf(obinFileName,sizeof(obinFileName),"%s.obin",octreeFileNameStem);
	return obinFileName;
	}

TempPointOctree::TempPointOctree(char* fileNameTemplate,unsigned int sMaxNumPointsPerNode,OctreePoint* points,unsigned int numPoints)
	:tempFileName(new char[strlen(fileNameTemplate)+1]),
	 file(mkstemp(fileNameTemplate),"w+b",PointOctreeFile::DontCare),
	 maxNumPointsPerNode(sMaxNumPointsPerNode)
	{
	/* Save the temporary file name: */
	strcpy(tempFileName,fileNameTemplate);
	
	/* Calculate the point set's bounding box: */
	pointBbox=Box::empty;
	for(unsigned int i=0;i<numPoints;++i)
		pointBbox.addPoint(points[i]);
	
	/* Set the root's domain to contain all points: */
	root.domain=Cube(pointBbox);
	
	/* Create the root node's subtree: */
	createSubTree(root,points,numPoints);
	}

TempPointOctree::TempPointOctree(const char* octreeFileNameStem,bool newOctreeFileFormat)
	:tempFileName(0),
	 file(getObinFileName(octreeFileNameStem),"rb",PointOctreeFile::LittleEndian),
	 maxNumPointsPerNode(~0x0)
	{
	/* Open the oct file: */
	char octFileName[1024];
	snprintf(octFileName,sizeof(octFileName),"%s.oct",octreeFileNameStem);
	PointOctreeFile octFile(octFileName,"rb",PointOctreeFile::LittleEndian);
	
	/* Read the oct file header: */
	PointOctreeFileHeader ofh(octFile);
	Point min,max;
	for(int i=0;i<3;++i)
		{
		min[i]=ofh.center[i]-ofh.size;
		max[i]=ofh.center[i]+ofh.size;
		}
	pointBbox=Box(min,max);
	maxNumPointsPerNode=ofh.maxNumPointsPerNode;
	
	/* Create the root node: */
	root.domain=Cube(pointBbox);
	
	/* Create the root node's subtree: */
	loadSubTree(root,octFile,octFile.tell());
	}

TempPointOctree::~TempPointOctree(void)
	{
	if(tempFileName!=0)
		{
		/* Delete the temporary octree file: */
		unlink(tempFileName);
		delete[] tempFileName;
		}
	}
