/***********************************************************************
PointAccumulator - Helper class to read point clouds from multiple input
files into a list of temporary out-of-core octree files.
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

#include <string.h>
#include <utility>
#include <iostream>

#include "TempOctree.h"

#include "PointAccumulator.h"

/*********************************
Methods of class PointAccumulator:
*********************************/

void PointAccumulator::savePoints(void)
	{
	/* Create a temporary octree for the current in-memory point set: */
	std::cout<<std::endl<<"Storing "<<points.size()<<" points as temporary octree..."<<std::flush;
	char tofnt[1024];
	strcpy(tofnt,tempOctreeFileNameTemplate.c_str());
	TempOctree* to=new TempOctree(tofnt,maxNumPointsPerNode,&points[0],points.size());
	std::cout<<" done"<<std::endl;
	tempOctrees.push_back(to);
	
	/* Clear the point set: */
	points.clear();
	}

PointAccumulator::PointAccumulator(void)
	:maxNumCacheablePoints(~0x0U),
	 maxNumPointsPerNode(4096),
	 tempOctreeFileNameTemplate("/tmp/LidarPreprocessorTempOctreeXXXXXX")
	{
	}

PointAccumulator::~PointAccumulator(void)
	{
	/* Delete all temporary octrees: */
	for(std::vector<TempOctree*>::iterator toIt=tempOctrees.begin();toIt!=tempOctrees.end();++toIt)
		delete *toIt;
	}

void PointAccumulator::setMemorySize(size_t memorySize,unsigned int newMaxNumPointsPerNode)
	{
	/* Set the memory limit: */
	maxNumCacheablePoints=(memorySize*1024U*1024U+sizeof(LidarPoint)-1)/sizeof(LidarPoint);
	maxNumPointsPerNode=newMaxNumPointsPerNode;
	
	/* Check if the current point set is already too large: */
	if(points.size()>maxNumCacheablePoints)
		{
		/* Save the current point set: */
		savePoints();
		}
	
	/* Allocate the maximum amount of memory for the in-memory point set: */
	points.reserve(maxNumCacheablePoints);
	}

void PointAccumulator::setTempOctreeFileNameTemplate(std::string newTempOctreeFileNameTemplate)
	{
	tempOctreeFileNameTemplate=newTempOctreeFileNameTemplate;
	}

void PointAccumulator::setTransform(const PointAccumulator::ONTransform& newTransform)
	{
	if(newTransform.getTranslation()!=ONTransform::Vector::zero||newTransform.getRotation()!=ONTransform::Rotation::identity)
		{
		/* Convert the transformation to an affine transformation: */
		transform=ATransform(newTransform);
		haveTransform=true;
		}
	else
		haveTransform=false;
	}

void PointAccumulator::finishReading(void)
	{
	if(!points.empty())
		{
		/* Write the leftover in-memory points into another temporary octree: */
		savePoints();
		}
	
	/* A hackety-hack to release the point vector's allocated memory: */
	std::vector<LidarPoint> empty;
	std::swap(points,empty);
	}

void PointAccumulator::deleteTempOctrees(void)
	{
	/* Delete all temporary octrees: */
	for(std::vector<TempOctree*>::iterator toIt=tempOctrees.begin();toIt!=tempOctrees.end();++toIt)
		delete *toIt;
	tempOctrees.clear();
	}
