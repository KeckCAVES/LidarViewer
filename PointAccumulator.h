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

#ifndef POINTACCUMULATOR_INCLUDED
#define POINTACCUMULATOR_INCLUDED

#include <string>
#include <vector>
#include <Geometry/OrthonormalTransformation.h>
#include <Geometry/AffineTransformation.h>

#include "LidarTypes.h"

/* Forward declarations: */
class TempOctree;

class PointAccumulator
	{
	/* Embedded classes: */
	public:
	typedef Geometry::OrthonormalTransformation<double,3> ONTransform; // Type for rigid body transformations
	typedef Geometry::AffineTransformation<double,3> ATransform; // Type for affine transformations
	
	/* Elements: */
	private:
	size_t maxNumCacheablePoints; // Maximum number of points allowed in memory at a time
	std::vector<LidarPoint> points; // Vector holding current in-memory point set
	unsigned int maxNumPointsPerNode; // Maximum number of points per node in the temporary octrees
	std::string tempOctreeFileNameTemplate; // File name template for temporary octrees
	std::vector<TempOctree*> tempOctrees; // List of temporary octrees holding out-of-memory point sets
	bool haveTransform; // Flag if there is a current point transformation
	ATransform transform; // The current point transformation as an affine transformation
	
	/* Private methods: */
	void savePoints(void); // Saves the current in-memory point set to a temporary octree file
	
	/* Constructors and destructors: */
	public:
	PointAccumulator(void); // Creates an empty out-of-core point accumulator
	~PointAccumulator(void); // Destroys the point accumulator
	
	/* Methods: */
	size_t getMaxNumCacheablePoints(void) const // Returns the maximum number of points to be held in memory
		{
		return maxNumCacheablePoints;
		}
	unsigned int getMaxNumPointsPerNode(void) const // Returns the maximum number of points in each temporary octree node
		{
		return maxNumPointsPerNode;
		}
	void setMemorySize(size_t memorySize,unsigned int newMaxNumPointsPerNode); // Limits the point accumulator to the given amount of memory in megabytes
	void setTempOctreeFileNameTemplate(std::string newTempOctreeFileNameTemplate); // Sets the template for temporary octree file names
	void setTransform(const ONTransform& newTransform); // Sets a new point transformation
	void addPoint(const LidarPoint& lp) // Pushes a LiDAR point into the current point set
		{
		/* Check if the current in-memory point set is too big: */
		if(points.size()==maxNumCacheablePoints)
			{
			/* Save the current point set: */
			savePoints();
			}
		
		/* Store the new point: */
		if(haveTransform)
			{
			LidarPoint tlp=LidarPoint(transform.transform(ATransform::Point(lp)),lp.value);
			points.push_back(tlp);
			}
		else
			points.push_back(lp);
		}
	void finishReading(void); // Finishes reading points from source files
	std::vector<TempOctree*>& getTempOctrees(void) // Returns the list of temporary octrees
		{
		return tempOctrees;
		}
	void deleteTempOctrees(void); // Deletes the temporary octrees
	};

#endif
