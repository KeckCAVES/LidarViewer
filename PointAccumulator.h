/***********************************************************************
PointAccumulator - Helper class to read point clouds from multiple input
files into a list of temporary out-of-core octree files.
Copyright (c) 2005-2011 Oliver Kreylos

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
#include <Geometry/Point.h>
#include <Geometry/Vector.h>
#include <Geometry/Box.h>
#include <Geometry/AffineTransformation.h>

#include "LidarTypes.h"

/* Forward declarations: */
class TempOctree;

class PointAccumulator
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Point<double,3> Point; // Type for double-valued points
	typedef Geometry::Vector<double,3> Vector; // Type for double-valued vectors
	typedef Geometry::Box<double,3> Box; // Type for double-valued axis-aligned boxes
	typedef Geometry::AffineTransformation<double,3> ATransform; // Type for affine transformations
	typedef Geometry::Point<float,3> Color; // Type for float-valued colors
	typedef Geometry::Box<float,3> ColorBox; // Type for float-valued axis-aligned boxes
	
	/* Elements: */
	private:
	size_t maxNumCacheablePoints; // Maximum number of points allowed in memory at a time
	std::vector<LidarPoint> points; // Vector holding current in-memory point set
	unsigned int maxNumPointsPerNode; // Maximum number of points per node in the temporary octrees
	std::string tempOctreeFileNameTemplate; // File name template for temporary octrees
	std::vector<TempOctree*> tempOctrees; // List of temporary octrees holding out-of-memory point sets
	bool havePointOffset; // Flag if there is a current point offset
	Vector pointOffset; // Offset vector applied to incoming points before the (optional) transformation is applied
	bool haveTransform; // Flag if there is a current point transformation
	ATransform transform; // The current point transformation as an affine transformation
	float colorMask[3]; // Color mask applied to incoming RGB color components
	Box bounds; // Spatial extents of currently added point set
	ColorBox colorBounds; // Color space extents of currently added point set
	
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
	void setPointOffset(const Vector& newPointOffset); // Sets the point offset
	void resetPointOffset(void); // Resets the point offset
	const Vector& getPointOffset(void) const // Returns the current point offset
		{
		if(havePointOffset)
			return pointOffset;
		else
			return Vector::zero;
		}
	template <class TransformationParam>
	void setTransform(const TransformationParam& newTransform) // Sets a new point transformation
		{
		/* Remember that there is a transformation now: */
		haveTransform=true;
		
		/* Convert the new transformation to an affine transformation: */
		transform=newTransform;
		}
	void resetTransform(void); // Turns off point transformations
	void setColorMask(const float newColorMask[3]); // Sets the current color mask
	void resetExtents(void); // Resets the accumulated spatial and color space extents
	void printExtents(void) const; // Prints the currently accumulated spatial and color space extents
	void addPoint(const Point& p,const Color& c) // Pushes a double-valued colored point into the current point set
		{
		/* Check if the current in-memory point set is too big: */
		if(points.size()==maxNumCacheablePoints)
			{
			/* Save the current point set: */
			savePoints();
			}
		
		/* Store the new point: */
		Point pt=p;
		if(havePointOffset)
			pt+=pointOffset;
		if(haveTransform)
			pt=transform.transform(pt);
		bounds.addPoint(pt);
		LidarPoint lp;
		lp=LidarPoint::Point(pt);
		
		/* Set the new point's color: */
		for(int i=0;i<3;++i)
			{
			float col=c[i]*colorMask[i];
			if(colorBounds.min[i]>col)
				colorBounds.min[i]=col;
			if(colorBounds.max[i]<col)
				colorBounds.max[i]=col;
			lp.value[i]=::Color::clampRound(col);
			}
		lp.value[3]=::Color::Scalar(255);
		
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
