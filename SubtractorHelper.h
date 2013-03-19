/***********************************************************************
SubtractorHelper - Helper functions and classes for point subtraction
utilities.
Copyright (c) 2009-2013 Oliver Kreylos

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

#ifndef SUBTRACTORHELPER_INCLUDED
#define SUBTRACTORHELPER_INCLUDED

#include <Math/Math.h>
#include <Geometry/ArrayKdTree.h>

#include "LidarTypes.h"

/* Forward declarations: */
namespace Geometry {
template <class ScalarParam,int dimensionParam>
class Vector;
}

/* Class to find points in a kd-tree: */
class MatchingPointFinder
	{
	/* Elements: */
	private:
	Point p; // Searched point
	Scalar epsilon,epsilon2; // Maximum search distance
	bool found; // Flag whether a matching point was found
	
	/* Constructors and destructors: */
	public:
	MatchingPointFinder(const Point& sP,Scalar sEpsilon)
		:p(sP),epsilon(sEpsilon),epsilon2(Math::sqr(epsilon)),
		 found(false)
		{
		}
	
	/* Methods: */
	const Point& getQueryPosition(void) const
		{
		return p;
		}
	bool operator()(const Point& node,int splitDimension)
		{
		/* Compare node's point to current closest point: */
		if(Geometry::sqrDist(node,p)<epsilon2)
			found=true;
		
		/* Stop traversal if split plane is farther away than epsilon: */
		return epsilon>Math::abs(node[splitDimension]-p[splitDimension]);
		};
	bool isFound(void) const
		{
		return found;
		}
	};

/* Function to load an ASCII or binary point file into a kd-tree for subtraction: */
typedef Geometry::ArrayKdTree<Point> PointKdTree;
PointKdTree* loadSubtractSet(const char* fileName,const Geometry::Vector<double,3>& offset,int numThreads =1);

#endif
