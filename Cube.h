/***********************************************************************
Cube - Class representing axis-aligned cubes.
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

#ifndef CUBE_INCLUDED
#define CUBE_INCLUDED

#include <utility>
#include <Math/Math.h>

#include "PointOctreeFile.h"

class Cube // Class representing a cube
	{
	/* Elements: */
	private:
	Point min,max; // Points on diagonally opposite corners of the cube
	
	/* Constructors and destructors: */
	public:
	Cube(void) // Dummy constructor
		{
		};
	Cube(const Point& sMin,const Point& sMax) // Elementwise constructor
		:min(sMin),max(sMax)
		{
		};
	Cube(const Box& box) // Creates a cube containing the given box
		{
		/* Calculate the largest size of the box: */
		float boxSize=box.getSize(0);
		for(int i=1;i<3;++i)
			if(boxSize<box.getSize(i))
				boxSize=box.getSize(i);
		
		/* Calculate the smallest cube completely containing the box: */
		for(int i=0;i<3;++i)
			{
			float sizeDiff=Math::div2(boxSize-box.getSize(i));
			min[i]=box.getMin(i)-sizeDiff;
			max[i]=box.getMax(i)+sizeDiff;
			}
		};
	Cube(const Cube& parentCube,int octantIndex) // Creates a cube describing one octant of a parent cube
		:min(parentCube.min),max(parentCube.max)
		{
		/* Shift the cube to the proper octant: */
		for(int i=0;i<3;++i)
			{
			if(octantIndex&(1<<i))
				min[i]=Math::mid(min[i],max[i]);
			else
				max[i]=Math::mid(min[i],max[i]);
			}
		};
	
	/* Methods: */
	const Point& getMin(void) const
		{
		return min;
		};
	const Point& getMax(void) const
		{
		return max;
		};
	float getCenter(int dimension) const
		{
		return Math::mid(min[dimension],max[dimension]);
		};
	std::pair<bool,bool> compareCube(const Cube& other) const // Compares the cube against the given cube; returns if it overlaps and/or is contained
		{
		bool overlapped=true;
		bool contained=true;
		for(int i=0;i<3;++i)
			{
			if(max[i]<=other.min[i]||min[i]>=other.max[i])
				overlapped=false;
			if(min[i]<other.min[i]||max[i]>other.max[i])
				contained=false;
			}
		return std::pair<bool,bool>(overlapped,contained);
		};
	bool contains(const Point& points) const // Returns true if the cube contains the given point
		{
		bool result=true;
		for(int i=0;i<3;++i)
			if(points[i]<min[i]||points[i]>=max[i])
				result=false;
		return result;
		};
	};

#endif
