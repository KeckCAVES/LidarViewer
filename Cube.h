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

#include <Math/Math.h>

#include "LidarTypes.h"

class Cube // Class representing a cube
	{
	/* Embedded classes: */
	public:
	enum CompareResult // Enumerated type for cube comparison results
		{
		SEPARATE=0x0,OVERLAPS=0x1,CONTAINS=0x2
		};
	
	/* Elements: */
	private:
	Point min,max; // Points on diagonally opposite corners of the cube
	
	/* Constructors and destructors: */
	public:
	Cube(void) // Dummy constructor
		{
		}
	Cube(const Point& sMin,const Point& sMax) // Elementwise constructor
		:min(sMin),max(sMax)
		{
		}
	Cube(const Box& box) // Creates a cube containing the given box
		{
		/* Calculate the largest size of the box: */
		Scalar boxSize=box.getSize(0);
		for(int i=1;i<3;++i)
			if(boxSize<box.getSize(i))
				boxSize=box.getSize(i);
		
		/* Calculate the smallest cube completely containing the box: */
		for(int i=0;i<3;++i)
			{
			Scalar sizeDiff=Math::div2(boxSize-box.getSize(i));
			min[i]=box.min[i]-sizeDiff;
			max[i]=box.max[i]+sizeDiff;
			}
		}
	Cube(Box& box) // Ditto
		{
		/* Calculate the largest size of the box: */
		Scalar boxSize=box.getSize(0);
		for(int i=1;i<3;++i)
			if(boxSize<box.getSize(i))
				boxSize=box.getSize(i);
		
		/* Calculate the smallest cube completely containing the box: */
		for(int i=0;i<3;++i)
			{
			Scalar sizeDiff=Math::div2(boxSize-box.getSize(i));
			min[i]=box.min[i]-sizeDiff;
			max[i]=box.max[i]+sizeDiff;
			}
		}
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
		}
	template <class SourcePipeParam>
	Cube(SourcePipeParam& source) // Reads a cube from a file or pipe
		{
		source.read(min.getComponents(),3);
		source.read(max.getComponents(),3);
		}
	
	/* Methods: */
	static size_t getFileSize(void) // Returns the size of a cube when written to a file or pipe
		{
		return sizeof(Scalar)*3*2;
		}
	template <class SourcePipeParam>
	Cube& read(SourcePipeParam& source) // Reads a cube from a file or pipe
		{
		source.read(min.getComponents(),3);
		source.read(max.getComponents(),3);
		return *this;
		}
	template <class SinkPipeParam>
	void write(SinkPipeParam& sink) const // Writes a cube to a file or pipe
		{
		sink.write(min.getComponents(),3);
		sink.write(max.getComponents(),3);
		}
	const Point& getMin(void) const
		{
		return min;
		}
	const Point& getMax(void) const
		{
		return max;
		}
	Point getCenter(void) const
		{
		return Geometry::mid(min,max);
		}
	Scalar getCenter(int dimension) const
		{
		return Math::mid(min[dimension],max[dimension]);
		}
	int compareCube(const Cube& other) const // Compares the cube against the given cube; returns if it overlaps and/or is contained
		{
		int result=OVERLAPS|CONTAINS;
		for(int i=0;i<3;++i)
			{
			if(max[i]<=other.min[i]||min[i]>=other.max[i])
				result=SEPARATE;
			if(min[i]<other.min[i]||max[i]>other.max[i])
				result&=~CONTAINS;
			}
		return result;
		}
	int compareBox(const Box& box) const // Compares the cube against the given box; returns if it overlaps and/or is contained
		{
		int result=OVERLAPS|CONTAINS;
		for(int i=0;i<3;++i)
			{
			if(max[i]<=box.min[i]||min[i]>=box.max[i])
				result=SEPARATE;
			if(min[i]<box.min[i]||max[i]>box.max[i])
				result&=~CONTAINS;
			}
		return result;
		}
	bool contains(const Point& point) const // Returns true if the cube contains the given point
		{
		bool result=true;
		for(int i=0;i<3&&result;++i)
			result=min[i]<=point[i]&&point[i]<max[i];
		return result;
		}
	int findChild(const Point& point) const // Returns the index of the cube's octant containing the given point
		{
		int result=0x0;
		for(int i=0;i<3;++i)
			if(point[i]>=Math::mid(min[i],max[i]))
				result|=0x1<<i;
		return result;
		}
	Scalar sqrDist(const Point& point) const // Returns the squared distance from the point to the cube, or 0 if cube contains point
		{
		Scalar result=Scalar(0);
		for(int i=0;i<3;++i)
			{
			Scalar d;
			if((d=point[i]-max[i])>Scalar(0))
				result+=Math::sqr(d);
			else if((d=point[i]-min[i])<Scalar(0))
				result+=Math::sqr(d);
			}
		return result;
		}
	};

#endif
