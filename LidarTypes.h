/***********************************************************************
LidarTypes - Declarations of common data types in LiDAR processing and
visualization.
Copyright (c) 2008 Oliver Kreylos

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

#ifndef LIDARTYPES_INCLDUDED
#define LIDARTYPES_INCLDUDED

#include <Geometry/Point.h>
#include <Geometry/Vector.h>
#include <Geometry/Box.h>
#include <Geometry/ValuedPoint.h>

typedef float Scalar; // Scalar type for 3D points
typedef Geometry::Point<Scalar,3> Point; // Type for 3D points
typedef Geometry::Vector<Scalar,3> Vector; // Type for vectors
typedef Geometry::Box<Scalar,3> Box; // Type for axis-aligned boxes

struct Color // Structure for RGBA colors (compatible with OpenGL)
	{
	/* Embedded classes: */
	public:
	typedef unsigned char Scalar; // Scalar type for color components
	
	/* Elements: */
	private:
	Scalar rgba[4]; // Color components
	
	/* Constructors and destructors: */
	public:
	Color(void)
		{
		}
	Color(Scalar r,Scalar g,Scalar b,Scalar a =255)
		{
		rgba[0]=r;
		rgba[1]=g;
		rgba[2]=b;
		rgba[3]=a;
		}
	
	/* Methods: */
	const Scalar* getRgba(void) const // Returns color components as an array
		{
		return rgba;
		}
	Scalar* getRgba(void) // Ditto
		{
		return rgba;
		}
	Scalar operator[](int index) const // Returns a single color component
		{
		return rgba[index];
		}
	Scalar& operator[](int index) // Ditto
		{
		return rgba[index];
		}
	template <class SourceScalarParam>
	static inline Scalar clampRound(SourceScalarParam value) // Rounds and clamps a source value
		{
		value=value<SourceScalarParam(0)?SourceScalarParam(0):value;
		value=value>SourceScalarParam(255)?SourceScalarParam(255):value;
		return Scalar(value+SourceScalarParam(0.5));
		}
	template <class SourceScalarParam>
	void setRgb(const SourceScalarParam newRgba[3]) // Sets a color from source array, with clamping and rounding
		{
		for(int i=0;i<3;++i)
			rgba[i]=clampRound(newRgba[i]);
		rgba[3]=Scalar(255);
		}
	template <class SourceScalarParam>
	void setRgba(const SourceScalarParam newRgba[4]) // Sets a color from source array, with clamping and rounding
		{
		for(int i=0;i<4;++i)
			rgba[i]=clampRound(newRgba[i]);
		}
	};

typedef Geometry::ValuedPoint<Point,Color> LidarPoint; // Type for basic LiDAR data points (3D position + RGB color)

#endif
