/***********************************************************************
Primitive - Base class for primitives (planes, spheres, ...) that can be
extracted from point clouds.
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

#ifndef PRIMITIVE_INCLUDED
#define PRIMITIVE_INCLUDED

#include <GL/gl.h>
#include <GL/GLColor.h>

/* Forward declarations: */
namespace Misc {
class File;
}
namespace Geometry {
template <class ScalarParam,int dimensionParam>
class Point;
template <class ScalarParam,int dimensionParam>
class Vector;
template <class ScalarParam,int dimensionParam>
class Ray;
}
class GLContextData;

class Primitive
	{
	/* Embedded classes: */
	public:
	typedef GLColor<GLfloat,4> Color; // Type for colors with opacity values
	typedef double Scalar; // Scalar type for primitive parameters
	typedef Geometry::Point<Scalar,3> Point; // Point type for primitives
	typedef Geometry::Vector<Scalar,3> Vector; // Vector type for primitives
	
	/* Elements: */
	protected:
	size_t numPoints; // Number of points used to construct the primitive
	Scalar rms; // Root-mean square residual of the primitive with respect to its source points
	Color surfaceColor; // Color to render the primitive's surface
	Color gridColor; // Color to render the primitive's grid
	unsigned int version; // Version number of the primitive's settings
	
	/* Constructors and destructors: */
	public:
	Primitive(void)
		:numPoints(0),rms(0),
		 surfaceColor(0.6f,0.6f,0.1f,0.5f),
		 gridColor(0.2f,0.2f,0.2f),
		 version(1)
		{
		};
	virtual ~Primitive(void) // Destroys the primitive
		{
		};
	
	/* Methods: */
	void setSurfaceColor(const Color& newSurfaceColor) // Sets the primitive's surface color
		{
		surfaceColor=newSurfaceColor;
		++version;
		};
	void setGridColor(const Color& newGridColor) // Sets the primitive's grid color
		{
		gridColor=newGridColor;
		++version;
		};
	virtual Scalar pick(const Point& pickPoint,Scalar maxPickDistance) const =0; // Returns the distance from the pick point to the primitive, or the maximum pick distance, whichever is smaller
	virtual void glRenderAction(GLContextData& contextData) const =0; // Draws the primitive
	virtual void write(Misc::File& file,const Vector& translation) const =0; // Writes a primitive to a binary file
	};

#endif
