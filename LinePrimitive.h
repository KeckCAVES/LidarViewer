/***********************************************************************
LinePrimitive - Class for lines extracted from point clouds by
intersecting two plane primitives.
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

#ifndef LINEPRIMITIVE_INCLUDED
#define LINEPRIMITIVE_INCLUDED

#include <Geometry/Point.h>
#include <Geometry/Vector.h>

#include "Primitive.h"

/* Forward declarations: */
namespace Comm {
class MulticastPipe;
}
class PlanePrimitive;

class LinePrimitive:public Primitive
	{
	/* Elements: */
	protected:
	Point center; // Point on the line, halfway between the lower and upper boundary
	Vector axis; // Normalized line direction
	Scalar length; // Line length
	
	/* Constructors and destructors: */
	protected:
	LinePrimitive(void) // Dummy constructor
		{
		};
	public:
	LinePrimitive(const PlanePrimitive* p1,const PlanePrimitive* p2,Comm::MulticastPipe* pipe); // Creates line primitive by intersecting the two given plane primitives; writes result to given pipe if !=0
	LinePrimitive(Comm::MulticastPipe* pipe); // Creates cylinder by reading cylinder data from given pipe
	LinePrimitive(Misc::File& file,const Vector& translation); // Reads a cylinder primitive from a binary file
	
	/* Methods: */
	virtual Scalar pick(const Point& pickPoint,Scalar maxPickDistance) const;
	virtual void glRenderAction(GLContextData& contextData) const;
	virtual void write(Misc::File& file,const Vector& translation) const;
	const Point& getCenter(void) const // Returns the line's center point
		{
		return center;
		};
	const Vector& getAxis(void) const // Returns the line's axis direction
		{
		return axis;
		};
	};

#endif
