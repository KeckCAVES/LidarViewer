/***********************************************************************
SpherePrimitive - Class for spheres extracted from point clouds.
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

#ifndef SPHEREPRIMITIVE_INCLUDED
#define SPHEREPRIMITIVE_INCLUDED

#include <Geometry/Point.h>
#include <GL/gl.h>
#include <GL/GLObject.h>

#include "PointPrimitive.h"

/* Forward declarations: */
namespace Comm {
class MulticastPipe;
}

class SpherePrimitive:public PointPrimitive,public GLObject
	{
	/* Embedded classes: */
	private:
	struct DataItem:public GLObject::DataItem
		{
		/* Elements: */
		public:
		GLuint displayListId; // Base ID of the display lists containing the sphere's surface and grid
		
		/* Constructors and destructors: */
		DataItem(void);
		virtual ~DataItem(void);
		};
	
	/* Elements: */
	private:
	Scalar radius; // Sphere radius
	
	/* Constructors and destructors: */
	public:
	SpherePrimitive(const LidarOctree* octree,Comm::MulticastPipe* pipe); // Creates sphere by processing selected points from the given octree; writes result to given pipe if !=0
	SpherePrimitive(Comm::MulticastPipe* pipe); // Creates sphere by reading sphere data from given pipe
	SpherePrimitive(Misc::File& file,const Vector& translation); // Reads a sphere primitive from a binary file
	
	/* Methods: */
	virtual Scalar pick(const Point& pickPoint,Scalar maxPickDistance) const;
	virtual void initContext(GLContextData& contextData) const;
	virtual void glRenderAction(GLContextData& contextData) const;
	virtual void write(Misc::File& file,const Vector& translation) const;
	};

#endif
