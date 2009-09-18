/***********************************************************************
BruntonPrimitive - Class for planes extracted from point clouds, with
additional direct visualization of strike and dip angles.
Copyright (c) 2009 Oliver Kreylos

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

#ifndef BRUNTONPRIMITIVE_INCLUDED
#define BRUNTONPRIMITIVE_INCLUDED

#include <SceneGraph/GroupNode.h>

#include "PlanePrimitive.h"

class BruntonPrimitive:public PlanePrimitive
	{
	/* Elements: */
	SceneGraph::GroupNodePointer root; // Root node of the brunton visualization
	
	/* Private methods: */
	void buildBrunton(void); // Creates the brunton visualization
	
	/* Constructors and destructors: */
	public:
	BruntonPrimitive(const LidarOctree* octree,Comm::MulticastPipe* pipe); // Creates brunton by processing selected points from the given octree; writes result to given pipe if !=0
	BruntonPrimitive(Comm::MulticastPipe* pipe); // Creates brunton by reading brunton data from given pipe
	BruntonPrimitive(Misc::File& file,const Vector& translation); // Reads a brunton primitive from a binary file
	virtual ~BruntonPrimitive(void);
	};

#endif
