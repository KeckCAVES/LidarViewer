/***********************************************************************
SceneGraph - Helper functions to manage and render a global scene graph.
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

#include "SceneGraph.h"

#include <Geometry/Point.h>
#include <Geometry/Vector.h>
#include <Geometry/OrthogonalTransformation.h>
#include <GL/gl.h>
#include <Vrui/Vrui.h>
#include <Vrui/SceneGraphSupport.h>

namespace {

/***************
Global elements:
***************/

SceneGraph::GroupNodePointer root; // The scene graph root node

}

/****************
Global functions:
****************/

void createSceneGraph(void)
	{
	root=new SceneGraph::GroupNode;
	}

SceneGraph::GroupNode& getSceneGraphRoot(void)
	{
	return *root;
	}

void renderSceneGraph(GLContextData& contextData)
	{
	/* Render the scene graph: */
	Vrui::renderSceneGraph(root.getPointer(),true,contextData);
	}

void destroySceneGraph(void)
	{
	root=0;
	}
