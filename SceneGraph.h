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

#ifndef SCENEGRAPH_INCLUDED
#define SCENEGRAPH_INCLUDED

#include <SceneGraph/GroupNode.h>

/* Forward declarations: */
class GLContextData;
namespace SceneGraph {
class GroupNode;
}

void createSceneGraph(void);
SceneGraph::GroupNode& getSceneGraphRoot(void);
void renderSceneGraph(GLContextData& contextData);
void destroySceneGraph(void);

#endif
