/***********************************************************************
LoadPointSet - Helper function to load a 2D point set into LiDAR Viewer
by elevating the points to the implicit LiDAR point cloud surface.
Copyright (c) 2010-2011 Oliver Kreylos

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

#ifndef LOADPOINTSET_INCLUDED
#define LOADPOINTSET_INCLUDED

#include <IO/File.h>

#include "LidarTypes.h"

void loadPointSet(const char* lidarFileName,unsigned int memCacheSize,IO::FilePtr pointFile,Scalar filterRadius,int numLobes,const Vector& pointOffset);

#endif
