/***********************************************************************
ProfileExtractor - Algorithm to extract straight-line profiles from 2.5D
LiDAR data.
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

#ifndef PROFILEEXTRACTOR_INCLUDED
#define PROFILEEXTRACTOR_INCLUDED

#include "LidarTypes.h"

/* Forward declarations: */
namespace Cluster {
class MulticastPipe;
}
class LidarOctree;

void extractProfile(const LidarOctree* octree,const Point& p0,const Point& p1,double segmentLength,int oversampling,double filterWidth,int numLobes,Cluster::MulticastPipe* pipe);

#endif
