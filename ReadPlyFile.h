/***********************************************************************
ReadPlyFile - Function to read 3D polygon files in PLY format.
Copyright (c) 2004-2014 Oliver Kreylos

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

#ifndef READPLYFILE_INCLUDED
#define READPLYFILE_INCLUDED

/* Forward declarations: */
class PointAccumulator;

void readPlyFile(PointAccumulator& pa,const char* fileName,const char* plyColorNames[3]);

#endif
