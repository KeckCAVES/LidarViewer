/***********************************************************************
LidarSelectionSaver - Point processor functor class to save the set of
selected points to a file
Copyright (c) 2005-2008 Oliver Kreylos

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

#ifndef LIDARSELECTIONSAVER_INCLUDED
#define LIDARSELECTIONSAVER_INCLUDED

#include <stdio.h>
#include <Misc/File.h>

#include "LidarTypes.h"
#include "LidarOctree.h"

class LidarSelectionSaver
	{
	/* Elements: */
	private:
	Misc::File selectionFile; // The file to save the selected points to
	Vector pointOffset; // Offset vector from octree coordinates to real coordinates
	
	/* Constructors and destructors: */
	public:
	LidarSelectionSaver(const char* selectionFileName,const Vector& sPointOffset)
		:selectionFile(selectionFileName,"wt"),
		 pointOffset(sPointOffset)
		{
		};
	
	/* Methods: */
	void operator()(const LidarPoint& lp) // Process the given LiDAR point
		{
		/* Write the point to the selection file: */
		Point op=lp+pointOffset;
		fprintf(selectionFile.getFilePtr(),"%.6f %.6f %.6f %u %u %u\n",op[0],op[1],op[2],(unsigned int)lp.value[0],(unsigned int)lp.value[1],(unsigned int)lp.value[2]);
		};
	};

#endif
