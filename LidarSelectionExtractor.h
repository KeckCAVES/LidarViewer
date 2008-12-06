/***********************************************************************
LidarSelectionExtractor - Point processor functor class to collect the
set of selected points into a vector.
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

#ifndef LIDARSELECTIONEXTRACTOR_INCLUDED
#define LIDARSELECTIONEXTRACTOR_INCLUDED

#include <vector>

template <class PointParam>
class LidarSelectionExtractor
	{
	/* Embedded classes: */
	public:
	typedef PointParam Point; // Type of extracted points
	
	/* Elements: */
	private:
	std::vector<Point> points; // Vector holding extracted points
	
	/* Methods: */
	public:
	void operator()(const LidarPoint& lp) // Process the given LiDAR point
		{
		/* Convert the point to the desired type and store it in the vector: */
		points.push_back(Point(lp[0],lp[1],lp[2])); // Workaround for stupid g++ 4.0.1 on Mac OS X; chooses plain copy constructor here
		};
	const std::vector<Point>& getPoints(void) const // Returns the vector of extracted points
		{
		return points;
		};
	};

#endif
