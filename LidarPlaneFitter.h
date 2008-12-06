/***********************************************************************
LidarPlaneFitter - Point processor functor class to extend a previously
extracted plane to include all points used in its extraction.
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

#ifndef LIDARPLANEFITTER_INCLUDED
#define LIDARPLANEFITTER_INCLUDED

#include <Geometry/Point.h>
#include <Geometry/Vector.h>
#include <Geometry/AffineTransformation.h>

class LidarPlaneFitter
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Point<double,3> Point; // Type for points
	typedef Geometry::Vector<double,3> Vector; // Type for vectors
	typedef LidarPoint::Scalar LScalar;
	typedef LidarPoint::Point LPoint;
	typedef Geometry::AffineTransformation<LScalar,3> ATransform; // Type for affine transformations
	
	/* Elements: */
	private:
	ATransform planeProjection; // Affine transformation to project points into the plane's coordinate system
	LScalar min[2],max[2]; // Bounding box of all processed points in plane's coordinate system
	size_t numPoints; // Number of accumulated points
	LScalar ms; // Accumulated RMS distance from points to plane
	
	/* Constructors and destructors: */
	public:
	LidarPlaneFitter(const Point& planeCentroid,const Vector planeFrame[3])
		{
		/* Assemble the plane coordinate system: */
		Geometry::AffineTransformation<double,3> planeCoord;
		for(int i=0;i<3;++i)
			planeCoord.setDirection(i,planeFrame[i]);
		planeCoord.setOrigin(planeCentroid);
		
		/* Calculate the plane projection transformation: */
		planeProjection=ATransform(Geometry::invert(planeCoord));
		
		/* Initialize the bounding box: */
		for(int i=0;i<2;++i)
			{
			min[i]=Math::Constants<LScalar>::max;
			max[i]=Math::Constants<LScalar>::min;
			}
		numPoints=0;
		ms=LScalar(0);
		};
	
	/* Methods: */
	void operator()(const LidarPoint& lp) // Process the given LiDAR point
		{
		/* Transform the point to plane coordinates: */
		LPoint p=planeProjection.transform(lp);
		
		/* Add the point to the bounding box: */
		for(int i=0;i<2;++i)
			{
			if(min[i]>p[i])
				min[i]=p[i];
			if(max[i]<p[i])
				max[i]=p[i];
			}
		
		/* Add the point to the RMS distance: */
		++numPoints;
		ms+=Math::sqr(p[2]);
		};
	double getMin(int dimension) const
		{
		return double(min[dimension]);
		};
	double getMax(int dimension) const
		{
		return double(max[dimension]);
		};
	double getRMS(void) const
		{
		return Math::sqrt(double(ms)/double(numPoints));
		};
	};

#endif
