/***********************************************************************
PointPCACalculator - Point processor functor class to calculate the
three-dimensional PCA of a point's neighborhood.
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

#ifndef POINTPCACALCULATOR_INCLUDED
#define POINTPCACALCULATOR_INCLUDED

#include <Geometry/Plane.h>
#include <Geometry/PCACalculator.h>

#include "LidarTypes.h"

class PointPCACalculator
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Plane<double,3> Plane; // Type for planes
	
	/* Elements: */
	private:
	Point point; // The query point
	Scalar radius2; // The squared radius of the neighborhood
	Geometry::PCACalculator<3> pca; // PCA calculation object
	
	/* Constructors and destructors: */
	public:
	PointPCACalculator(const Point& sPoint,Scalar sRadius2)
		:point(sPoint),radius2(sRadius2)
		{
		}
	
	/* Methods: */
	const Point& getQueryPoint(void) const
		{
		return point;
		}
	Scalar getQueryRadius2(void) const
		{
		return radius2;
		}
	void operator()(const LidarPoint& lp)
		{
		if(Geometry::sqrDist(point,lp)<=radius2)
			pca.accumulatePoint(lp);
		}
	size_t getNumPoints(void) const
		{
		return pca.getNumPoints();
		}
	unsigned int getEigenvalues(double eigenvalues[3])
		{
		/* Calculate the accumulated point set's eigenvalues: */
		pca.calcCovariance();
		return pca.calcEigenvalues(eigenvalues);
		}
	Plane getPlane(void) // Returns the best-fitting plane for the point neighborhood
		{
		/* Calculate the accumulated point set's eigenvalues: */
		pca.calcCovariance();
		double evs[3];
		pca.calcEigenvalues(evs);
		
		/* Calculate the eigenvector of the smallest eigenvalue: */
		Plane::Vector normal=pca.calcEigenvector(evs[2]);
		
		/* Calculate the centroid of the point neighborhood: */
		Plane::Point centroid=pca.calcCentroid();
		
		/* Return the resulting plane equation: */
		return Plane(normal,centroid);
		}
	};

#endif
