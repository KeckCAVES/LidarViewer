/***********************************************************************
LidarPlaneExtractor - Point processor functor class to extract least-
squares planes from sets of selected LiDAR points.
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

#ifndef LIDARPLANEEXTRACTOR_INCLUDED
#define LIDARPLANEEXTRACTOR_INCLUDED

#include <Geometry/Point.h>
#include <Geometry/Vector.h>
#include <Geometry/Box.h>
#include <Geometry/PCACalculator.h>

class LidarPlaneExtractor
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Point<double,3> Point; // Type for points
	typedef Geometry::Vector<double,3> Vector; // Type for vectors
	typedef Geometry::Box<double,3> Box; // Type for bounding boxes
	
	/* Elements: */
	private:
	Box bb; // Bounding box of all processed points
	Geometry::PCACalculator<3> pca; // Helper object to accumulate the points' covariance matrix and calculate their PCA
	
	/* Constructors and destructors: */
	public:
	LidarPlaneExtractor(void)
		:bb(Box::empty)
		{
		};
	
	/* Methods: */
	void operator()(const LidarPoint& lp) // Process the given LiDAR point
		{
		/* Add the node point to the bounding box: */
		bb.addPoint(lp);
		
		/* Add the point to the PCA calculator: */
		pca.accumulatePoint(lp);
		};
	size_t getNumPoints(void) const // Returns the number of processed points
		{
		return pca.getNumPoints();
		}
	const Box& getBB(void) const // Returns the processed points' bounding box
		{
		return bb;
		};
	void calcPlane(Point& centroid,Vector plane[3],double lengths[3]) // Returns the least-squares plane and its aligned normalized coordinate frame and the lengths of the eigenvectors
		{
		/* Calculate the point set's covariance matrix: */
		pca.calcCovariance();
		
		/* Calculate the covariance matrix' eigenvalues: */
		pca.calcEigenvalues(lengths);
		
		/* Calculate all eigenvectors: */
		for(int i=0;i<3;++i)
			plane[i]=pca.calcEigenvector(lengths[i]);
		
		/* Calculate the processed points' centroid: */
		centroid=pca.calcCentroid();
		};
	};

#endif
