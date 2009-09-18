/***********************************************************************
PointClassifier - Point processor functor class to classify selected
points into surface, edge, and corner classes using Patric Keller's PCA-
based point classification algorithm.
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

#ifndef POINTCLASSIFIER_INCLUDED
#define POINTCLASSIFIER_INCLUDED

#include <Math/Math.h>
#include <Geometry/PCACalculator.h>
#include <GL/gl.h>
#include <GL/GLColor.h>

#include "LidarTypes.h"
#include "LidarOctree.h"
#include "PointPCACalculator.h"

class PointClassifier
	{
	/* Elements: */
	private:
	const LidarOctree* lidarOctree; // Pointer to the LiDAR octree to traverse
	Scalar radius2; // Squared radius of the neighborhood for each point
	
	/* Constructors and destructors: */
	public:
	PointClassifier(const LidarOctree* sLidarOctree,Scalar sRadius) // Creates a point classifier for the given LiDAR octree with the given neighborhood radius
		:lidarOctree(sLidarOctree),
		 radius2(Math::sqr(sRadius))
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& lp,LidarOctree::Color& color) // Colors the given LiDAR point based on its classification
		{
		/* Calculate the PCA of the point's neighborhood: */
		PointPCACalculator ppca(lp,radius2);
		lidarOctree->processPointsDirected(ppca);
		double eigenvalues[3];
		ppca.getEigenvalues(eigenvalues);
		
		/* Color the point by its similarity to surface-like, curve-like, and point-like features: */
		double surfacity=(eigenvalues[0]-eigenvalues[2])/eigenvalues[0]; // High for surfaces and curves
		double curvacity=(eigenvalues[0]-eigenvalues[1])/eigenvalues[0]; // High for curves
		double pointicity=(eigenvalues[1]/eigenvalues[0])*(eigenvalues[2]/eigenvalues[0]); // High for points
		
		color=GLColor<GLfloat,3>(surfacity,curvacity,pointicity);
		};
	};

#endif
