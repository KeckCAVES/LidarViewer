/***********************************************************************
RidgeFinder - Functor classes to find ridges in LiDAR data sets.
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

#ifndef RIDGEFINDER_INCLUDED
#define RIDGEFINDER_INCLUDED

#include <Geometry/AffineTransformation.h>
#include <Geometry/PCACalculator.h>

#include "LidarTypes.h"
#include "LidarOctree.h"
#include "PointPCACalculator.h"

class RidgeFinder
	{
	/* Embedded classes: */
	private:
	class PointClassifier // Class to classify a point as ridge or else
		{
		/* Embedded classes: */
		public:
		typedef Geometry::AffineTransformation<double,3> ATransform; // Type for transformations into a plane's local coordinate system
		
		/* Elements: */
		private:
		Point queryPoint; // The query point for which to calculate a plane equation
		Scalar radius,radius2; // (Squared) search radius around query point
		ATransform toPlane; // Transformation to plane's local coordinates
		Geometry::PCACalculator<2> pcas[2]; // PCA calculators for points below and above the plane, respectively
		
		/* Constructors and destructors: */
		public:
		PointClassifier(const Point& sQueryPoint,Scalar sRadius2,const PointPCACalculator::Plane& sPlane)
			:queryPoint(sQueryPoint),radius2(sRadius2)
			{
			/* Calculate the transformation into the plane's local coordinate system: */
			Vector z=sPlane.getNormal();
			z.normalize();
			Vector x=Geometry::normal(z);
			x.normalize();
			Vector y=Geometry::cross(z,x);
			z.normalize();
			Point o=sPlane.project(queryPoint);
			ATransform::Matrix m;
			for(int i=0;i<3;++i)
				{
				m(i,0)=x[i];
				m(i,1)=y[i];
				m(i,2)=z[i];
				m(i,3)=o[i];
				}
			toPlane=ATransform(m);
			}
		
		/* Methods: */
		void operator()(const LidarPoint& lp) // Process the given LiDAR point
			{
			/* Convert the point to the plane's local coordinate system: */
			ATransform::Point pp=toPlane.transform(lp);
			
			/* Accumulate the point into the lower or upper PCA calculator: */
			if(pp[2]<0.0)
				pcas[0].accumulatePoint(pp);
			else
				pcas[1].accumulatePoint(pp);
			}
		const Point& getQueryPoint(void) const
			{
			return queryPoint;
			}
		Scalar getQueryRadius2(void) const
			{
			return radius2;
			}
		bool isRidge(void) // Returns true if the query point is a ridge point
			{
			if(pcas[0].getNumPoints()>=2&&pcas[1].getNumPoints()>=2)
				{
				double radius=Math::sqrt(radius2);
				
				/* Calculate eigenvalues of both PCAs: */
				double evs[2][2];
				unsigned int numEvs[2];
				bool isRidges[2];
				for(int i=0;i<2;++i)
					{
					pcas[i].calcCovariance();
					numEvs[i]=pcas[i].calcEigenvalues(evs[i]);
					isRidges[i]=numEvs[i]==2&&Math::abs(evs[i][0])>=radius*0.75&&Math::abs(evs[i][1])<=radius*0.333;
					}
				
				if(isRidges[0]&&!isRidges[1])
					{
					ATransform::Point qp=toPlane.transform(queryPoint);
					return pcas[0].calcEigenvector(evs[0][1])*(Geometry::PCACalculator<2>::Point(qp[0],qp[1])-pcas[0].calcCentroid())<=1.0;
					}
				else if(!isRidges[0]&&isRidges[1])
					{
					ATransform::Point qp=toPlane.transform(queryPoint);
					return pcas[1].calcEigenvector(evs[1][1])*(Geometry::PCACalculator<2>::Point(qp[0],qp[1])-pcas[1].calcCentroid())<=1.0;
					}
				else
					return false;
				}
			else
				return false;
			}
		};
	
	/* Elements: */
	private:
	const LidarOctree* lidarOctree; // Pointer to the LiDAR octree to traverse
	Scalar radius2; // Squared radius of the neighborhood for each point
	
	/* Constructors and destructors: */
	public:
	RidgeFinder(const LidarOctree* sLidarOctree,Scalar sRadius) // Creates a ridge finder for the given LiDAR octree with the given neighborhood radius
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
		
		if(ppca.getNumPoints()>=3)
			{
			/* Run a ridge classifier on the point: */
			PointClassifier pc(lp,radius2,ppca.getPlane());
			lidarOctree->processPointsDirected(pc);
			
			/* Color the point: */
			if(pc.isRidge())
				color=GLColor<GLfloat,3>(1.0,1.0,0.0);
			else
				color=GLColor<GLfloat,3>(0.0,0.0,1.0);
			}
		else
			color=GLColor<GLfloat,3>(0.0,0.0,1.0);
		};
	};

#if 0

#include <Geometry/Plane.h>
#include <Geometry/PCACalculator.h>

class PlaneCalculator // Class to find the best-fitting plane for a point neighborhood
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Plane<double,3> Plane; // Type for planes
	
	/* Elements: */
	private:
	Point queryPoint; // The query point for which to calculate a plane equation
	Scalar radius2; // Squared search radius around query point
	Geometry::PCACalculator<3> pca; // Helper structure to calculate principal components of traversed point set
	
	/* Constructors and destructors: */
	public:
	PlaneCalculator(const Point& sQueryPoint,Scalar sRadius2) // Creates an empty plane calculator
		:queryPoint(sQueryPoint),radius2(sRadius2)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& lp) // Process the given LiDAR point
		{
		/* Accumulate the point into the PCA calculator: */
		pca.accumulatePoint(lp);
		}
	const Point& getQueryPoint(void) const
		{
		return queryPoint;
		}
	Scalar getQueryRadius2(void) const
		{
		return radius2;
		}
	size_t getNumPoints(void) const // Returns the number of processed points
		{
		return pca.getNumPoints();
		}
	Plane calcPlane(void) const; // Returns the least-squares plane fitting the processed points
	};

class PointClassifier // Class to classify a point as ridge or else
	{
	/* Embedded classes: */
	public:
	typedef Geometry::AffineTransformation<double,3> ATransform; // Type for transformations into a plane's local coordinate system
	
	/* Elements: */
	private:
	Point queryPoint; // The query point for which to calculate a plane equation
	Scalar radius,radius2; // (Squared) search radius around query point
	PlaneCalculator::Plane plane; // The neighborhood's best-fitting plane
	ATransform toPlane; // Transformation to plane's local coordinates
	PCACalculator<2> pcas[2]; // PCA calculators for points below and above the plane, respectively
	
	/* Constructors and destructors: */
	PointClassifier(const Point& sQueryPoint,Scalar sRadius2,const PlaneCalculator::Plane& sPlane);
	
	/* Methods: */
	void operator()(const LidarPoint& lp) // Process the given LiDAR point
		{
		/* Convert the point to the plane's local coordinate system: */
		ATransform::Point pp=toPlane.transform(lp);
		
		/* Accumulate the point into the lower or upper PCA calculator: */
		if(pp[2]<0.0)
			pcas[0].accumulatePoint(pp);
		else
			pcas[1].accumulatePoint(pp);
		}
	const Point& getQueryPoint(void) const
		{
		return queryPoint;
		}
	Scalar getQueryRadius2(void) const
		{
		return radius2;
		}
	bool isRidge(void) const; // Returns true if the query point is a ridge point
	};

class RidgeFinder // Class to classify points in a LiDAR data set as ridges or else
	{
	/* Elements: */
	private:
	LidarProcessOctree& lpo; // The processed LiDAR octree
	Scalar radius2; // The squared search radius around each LiDAR point
	Color* colorBuffer; // Buffer to collect a node's point colors
	LidarFile::Offset colorDataSize; // Size of each record in the color file
	LidarFile colorFile; // The file to which to write the classification color data
	size_t numProcessedNodes; // Number of already processed nodes
	
	/* Constructors and destructors: */
	public:
	RidgeFinder(LidarProcessOctree& sLpo,Scalar sRadius,const char* colorFileName); // Creates a ridge finder with the given parameters
	~RidgeFinder(void);
	
	/* Methods: */
	void operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel);
	};

#endif

#endif
