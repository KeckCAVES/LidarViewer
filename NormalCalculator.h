/***********************************************************************
NormalCalculator - Functor classes to calculate a normal vector for each
point in a LiDAR data set.
Copyright (c) 2008-2014 Oliver Kreylos

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

#ifndef NORMALCALCULATOR_INCLUDED
#define NORMALCALCULATOR_INCLUDED

#include <Threads/Thread.h>
#include <Threads/Barrier.h>
#include <Geometry/ComponentArray.h>
#include <Geometry/Matrix.h>
#include <Geometry/Plane.h>

#include "LidarTypes.h"
#include "LidarFile.h"
#include "LidarProcessOctree.h"

/* Forward declarations: */
namespace Misc {
class File;
}

class NormalCalculator // Base class for functors calculating normal vectors based on point neighborhoods
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Plane<double,3> Plane; // Type for planes
	protected:
	typedef Geometry::ComponentArray<double,3> CA;
	typedef Geometry::Matrix<double,3,3> Matrix;
	
	/* Protected methods: */
	protected:
	static Plane::Vector calcEigenvector(const Matrix& cov,double eigenvalue); // Returns the eigenvector of the given matrix for the given eigenvalue
	static Plane::Vector calcNormal(const Matrix& cov); // Returns the normal vector of the plane defined by the given covariance matrix
	};

class RadiusNormalCalculator:public NormalCalculator // Class to calculate normal vectors by accumulating points from a neighborhood of fixed size
	{
	/* Elements: */
	private:
	Scalar radius2; // Squared search radius around query point
	Point queryPoint; // The query point for which to calculate a plane equation
	double pxpxs,pxpys,pxpzs,pypys,pypzs,pzpzs,pxs,pys,pzs; // Accumulated components of covariance matrix
	size_t numPoints; // Number of accumulated points
	Scalar closestDist2; // Squared distance to query point's closest non-identical neighbor
	
	/* Constructors and destructors: */
	public:
	RadiusNormalCalculator(Scalar sRadius); // Creates normal calculator for the given search radius
	
	/* Methods: */
	void operator()(const LidarPoint& lp) // Process the given LiDAR point
		{
		/* Check if the point is inside the search radius: */
		Scalar dist2=Geometry::sqrDist(lp,queryPoint);
		if(dist2<=radius2)
			{
			/* Accumulate the node point: */
			pxpxs+=double(lp[0])*double(lp[0]);
			pxpys+=double(lp[0])*double(lp[1]);
			pxpzs+=double(lp[0])*double(lp[2]);
			pypys+=double(lp[1])*double(lp[1]);
			pypzs+=double(lp[1])*double(lp[2]);
			pzpzs+=double(lp[2])*double(lp[2]);
			pxs+=double(lp[0]);
			pys+=double(lp[1]);
			pzs+=double(lp[2]);
			
			++numPoints;
			
			/* Check if this is the closest neighbor yet: */
			if(dist2>Scalar(0)&&closestDist2>=dist2)
				closestDist2=dist2;
			}
		}
	const Point& getQueryPoint(void) const
		{
		return queryPoint;
		}
	Scalar getQueryRadius2(void) const
		{
		return radius2;
		}
	
	void prepare(const Point& newQueryPoint); // Prepares the normal calculator for traversal around the given point
	size_t getNumPoints(void) const // Returns the number of processed points
		{
		return numPoints;
		}
	Plane calcPlane(void) const; // Returns the least-squares plane fitting the processed points
	Scalar getClosestDist(void) const // Returns the distance to the closest non-identical neighbor
		{
		return Math::sqrt(closestDist2);
		}
	};

class NumberRadiusNormalCalculator:public NormalCalculator // Class to calculate normal vectors by accumulating points from a neighborhood of maximum size and maximum number of neighbors
	{
	/* Embedded classes: */
	private:
	struct Neighbor // Structure to store a neighbor
		{
		/* Elements: */
		public:
		Point point; // Copy of neighbor's position
		Scalar dist2; // Squared distance from query position to neighbor
		};
	
	/* Elements: */
	private:
	unsigned int maxNumNeighbors; // Maximum number of neighbors
	Scalar maxDist2; // Squared maximum distance to neighbors
	Point queryPoint; // The query point position
	Neighbor* neighbors; // Array of current neighbor candidates, organized as a heap
	unsigned int currentNumNeighbors; // Current number of neighbor candidates
	Scalar currentMaxDist2; // Current maximum distance to any neighbor candidate
	
	/* Constructors and destructors: */
	public:
	NumberRadiusNormalCalculator(unsigned int sMaxNumNeighbors); // Creates normal calculator for the given number of neighbors and unlimited neighborhood size
	NumberRadiusNormalCalculator(unsigned int sMaxNumNeighbors,Scalar sMaxDist); // Creates normal calculator for the given number of neighbors and given neighborhood size
	NumberRadiusNormalCalculator(const NumberRadiusNormalCalculator& source); // Copy constructor
	NumberRadiusNormalCalculator& operator=(const NumberRadiusNormalCalculator& source); // Assignment operator
	~NumberRadiusNormalCalculator(void);
	
	/* Methods: */
	void operator()(const LidarPoint& point);
	const Point& getQueryPoint(void) const
		{
		return queryPoint;
		}
	Scalar getQueryRadius2(void) const
		{
		return currentMaxDist2;
		}
	
	void prepare(const Point& newQueryPoint); // Prepares the normal calculator for traversal around the given point
	unsigned int getNumPoints(void) const // Returns the number of neighbors in the neighborhood
		{
		return currentNumNeighbors;
		}
	Plane calcPlane(void) const; // Returns the least-squares plane fitting the processed points
	Scalar getClosestDist(void) const; // Returns the distance to the closest non-identical neighbor
	};

template <class NormalCalculatorParam>
class NodeNormalCalculator
	{
	/* Embedded classes: */
	public:
	typedef NormalCalculatorParam NormalCalculator; // Normal calculator class
	
	/* Elements: */
	private:
	LidarProcessOctree& lpo; // The processed LiDAR octree
	const NormalCalculator& normalCalculator; // The "prototype" normal calculator
	Vector* normalBuffer; // Array to hold normal vectors for a node during processing
	Vector* childNormalBuffers[8]; // Array of normal arrays for a node's children during subsampling
	LidarFile::Offset normalDataSize; // Size of each record in the normal file
	LidarFile normalFile; // The file to which to write the normal vector data
	bool saveOutliers; // Flag whether to save points with undefined normal vectors to an outlier file
	Threads::Mutex outlierMutex; // Mutex serializing access to the outlier array
	size_t numOutliers; // Number of outliers encountered in the current node
	Point* outliers; // Array of outlier positions
	Misc::File* outlierFile; // If outliers are to be saved, pointer to the file to which to write them
	unsigned int numThreads; // Number of processing threads
	bool shutdownThreads; // Flag to shut down threads at the end
	Threads::Thread* calcThreads; // Array of threads to calculate normal vectors from point neighborhoods
	Threads::Barrier calcBarrier;
	Threads::Thread* subsampleThreads; // Array of threads to subsample normal vectors from node children
	Threads::Barrier subsampleBarrier;
	LidarProcessOctree::Node* currentNode; // Pointer to currently processed node
	LidarProcessOctree::Node* currentChildren[8]; // Pointer to children of currently processed node
	size_t numProcessedNodes; // Number of already processed nodes
	
	/* Private methods: */
	void* calcThreadMethod(unsigned int threadIndex);
	void* subsampleThreadMethod(unsigned int threadIndex);
	
	/* Constructors and destructors: */
	public:
	NodeNormalCalculator(LidarProcessOctree& sLpo,const NormalCalculator& sNormalCalculator,const char* normalFileName,unsigned int sNumThreads =1); // Creates a node normal calculator with the given point normal calculator and parameters
	~NodeNormalCalculator(void);
	
	/* Methods: */
	void saveOutlierPoints(const char* outlierFileName); // Saves outlier points to the given file
	void operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel);
	};

#ifndef NORMALCALCULATOR_IMPLEMENTATION
#include "NormalCalculator.icpp"
#endif

#endif
