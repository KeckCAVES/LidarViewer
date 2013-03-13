/***********************************************************************
NormalCalculator - Functor classes to calculate a normal vector for each
point in a LiDAR data set.
Copyright (c) 2008-2010 Oliver Kreylos

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

class NormalCalculator
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Plane<double,3> Plane; // Type for planes
	private:
	typedef Geometry::ComponentArray<double,3> CA;
	typedef Geometry::Matrix<double,3,3> Matrix;
	
	/* Elements: */
	private:
	Point queryPoint; // The query point for which to calculate a plane equation
	Scalar radius2; // Squared search radius around query point
	double pxpxs,pxpys,pxpzs,pypys,pypzs,pzpzs,pxs,pys,pzs; // Accumulated components of covariance matrix
	size_t numPoints; // Number of accumulated points
	
	/* Private methods: */
	Plane::Vector calcEigenvector(const Matrix& cov,double eigenvalue) const; // Returns the eigenvector of the given matrix for the given eigenvalue
	
	/* Constructors and destructors: */
	public:
	NormalCalculator(const Point& sQueryPoint,Scalar sRadius2); // Creates an empty normal calculator
	
	/* Methods: */
	void operator()(const LidarPoint& lp) // Process the given LiDAR point
		{
		/* Check if the point is inside the search radius: */
		if(Geometry::sqrDist(lp,queryPoint)<=radius2)
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
	size_t getNumPoints(void) const // Returns the number of processed points
		{
		return numPoints;
		}
	Plane calcPlane(void) const; // Returns the least-squares plane fitting the processed points
	};

class NodeNormalCalculator
	{
	/* Elements: */
	private:
	LidarProcessOctree& lpo; // The processed LiDAR octree
	Scalar radius2; // The squared search radius around each LiDAR point
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
	NodeNormalCalculator(LidarProcessOctree& sLpo,Scalar sRadius,const char* normalFileName,unsigned int sNumThreads =1); // Creates a normal calculator with the given parameters
	~NodeNormalCalculator(void);
	
	/* Methods: */
	void saveOutlierPoints(const char* outlierFileName); // Saves outlier points to the given file
	void operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel);
	};

#endif
