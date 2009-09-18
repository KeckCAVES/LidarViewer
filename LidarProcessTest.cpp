/***********************************************************************
LidarProcessTest - Test program for the LiDAR processing octree data
structure.
Copyright (c) 2008-2009 Oliver Kreylos

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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <Misc/Timer.h>
#include <Threads/Thread.h>
#include <Math/Constants.h>
#include <Math/Random.h>
#include <Geometry/ArrayKdTree.h>

#include "LidarTypes.h"
#include "LidarProcessOctree.h"

class PointCounter // Functor class to count the number of points stored in an octree file
	{
	/* Elements: */
	private:
	size_t numPoints;
	
	/* Constructors and destructors: */
	public:
	PointCounter(void)
		:numPoints(0)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& p)
		{
		++numPoints;
		}
	size_t getNumPoints(void) const
		{
		return numPoints;
		}
	};

class DirectedProcessFunctor
	{
	/* Elements: */
	protected:
	Point queryPoint;
	Scalar queryRadius2;
	
	/* Constructors and destructors: */
	public:
	DirectedProcessFunctor(const Point& sQueryPoint,Scalar sQueryRadius2)
		:queryPoint(sQueryPoint),queryRadius2(sQueryRadius2)
		{
		}
	
	/* Methods: */
	const Point& getQueryPoint(void) const
		{
		return queryPoint;
		}
	Scalar getQueryRadius2(void) const
		{
		return queryRadius2;
		}
	};

class NeighborCounter:public DirectedProcessFunctor
	{
	public:
	size_t numPoints;
	
	NeighborCounter(const Point& queryPoint,Scalar queryRadius2)
		:DirectedProcessFunctor(queryPoint,queryRadius2),
		 numPoints(0)
		{
		}
	
	void operator()(const LidarPoint& point)
		{
		++numPoints;
		}
	};

class PointDensityCalculator
	{
	struct ThreadArgs
		{
		unsigned int numPoints;
		const LidarPoint* points;
		size_t numNeighbors;
		};
	
	LidarProcessOctree& lpo;
	Scalar neighborhoodRadius;
	int numThreads;
	Threads::Thread* threads;
	ThreadArgs* threadArgs;
	public:
	size_t numPoints;
	size_t totalNumNeighbors;
	
	void* processThreadMethod(ThreadArgs* args)
		{
		args->numNeighbors=0;
		for(unsigned int i=0;i<args->numPoints;++i)
			{
			NeighborCounter nc(args->points[i],neighborhoodRadius*neighborhoodRadius);
			lpo.processPointsDirected(nc);
			args->numNeighbors+=nc.numPoints;
			}
		}
	
	PointDensityCalculator(LidarProcessOctree& sLpo,Scalar sNeighborhoodRadius,int sNumThreads)
		:lpo(sLpo),
		 neighborhoodRadius(sNeighborhoodRadius),
		 numThreads(sNumThreads),threads(new Threads::Thread[numThreads]),threadArgs(new ThreadArgs[numThreads]),
		 numPoints(0),totalNumNeighbors(0)
		{
		}
	~PointDensityCalculator(void)
		{
		delete[] threads;
		delete[] threadArgs;
		}
	
	void operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel)
		{
		if(!node.isLeaf()||node.getNumPoints()==0)
			return;
		
		for(int i=0;i<numThreads;++i)
			{
			unsigned int firstPoint=(node.getNumPoints()*i)/numThreads;
			unsigned int lastPoint=(node.getNumPoints()*(i+1))/numThreads;
			threadArgs[i].numPoints=lastPoint-firstPoint;
			threadArgs[i].points=node.getPoints()+firstPoint;
			threads[i].start(this,&PointDensityCalculator::processThreadMethod,&threadArgs[i]);
			}
		for(int i=0;i<numThreads;++i)
			{
			threads[i].join();
			totalNumNeighbors+=threadArgs[i].numNeighbors;
			}
		numPoints+=node.getNumPoints();
		}
	void operator()(const LidarPoint& point)
		{
		++numPoints;
		NeighborCounter nc(point,neighborhoodRadius*neighborhoodRadius);
		lpo.processPointsDirected(nc);
		totalNumNeighbors+=nc.numPoints;
		}
	};

int main(int argc,char* argv[])
	{
	LidarProcessOctree lpo(argv[1],512*1024*1024);
	
	PointCounter pc;
	Scalar radius=Scalar(0.1);
	PointDensityCalculator pdc(lpo,radius,4);
	#if 0
	Box box;
	box.min=Point(930,530,0);
	box.max=Point(970,570,100);
	lpo.processPointsInBox(box,pdc);
	#else
	// lpo.processPoints(pdc);
	lpo.processNodesPostfix(pdc);
	#endif
	
	std::cout<<"Number of processed points: "<<pdc.numPoints<<std::endl;
	std::cout<<"Total number of found neighbors: "<<pdc.totalNumNeighbors<<std::endl;
	std::cout<<"Total number of loaded octree nodes: "<<lpo.getNumSubdivideCalls()<<", "<<lpo.getNumLoadedNodes()<<std::endl;
	std::cout<<"Average point density in 1/m^3: "<<pdc.totalNumNeighbors/(pdc.numPoints*4.0/3.0*3.141592654*radius*radius*radius)<<std::endl;
	
	return 0;
	}

#if 0

class NearestNeighborFinder // Functor class to find the nearest non-duplicate neighbor of a point in an octree file
	{
	/* Elements: */
	private:
	Point queryPoint; // The query point position
	Scalar dist2; // Squared distance to current nearest neighbor candidate
	const LidarPoint* nearest; // Pointer to nearest neighbor
	
	/* Constructors and destructors: */
	public:
	NearestNeighborFinder(const Point& sQueryPoint)
		:queryPoint(sQueryPoint),
		 dist2(Math::Constants<Scalar>::max),
		 nearest(0)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& point)
		{
		Scalar pDist2=Geometry::sqrDist(point,queryPoint);
		if(pDist2>Scalar(0)&&pDist2<dist2)
			{
			nearest=&point;
			dist2=pDist2;
			}
		}
	const Point& getQueryPoint(void) const
		{
		return queryPoint;
		}
	Scalar getQueryRadius2(void) const
		{
		return dist2;
		}
	const LidarPoint* getNearest(void) const
		{
		return nearest;
		}
	Scalar getNearestDistance2(void) const
		{
		return dist2;
		}
	};

class KNearestNeighborFinder // Functor class to find the k nearest neighbors of a point in an octree file
	{
	/* Embedded classes: */
	private:
	struct Neighbor // Structure to store a neighbor
		{
		/* Elements: */
		public:
		const LidarPoint* point; // Pointer to neighbor
		Scalar dist2; // Squared distance from query position to neighbor
		};
	
	/* Elements: */
	private:
	Point queryPoint; // The query point position
	int maxNumNeighbors; // Maximum number of neighbors to find
	Neighbor* neighbors; // Array of current neighbor candidates
	int numNeighbors; // Current number of neighbor candidates
	Scalar maxDist2; // Current maximum distance to any neighbor candidate
	
	/* Constructors and destructors: */
	public:
	KNearestNeighborFinder(const Point& sQueryPoint,int sMaxNumNeighbors)
		:queryPoint(sQueryPoint),
		 maxNumNeighbors(sMaxNumNeighbors),
		 neighbors(new Neighbor[maxNumNeighbors]),
		 numNeighbors(0),
		 maxDist2(Math::Constants<Scalar>::max)
		{
		}
	~KNearestNeighborFinder(void)
		{
		delete[] neighbors;
		}
	
	/* Methods: */
	void operator()(const LidarPoint& point)
		{
		Scalar dist2=Geometry::sqrDist(point,queryPoint);
		if(numNeighbors<maxNumNeighbors)
			{
			/* Insert the new point into the heap: */
			int insertionPos=numNeighbors;
			while(insertionPos>0)
				{
				int parent=(insertionPos-1)>>1;
				if(neighbors[parent].dist2>=dist2)
					break;
				neighbors[insertionPos]=neighbors[parent];
				insertionPos=parent;
				}
			neighbors[insertionPos].point=&point;
			neighbors[insertionPos].dist2=dist2;
			
			++numNeighbors;
			if(numNeighbors==maxNumNeighbors)
				maxDist2=neighbors[0].dist2;
			}
		else if(dist2<maxDist2)
			{
			/* Replace the currently farthest-away neighbor in the heap: */
			int insertionPos=0;
			while(true)
				{
				int biggestIndex=insertionPos;
				Scalar biggest=dist2;
				int child=(insertionPos<<1);
				for(int i=0;i<2;++i)
					{
					++child;
					if(child<maxNumNeighbors&&neighbors[child].dist2>biggest)
						{
						biggestIndex=child;
						biggest=neighbors[child].dist2;
						}
					}
				if(biggestIndex==insertionPos)
					break;
				neighbors[insertionPos]=neighbors[biggestIndex];
				insertionPos=biggestIndex;
				}
			neighbors[insertionPos].point=&point;
			neighbors[insertionPos].dist2=dist2;
			
			maxDist2=neighbors[0].dist2;
			}
		}
	const Point& getQueryPoint(void) const
		{
		return queryPoint;
		}
	Scalar getQueryRadius2(void) const
		{
		return maxDist2;
		}
	int getNumNeighbors(void) const
		{
		return numNeighbors;
		}
	const LidarPoint* getNeighbor(int index) const
		{
		return neighbors[index].point;
		}
	Scalar getNeighborDistance2(int index) const
		{
		return neighbors[index].dist2;
		}
	};

class PointNearestNeighborFinder // Functor class to find the nearest neighbors of all points in an octree file
	{
	/* Elements: */
	private:
	LidarProcessOctree& lpo;
	size_t numPoints;
	double totalDistance2;
	
	/* Constructors and destructors: */
	public:
	PointNearestNeighborFinder(LidarProcessOctree& sLpo)
		:lpo(sLpo),
		 numPoints(0),totalDistance2(0.0)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& p)
		{
		/* Find the nearest neighbor of this point: */
		NearestNeighborFinder nnf(p);
		lpo.processPointsDirected(nnf);
		++numPoints;
		totalDistance2+=double(nnf.getNearestDistance2());
		}
	size_t getNumPoints(void) const
		{
		return numPoints;
		}
	Scalar getAverageDist(void) const
		{
		return Scalar(Math::sqrt(totalDistance2/double(numPoints)));
		}
	};

class PointKNearestNeighborFinder // Functor class to find the k nearest neighbors of all points in an octree file
	{
	/* Elements: */
	private:
	LidarProcessOctree& lpo;
	int maxNumNeighbors;
	size_t numPoints;
	double totalDistance2;
	
	/* Constructors and destructors: */
	public:
	PointKNearestNeighborFinder(LidarProcessOctree& sLpo,int sMaxNumNeighbors)
		:lpo(sLpo),
		 maxNumNeighbors(sMaxNumNeighbors),
		 numPoints(0),totalDistance2(0.0)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& p)
		{
		/* Find the k nearest neighbors of this point: */
		KNearestNeighborFinder knnf(p,maxNumNeighbors);
		lpo.processPointsDirected(knnf);
		++numPoints;
		totalDistance2+=double(knnf.getNeighborDistance2(0));
		}
	size_t getNumPoints(void) const
		{
		return numPoints;
		}
	Scalar getAverageDist(void) const
		{
		return Scalar(Math::sqrt(totalDistance2/double(numPoints)));
		}
	};

class PointExtractor // Functor class to load all points from an octree file into a std::vector
	{
	/* Elements: */
	private:
	std::vector<LidarPoint> points; // Point array
	
	/* Constructors and destructors: */
	public:
	PointExtractor(void)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& p)
		{
		points.push_back(p);
		}
	const std::vector<LidarPoint>& getPoints(void) const
		{
		return points;
		}
	};

class NearestNeighborFinderKdtree // Functor class to find the non-duplicate nearest neighbor of a point in an in-memory kd-tree
	{
	/* Elements: */
	private:
	Point queryPoint; // The query point position
	Scalar dist2; // Squared distance to current nearest neighbor candidate
	const LidarPoint* nearest; // Pointer to nearest neighbor
	
	/* Constructors and destructors: */
	public:
	NearestNeighborFinderKdtree(const Point& sQueryPoint)
		:queryPoint(sQueryPoint),
		dist2(Math::Constants<Scalar>::max),
		nearest(0)
		{
		}
	
	/* Methods: */
	bool operator()(const LidarPoint& point,int splitDimension)
		{
		Scalar pDist2=Geometry::sqrDist(point,queryPoint);
		if(pDist2>Scalar(0)&&pDist2<dist2)
			{
			nearest=&point;
			dist2=pDist2;
			}
		
		/* Stop traversal if split plane is farther away than closest point: */
		return dist2>Math::sqr(point[splitDimension]-queryPoint[splitDimension]);
		}
	const Point& getQueryPosition(void) const
		{
		return queryPoint;
		}
	Scalar getQueryRadius2(void) const
		{
		return dist2;
		}
	const LidarPoint* getNearest(void) const
		{
		return nearest;
		}
	Scalar getNearestDistance2(void) const
		{
		return dist2;
		}
	};

int main(int argc,char* argv[])
	{
	const char* fileName=0;
	int maxNumNeighbors=10;
	int cacheSize=512;
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"numNeighbors")==0)
				{
				++i;
				maxNumNeighbors=atoi(argv[i]);
				}
			else if(strcasecmp(argv[i]+1,"cache")==0)
				{
				++i;
				cacheSize=atoi(argv[i]);
				}
			}
		else if(fileName==0)
			fileName=argv[i];
		}
	if(fileName==0)
		{
		std::cerr<<"No file name provided"<<std::endl;
		return 1;
		}
	
	Misc::Timer t;
	
	/* Create a processing octree: */
	LidarProcessOctree lpo(fileName,size_t(cacheSize)*size_t(1024*1024));
	
	#if 1
	
	/* Count the number of points in the octree: */
	PointCounter pc;
	const Cube& c=lpo.getDomain();
	// Box box(c.getMin(),c.getMax());
	Box box(Point(700.0,300.0,-200.0),Point(1200.0,800.0,300.0));
	lpo.processPointsInBox(box,pc);
	std::cout<<"Octree contains "<<pc.getNumPoints()<<" points"<<std::endl;
	
	#elif 0
	
	if(maxNumNeighbors==1)
		{
		/* Find the nearest neighbor of each point in the octree: */
		PointNearestNeighborFinder pnnf(lpo);
		lpo.processPoints(pnnf);
		std::cout<<"Octree contains "<<pnnf.getNumPoints()<<" points, average point density is "<<pnnf.getAverageDist()<<std::endl;
		}
	else
		{
		/* Find the k nearest neighbors of each point in the octree: */
		PointKNearestNeighborFinder pknnf(lpo,maxNumNeighbors);
		lpo.processPoints(pknnf);
		std::cout<<"Octree contains "<<pknnf.getNumPoints()<<" points, average point density is "<<pknnf.getAverageDist()<<std::endl;
		}
	
	#elif 0
	
	/* Extract the points from the octree: */
	PointExtractor pe;
	lpo.processPoints(pe);
	
	/* Build a kd-tree from the points: */
	Geometry::ArrayKdTree<LidarPoint> kdTree;
	LidarPoint* points=kdTree.createTree(pe.getPoints().size());
	for(size_t i=0;i<pe.getPoints().size();++i)
		points[i]=pe.getPoints()[i];
	kdTree.releasePoints();
	
	#if 0
	/* Find the nearest neighbor of all points in the kd-tree: */
	double nearestDistance2=0.0;
	for(size_t i=0;i<pe.getPoints().size();++i)
		{
		NearestNeighborFinderKdtree nnfkd(points[i]);
		kdTree.traverseTreeDirected(nnfkd);
		nearestDistance2+=double(nnfkd.getNearestDistance2());
		}
	std::cout<<"Octree contains "<<pe.getPoints().size()<<" points, average point density is "<<Math::sqrt(nearestDistance2/double(pe.getPoints().size()))<<std::endl;
	#else
	/* Search points randomly to test for correctness: */
	for(int index=0;index<pe.getPoints().size();++index)
		{
		/* Use the kd-tree: */
		NearestNeighborFinderKdtree nnfkd(points[index]);
		kdTree.traverseTreeDirected(nnfkd);
		
		/* Use the out-of-core search algorithm: */
		NearestNeighborFinder nnf(points[index]);
		lpo.processPointsDirected(nnf);
		
		if(nnfkd.getNearestDistance2()!=nnf.getNearestDistance2())
			{
			/* Whoopsie! */
			std::cout<<"Mismatch for search point "<<points[index][0]<<", "<<points[index][1]<<", "<<points[index][2]<<":"<<std::endl;
			std::cout<<"Kd-tree: "<<nnfkd.getNearestDistance2()<<" "<<(*nnfkd.getNearest())[0]<<", "<<(*nnfkd.getNearest())[1]<<", "<<(*nnfkd.getNearest())[0]<<std::endl;
			std::cout<<"Octree : "<<nnf.getNearestDistance2()<<" "<<(*nnf.getNearest())[0]<<", "<<(*nnf.getNearest())[1]<<", "<<(*nnf.getNearest())[0]<<std::endl;
			}
		}
	#endif
	
	#else
	
	/* Extract the points from the octree: */
	PointExtractor pe;
	lpo.processPoints(pe);
	
	/* Compare a brute-force search against the octree-based algorithm: */
	const std::vector<LidarPoint>& points=pe.getPoints();
	for(int i=0;i<100000;++i)
		{
		size_t index=Math::randUniformCO(0,points.size());
		
		/* Use the out-of-core search algorithm: */
		NearestNeighborFinder nnf(points[index]);
		lpo.processPointsDirected(nnf);
		
		/* Use brute force: */
		const LidarPoint* nearest=0;
		Scalar minDist2=Math::Constants<Scalar>::max;
		for(size_t i=0;i<points.size();++i)
			{
			Scalar dist2=Geometry::sqrDist(points[index],points[i]);
			if(dist2>Scalar(0)&&dist2<minDist2)
				{
				nearest=&points[i];
				minDist2=dist2;
				}
			}
		
		if(nnf.getNearestDistance2()!=minDist2)
			{
			/* Whoopsie! */
			std::cout<<"Mismatch for search point "<<points[index][0]<<", "<<points[index][1]<<", "<<points[index][2]<<":"<<std::endl;
			std::cout<<"Octree     : "<<nnf.getNearestDistance2()<<" "<<(*nnf.getNearest())[0]<<", "<<(*nnf.getNearest())[1]<<", "<<(*nnf.getNearest())[0]<<std::endl;
			std::cout<<"Brute force: "<<minDist2<<" "<<(*nearest)[0]<<", "<<(*nearest)[1]<<", "<<(*nearest)[2]<<std::endl;
			}
		}
	#endif
	
	t.elapse();
	
	/* Clean up and exit: */
	std::cout<<"Loaded "<<lpo.getNumLoadedNodes()<<" nodes"<<std::endl;
	std::cout<<"Total runtime: "<<t.getTime()<<" s"<<std::endl;
	return 0;
	}

#endif
