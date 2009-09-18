/***********************************************************************
LidarProcessOctree - Class to process multiresolution LiDAR point sets.
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

#ifndef LIDARPROCESSOCTREE_INCLUDED
#define LIDARPROCESSOCTREE_INCLUDED

#define OPTIMIZED_TRAVERSAL 1
#define ALLOW_THREADING 1

#include <Math/Math.h>
#if ALLOW_THREADING
#include <Threads/Mutex.h>
#endif

#include "LidarTypes.h"
#include "Cube.h"
#include "LidarFile.h"

class LidarProcessOctree
	{
	/* Embedded classes: */
	public:
	class Node // Class for in-memory octree nodes
		{
		friend class LidarProcessOctree;
		
		/* Elements: */
		private:
		#if ALLOW_THREADING
		Threads::Mutex mutex; // Mutex to serialize changes to a node's state
		#endif
		Node* parent; // Pointer to node's parent node; 0 for root
		LidarFile::Offset childrenOffset; // Offset of the node's children in the octree file (0 if node is leaf node)
		Node* children; // Pointer to an array of eight child nodes (0 if node is not subdivided)
		Cube domain; // Node's domain
		unsigned int numPoints; // Number of LiDAR points belonging to this node
		Scalar detailSize; // Detail size of this node, for proper LOD computation
		LidarFile::Offset dataOffset; // Offset of the node's data in the octree's data file(s) in units of record size
		LidarPoint* points; // Pointer to the LiDAR points belonging to this node
		unsigned int processCounter; // Number of processes who have entered this node's subtree
		Node* lruPred; // Pointer to previous leaf parent node
		Node* lruSucc; // Pointer to next leaf parent node
		
		/* Private methods: */
		void enter(void)
			{
			#if ALLOW_THREADING
			Threads::Mutex::Lock lock(mutex);
			#endif
			++processCounter;
			}
		void leave(void)
			{
			#if ALLOW_THREADING
			Threads::Mutex::Lock lock(mutex);
			#endif
			--processCounter;
			}
		
		/* Constructors and destructors: */
		private:
		Node(void); // Creates a leaf node without points
		~Node(void); // Destroys a node and its subtree
		
		/* Methods: */
		public:
		bool isLeaf(void) const // Returns true if the node is a leaf in the out-of-core octree
			{
			return childrenOffset==LidarFile::Offset(0);
			}
		const Cube& getDomain(void) const // Returns the node's domain
			{
			return domain;
			}
		unsigned int getNumPoints(void) const // Returns the node's number of points
			{
			return numPoints;
			}
		Scalar getDetailSize(void) const // Returns the node's detail size
			{
			return detailSize;
			}
		LidarFile::Offset getDataOffset(void) const // Returns the node's data offset
			{
			return dataOffset;
			}
		const LidarPoint* getPoints(void) const // Returns the node's point array
			{
			return points;
			}
		const LidarPoint& operator[](unsigned int index) const // Returns one of the node's points
			{
			return points[index];
			}
		};
	
	/* Elements: */
	private:
	#if ALLOW_THREADING
	Threads::Mutex fileMutex; // Mutex to serialize access to the underlying LiDAR files
	#endif
	LidarFile indexFile; // The file containing the octree's structural data
	LidarFile pointsFile; // The file containing the octree's point data
	LidarFile::Offset pointsRecordSize; // Record size of points file
	size_t numNodes; // Total number of nodes in the octree
	unsigned int maxNumPointsPerNode; // Maximum number of points per node
	Node root; // Root node
	unsigned int cacheSize; // Maximum number of nodes that can be held in memory
	unsigned int numCachedNodes; // Number of nodes currently held in memory
	size_t numSubdivideCalls; // Total number of calls to the subdivide method
	size_t numLoadedNodes; // Total number of nodes that had to be loaded from file
	#if ALLOW_THREADING
	Threads::Mutex lruMutex; // Mutex to serialize access to the LRU list
	#endif
	Node* lruHead; // Head of leaf parent node list
	Node* lruTail; // Tail of leaf parent node list
	
	/* Private methods: */
	void subdivide(Node& node); // Subdivides the given parent node
	template <class NodeProcessFunctorParam>
	void processNodesPrefix(Node& node,unsigned int nodeLevel,NodeProcessFunctorParam& processFunctor); // Recursive node processor for octree in pre-fix order
	template <class NodeProcessFunctorParam>
	void processNodesPostfix(Node& node,unsigned int nodeLevel,NodeProcessFunctorParam& processFunctor); // Recursive node processor for octree in post-fix order
	template <class ProcessFunctorParam>
	void processPoints(Node& node,ProcessFunctorParam& processFunctor); // Recursive point processor
	template <class ProcessFunctorParam>
	void processPointsInBox(Node& node,const Box& box,ProcessFunctorParam& processFunctor); // Recursive point processor
	template <class DirectedProcessFunctorParam>
	static void processPointsDirectedKdtree(const LidarPoint* points,unsigned int left,unsigned int right,unsigned int splitDimension,DirectedProcessFunctorParam& processFunctor); // Recursive directed point processor for intra-node kd-tree
	template <class DirectedProcessFunctorParam>
	void processPointsDirectedOctree(Node& node,DirectedProcessFunctorParam& processFunctor); // Recursive directed point processor for octree
	
	/* Constructors and destructors: */
	public:
	LidarProcessOctree(const char* lidarFileName,size_t sCacheSize); // Creates octree from the given octree file; cache size is in bytes
	~LidarProcessOctree(void);
	
	/* Methods: */
	unsigned int getMaxNumPointsPerNode(void) const // Returns the maximum number of points stored in each octree node
		{
		return maxNumPointsPerNode;
		}
	const Cube& getDomain(void) const // Returns the LiDAR data set's domain
		{
		return root.domain;
		}
	Point getRootCenter(void) const; // Returns the center of the root node's domain
	Scalar getRootSize(void) const; // Returns the size of the root node's domain
	size_t getNumNodes(void) const // Returns the total number of nodes in the octree
		{
		return numNodes;
		}
	Node* getChild(Node* node,int childIndex) // Returns a child node of the given node
		{
		/* Check if the node is an interior node: */
		if(node->childrenOffset!=LidarFile::Offset(0))
			{
			/* Subdivide the node if necessary: */
			if(node->children==0)
				subdivide(*node);
			
			/* Return the requested child node: */
			return &(node->children[childIndex]);
			}
		else
			return 0;
		}
	template <class DirectedProcessFunctorParam>
	void processNodePointsDirected(Node* node,DirectedProcessFunctorParam& processFunctor) // Processes the points inside the given node in order of distance from a query point
		{
		if(node->numPoints>0)
			{
			/* Call recursive point processor on node's points: */
			processPointsDirectedKdtree(node->points,0,node->numPoints-1,0,processFunctor);
			}
		}
	template <class NodeProcessFunctorParam>
	void processNodesPrefix(NodeProcessFunctorParam& processFunctor) // Calls given functor for each octree node in the octree in pre-fix order
		{
		/* Call recursive node processor on the root node: */
		processNodesPrefix(root,0,processFunctor);
		}
	template <class NodeProcessFunctorParam>
	void processNodesPostfix(NodeProcessFunctorParam& processFunctor) // Calls given functor for each octree node in the octree in post-fix order
		{
		/* Call recursive node processor on the root node: */
		processNodesPostfix(root,0,processFunctor);
		}
	template <class ProcessFunctorParam>
	void processPoints(ProcessFunctorParam& processFunctor) // Calls given functor for each LiDAR point stored in the octree
		{
		/* Call recursive point processor on the root node: */
		processPoints(root,processFunctor);
		}
	template <class ProcessFunctorParam>
	void processPointsInBox(const Box& box,ProcessFunctorParam& processFunctor) // Calls given functor for each LiDAR point stored in the octree that is inside the given box
		{
		/* Call recursive point processor on the root node: */
		processPointsInBox(root,box,processFunctor);
		}
	template <class DirectedProcessFunctorParam>
	void processPointsDirected(DirectedProcessFunctorParam& processFunctor) // Calls given functor for each LiDAR point stored in the octree in order of distance from a query point
		{
		/* Call recursive directed point processor on the root node: */
		processPointsDirectedOctree(root,processFunctor);
		}
	size_t getNumSubdivideCalls(void) const // Returns total number of calls to the subdivide method
		{
		return numSubdivideCalls;
		}
	size_t getNumLoadedNodes(void) const // Returns total number of nodes loaded from file
		{
		return numLoadedNodes;
		}
	};

#if 0

/***************************************
Templates for point processing functors:
***************************************/

class NodeProcessFunctor
	{
	/* Methods: */
	public:
	void operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel); // Processes the given octree node
	};

class ProcessFunctor
	{
	/* Methods: */
	public:
	void operator()(const LidarPoint& point); // Processes the given LiDAR point
	};

class DirectedProcessFunctor
	{
	/* Methods: */
	public:
	void operator()(const LidarPoint& point); // Processes the given LiDAR point
	const Point& getQueryPoint(void) const; // Returns the query point directing the process
	Scalar getQueryRadius2(void) const; // Returns the square of the query sphere radius
	};

#endif

#ifndef LIDARPROCESSOCTREE_IMPLEMENTATION
#include "LidarProcessOctree.icpp"
#endif

#endif
