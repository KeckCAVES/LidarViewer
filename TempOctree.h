/***********************************************************************
TempOctree - Class to store points in a temporary octree for out-of-core
preprocessing of large point clouds.
Copyright (c) 2007-2013 Oliver Kreylos

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

#ifndef TEMPOCTREE_INCLUDED
#define TEMPOCTREE_INCLUDED

#include <deque>
#include <Threads/MutexCond.h>
#include <Threads/Thread.h>
#include <IO/StandardFile.h>

#include "LidarTypes.h"
#include "Cube.h"
#include "LidarFile.h"

class TempOctree
	{
	/* Embedded classes: */
	private:
	typedef IO::StandardFile File; // Type representing temporary octree files
	typedef File::Offset Offset; // Type for offsets in temporary octree files
	
	struct Node // Node of the temporary octree
		{
		/* Elements: */
		public:
		Cube domain; // Node's domain
		size_t numPoints; // Total number of points contained in the node's subtree
		union
			{
			LidarPoint* points; // Pointer to array of points contained in the node's subtree; used only during tree creation
			Offset pointsOffset; // Offset of node's points in temporary octree file if node is a leaf (0 for interior nodes)
			};
		Node* children; // Pointer to array of eight children if node is interior (0 for leaf nodes)
		
		/* Constructors and destructors: */
		Node(void) // Creates a leaf node
			:children(0)
			{
			};
		~Node(void) // Destroys node and its subtree
			{
			delete[] children;
			};
		
		/* Methods: */
		size_t estimateNumPointsInCube(const Cube& cube) const;
		size_t boundNumPointsInCube(const Cube& cube) const;
		};
	
	/* Elements: */
	char* tempFileName; // Name of the temporary octree file
	File file; // Handle of temporary octree file
	unsigned int maxNumPointsPerNode; // The maximum number of points that can be stored in each leaf node
	Box pointBbox; // Bounding box containing all points in this octree
	Node root; // Root of the temporary octree
	Threads::MutexCond writeQueueCond; // Condition variable to signal new entries in the write queue
	std::deque<Node*> writeQueue; // Queue of octree leaf nodes to be written to the octree file
	volatile bool writerThreadRun; // Flag to shut down the writer thread after the octree has been built
	Threads::Thread writerThread; // Thread to write leaf nodes' point sets to the octree file in the background
	
	/* Private methods: */
	void readLidarSubtree(Node& node,LidarFile::Offset childrenOffset,LidarFile& indexFile);
	void* writerThreadMethod(void);
	void createSubTree(Node& node);
	LidarPoint* getPointsInCube(Node& node,const Cube& cube,LidarPoint* points);
	
	/* Constructors and destructors: */
	public:
	TempOctree(char* fileNameTemplate,unsigned int sMaxNumPointsPerNode,LidarPoint* points,size_t numPoints); // Creates a temporary octree for the given array of points; shuffles point array in the process
	TempOctree(const char* lidarFileName); // Creates a TempOctree wrapper around an existing octree in .LiDAR format
	~TempOctree(void);
	
	/* Methods: */
	size_t getTotalNumPoints(void) const // Returns the total number of points in this octree
		{
		return root.numPoints;
		};
	const Box& getPointBbox(void) const // Returns the bounding box of all points in this octree
		{
		return pointBbox;
		};
	size_t estimateNumPointsInCube(const Cube& cube) const // Returns a lower bound on the number of points contained in the given cube
		{
		return root.estimateNumPointsInCube(cube);
		};
	size_t boundNumPointsInCube(const Cube& cube) const // Returns an upper bound on the number of points contained in the given cube
		{
		return root.boundNumPointsInCube(cube);
		};
	LidarPoint* getPointsInCube(const Cube& cube,LidarPoint* points) // Returns exactly the points contained in the given cube and stores them in the given point array; returns number of points written
		{
		return getPointsInCube(root,cube,points);
		};
	};

#endif
