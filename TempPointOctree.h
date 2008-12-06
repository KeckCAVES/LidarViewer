/***********************************************************************
TempPointOctree - Class to store points in a temporary octree for
out-of-core preprocessing of large point clouds.
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

#ifndef TEMPPOINTOCTREE_INCLUDED
#define TEMPPOINTOCTREE_INCLUDED

#include <vector>

#include "PointOctreeFile.h"
#include "Cube.h"

class TempPointOctree
	{
	/* Embedded classes: */
	private:
	struct Node // Node of the temporary octree
		{
		/* Elements: */
		public:
		Cube domain; // Node's domain
		unsigned int numPoints; // Total number of points contained in the node's subtree
		FileOffset pointsOffset; // Offset of node's points in temporary octree file if node is a leaf (0 for interior nodes)
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
		unsigned int estimateNumPointsInCube(const Cube& cube) const;
		unsigned int boundNumPointsInCube(const Cube& cube) const;
		};
	
	/* Elements: */
	char* tempFileName; // Name of the temporary octree file
	PointOctreeFile file; // Handle of temporary octree file
	unsigned int maxNumPointsPerNode; // The maximum number of points that can be stored in each leaf node
	Box pointBbox; // Bounding box containing all points in this octree
	Node root; // Root of the temporary octree
	
	/* Private methods: */
	void createSubTree(Node& node,OctreePoint* points,unsigned int numPoints);
	void loadSubTree(Node& node,PointOctreeFile& octFile,FileOffset nodeOffset);
	void getPointsInCube(Node& node,const Cube& cube,std::vector<OctreePoint>& points);
	static const char* getObinFileName(const char* octreeFileNameStem);
	
	/* Constructors and destructors: */
	public:
	TempPointOctree(char* fileNameTemplate,unsigned int sMaxNumPointsPerNode,OctreePoint* points,unsigned int numPoints); // Creates a temporary octree for the given array of points; shuffles point array in the process
	TempPointOctree(const char* octreeFileNameStem,bool newOctreeFileFormat); // Creates a temporary octree from an existing LiDAR octree file pair
	~TempPointOctree(void);
	
	/* Methods: */
	unsigned int getTotalNumPoints(void) const // Returns the total number of points in this octree
		{
		return root.numPoints;
		};
	const Box& getPointBbox(void) const // Returns the bounding box of all points in this octree
		{
		return pointBbox;
		};
	unsigned int estimateNumPointsInCube(const Cube& cube) const // Gives a lower limit on the number of points contained in the given cube
		{
		return root.estimateNumPointsInCube(cube);
		};
	unsigned int boundNumPointsInCube(const Cube& cube) const // Returns an upper bound on the number of points contained in the given cube
		{
		return root.boundNumPointsInCube(cube);
		};
	void getPointsInCube(const Cube& cube,std::vector<OctreePoint>& points) // Returns exactly the points contained in the given cube
		{
		getPointsInCube(root,cube,points);
		};
	};

#endif
