/***********************************************************************
LidarOctreeCreator - Class to represent Octrees of LiDAR data points to
generate multiresolution point sets.
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

#ifndef LIDAROCTREECREATOR_INCLUDED
#define LIDAROCTREECREATOR_INCLUDED

#include "LidarOctreeFile.h"

class LidarOctreeCreator
	{
	/* Embedded classes: */
	private:
	struct Node // Structure for Octree nodes
		{
		/* Elements: */
		public:
		Node* children; // Pointer to array of eight child nodes (0 if node is leaf node)
		unsigned int numPoints; // Number of LiDAR points belonging to this node
		LidarPoint* points; // Pointer to the LiDAR points belonging to this node (part of larger array for leaf nodes; private array for interior nodes)
		Misc::LargeFile::Offset nodeOffset; // Offset of node's structure data in the octree file
		Misc::LargeFile::Offset pointsOffset; // Offset of point data in the points file
		
		/* Constructors and destructors: */
		Node(void) // Creates a leaf node without points
			:children(0),
			 numPoints(0),points(0)
			{
			};
		~Node(void); // Destroys a node and its subtree
		};
	
	struct Traversal // Structure to traverse an octree recursively
		{
		/* Elements: */
		public:
		Node* node; // Pointer to the traversed node
		Point center; // Center point of the traversed node
		float radius; // Radius (half side length) of the traversed node's domain
		};
	
	/* Elements: */
	Node* root; // Pointer to the octree's root node
	Traversal rootTraversal; // A traversal structure for the octree's root node
	unsigned int maxPointsPerNode; // Maximum number of points per leaf node
	unsigned int subsamplingFactor; // Subsampling factor for interior nodes
	unsigned int numNodes; // Total number of nodes in the tree
	unsigned int maxNumPointsPerInteriorNode; // Maximum number of points per interior node
	int numLevels; // Number of levels in the tree (including the root node level)
	
	/* Private methods: */
	unsigned int splitPoints(LidarPoint* points,unsigned int numPoints,int dimension,float split); // Splits an array of points along the center position in the given dimension
	unsigned int createSubtree(const Traversal& t,LidarPoint* points,unsigned int numPoints,unsigned int treeLevel); // Creates the subtree for the traversed node
	bool calcFileOffsets(Node* node,int level,Misc::LargeFile::Offset& octreeFilePos,Misc::LargeFile::Offset& pointsFilePos); // Calculates the file offsets of a node's subtree
	void writeSubtree(const Node* node,int level,Misc::LargeFile& octreeFile,Misc::LargeFile& pointsFile) const; // Writes a node's subtrees' structure and data to two files
	
	/* Constructors and destructors: */
	public:
	LidarOctreeCreator(LidarPoint* points,unsigned int numPoints,unsigned int sMaxPointsPerNode,unsigned int sSubsamplingFactor); // Creates octree for given point set and max node size
	~LidarOctreeCreator(void);
	
	/* Methods: */
	void write(const char* octreeFileName,const char* pointsFileName) const; // Writes the octree structure and data to two files
	};

#endif
