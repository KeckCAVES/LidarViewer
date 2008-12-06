/***********************************************************************
PointOctreeCreator - Class to create point octrees from point clouds
using an out-of-core algorithm.
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

#ifndef POINTOCTREECREATOR_INCLUDED
#define POINTOCTREECREATOR_INCLUDED

#include <vector>
#include <string>

#include "PointOctreeFile.h"
#include "Cube.h"

/* Forward declarations: */
class TempPointOctree;

class PointOctreeCreator
	{
	/* Embedded classes: */
	public:
	typedef std::vector<TempPointOctree*> OctreeList; // Type for lists of temporary octrees
	
	private:
	struct Node // Structure for octree nodes during creation
		{
		/* Elements: */
		public:
		Node* children; // Pointer to array of eight child nodes (0 if node is leaf node)
		float detailSize; // Detail size of this node (distance between nearest neighbors)
		unsigned int numPoints; // Number of points belonging to this node
		OctreePoint* points; // Pointer to the points belonging to this node
		FileOffset pointsOffset; // Offset of the node's point set in the temporary point file
		FileOffset octreeNodeOffset; // Offset of node's structure data in the octree file
		FileOffset octreePointsOffset; // Offset of point data in the points file
		
		/* Constructors and destructors: */
		Node(void) // Creates an empty leaf node
			:children(0),
			 numPoints(0),points(0)
			{
			};
		~Node(void) // Destroys a node and its subtree
			{
			delete[] children;
			};
		};
	
	struct PointFile // Structure to store names and handles for temporary files holding created points for each tree level
		{
		/* Elements: */
		public:
		PointOctreeFile* file; // Pointer to the temporary point file
		std::string fileName; // Name of the temporary point file
		
		/* Constructors and destructors: */
		PointFile(void)
			:file(0)
			{
			};
		};
	
	/* Elements: */
	const OctreeList& octrees; // List of octrees containing the point set
	unsigned int maxNumCachablePoints; // Maximum number of points to be held in memory at any time
	unsigned int maxNumPointsPerNode; // Maximum number of points per leaf node
	unsigned int subsamplingFactor; // Subsampling factor for interior nodes
	Box domainBox; // Bounding box of all points in the point set
	Cube rootDomain; // The domain of the root node
	Node root; // The octree's root node
	std::string fileNameTemplate; // Template to generate temporary point file names
	std::vector<PointFile> pointFiles; // Vector of temporary point files for each tree level
	unsigned int totalNumPoints; // Total number of points in all input files
	unsigned int totalNumReadPoints; // Total number of points read from the temporary octrees; should equal the total number of points in the end
	unsigned int totalNumNodes; // Total number of nodes in the octree
	unsigned int maxLevel; // Number of the highest level in the octree
	unsigned int maxNumPointsPerInteriorNode; // Actual maximum number of points in any node in the octree
	unsigned int numWrittenNodes; // Number of nodes already written to the final octree files
	unsigned int nextNumWrittenNodesUpdate; // Next number of written nodes at which to print a percentage update
	
	/* Private methods: */
	void calcDetailSize(Node& node); // Calculates the detail size of the given node
	void writeNodePoints(Node& node,unsigned int level); // Writes a node's points to a temporary point file
	void subsample(Node& node,bool deleteAllPoints);
	void createSubTree(Node& node,const Cube& nodeDomain,unsigned int level);
	void createSubTree(Node& node,const Cube& nodeDomain,OctreePoint* points,unsigned int numPoints,unsigned int level);
	void calcFileOffsets(Node& node,unsigned int level,FileOffset& octreeFilePos,FileOffset& pointsFilePos);
	void writeSubtree(const Node& node,unsigned int level,PointOctreeFile& tempPointFile,PointOctreeFile& octreeFile,PointOctreeFile& pointsFile,OctreePoint* pointBuffer);
	
	/* Constructors and destructors: */
	public:
	PointOctreeCreator(const char* sFileNameTemplate,const OctreeList& sOctrees,unsigned int sMaxNumCachablePoints,unsigned int sMaxNumPointsPerNode,unsigned int sSubsamplingFactor); // Creates point octree for the union of the given point sets and the given node parameters
	~PointOctreeCreator(void);
	
	/* Methods: */
	void write(const char* octreeFileName,const char* pointsFileName); // Writes the octree structure and data to two files
	};

#endif
