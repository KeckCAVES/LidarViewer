/***********************************************************************
LidarOctreeCreator - Class to create LiDAR octrees from point clouds
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

#ifndef LIDAROCTREECREATOR_INCLUDED
#define LIDAROCTREECREATOR_INCLUDED

#include <vector>
#include <string>
#include <Misc/LargeFile.h>

#include "LidarTypes.h"
#include "Cube.h"
#include "LidarFile.h"

/* Forward declarations: */
class TempOctree;

class LidarOctreeCreator
	{
	/* Embedded classes: */
	public:
	typedef std::vector<TempOctree*> TempOctreeList; // Type for lists of temporary octrees
	private:
	typedef Misc::LargeFile File; // Type for temporary files used during octree creation
	
	struct Node // Structure for octree nodes during creation
		{
		/* Elements: */
		public:
		Node* children; // Pointer to array of eight child nodes (0 if node is leaf node)
		Scalar detailSize; // Detail size of this node (distance between nearest neighbors)
		unsigned int numPoints; // Number of points belonging to this node
		LidarPoint* points; // Pointer to the points belonging to this node
		File::Offset pointsOffset; // Offset of the node's point set in the temporary point file
		LidarFile::Offset octreeNodeOffset; // Offset of node's structure data in the resulting octree file
		LidarFile::Offset octreeDataOffset; // Offset of node's point and ancillary data in the resulting octree file, in units of data record size
		
		/* Constructors and destructors: */
		Node(void) // Creates an empty leaf node
			:children(0),
			 numPoints(0),points(0)
			{
			};
		~Node(void) // Destroys a node and its subtree
			{
			delete[] children;
			delete[] points;
			};
		};
	
	struct TempPointFile // Structure to store names and handles for temporary files holding created points for each tree level
		{
		/* Elements: */
		public:
		File* file; // Pointer to the temporary point file
		std::string fileName; // Name of the temporary point file
		
		/* Constructors and destructors: */
		TempPointFile(void)
			:file(0)
			{
			};
		};
	
	typedef std::vector<TempPointFile> TempPointFileList; // Type for lists of temporary point files
	
	/* Elements: */
	size_t maxNumCachablePoints; // Maximum number of points to be held in memory at any time
	unsigned int maxNumPointsPerNode; // Maximum number of points per leaf node
	const TempOctreeList& tempOctrees; // List of temporary octrees containing the point set
	Box domainBox; // Bounding box of all points in the point set
	Cube rootDomain; // The domain of the root node
	Node root; // The octree's root node
	std::string tempPointFileNameTemplate; // Template to generate temporary point file names
	TempPointFileList tempPointFiles; // Vector of temporary point files for each tree level
	size_t totalNumPoints; // Total number of points in all input files
	size_t totalNumReadPoints; // Total number of points read from the temporary octrees; should equal the total number of points in the end
	unsigned int totalNumNodes; // Total number of nodes in the octree
	unsigned int maxLevel; // Number of the highest level in the octree
	unsigned int maxNumPointsPerInteriorNode; // Actual maximum number of points in any node in the octree
	unsigned int numWrittenNodes; // Number of nodes already written to the final octree files
	unsigned int nextNumWrittenNodesUpdate; // Next number of written nodes at which to print a percentage update
	
	/* Private methods: */
	void writeNodePoints(Node& node,unsigned int level); // Writes a node's points to a temporary point file
	void subsample(Node& node,bool deleteAllPoints);
	void createSubTree(Node& node,const Cube& nodeDomain,unsigned int level);
	void createSubTreeWithPoints(Node& node,const Cube& nodeDomain,LidarPoint* points,size_t numPoints,unsigned int level);
	void calcFileOffsets(Node& node,unsigned int level,LidarFile::Offset& octreeFilePos,LidarFile::Offset& dataFilePos);
	void writeSubtree(const Node& node,unsigned int level,File& tempPointFile,LidarFile& octreeFile,LidarFile& pointsFile,LidarPoint* pointBuffer);
	
	/* Constructors and destructors: */
	public:
	LidarOctreeCreator(size_t sMaxNumCachablePoints,unsigned int sMaxNumPointsPerNode,const TempOctreeList& sTempOctrees,std::string sTempPointFileNameTemplate); // Creates point octree for the union of the given point sets and the given node parameters
	~LidarOctreeCreator(void);
	
	/* Methods: */
	void write(const char* lidarFileName); // Writes the octree structure and point data to the given LiDAR file directory
	};

#endif
