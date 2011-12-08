/***********************************************************************
LidarFile - Structures defining octrees of LiDAR point data.
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

#ifndef LIDARFILE_INCLUDED
#define LIDARFILE_INCLUDED

#include <IO/StandardFile.h>

#include "LidarTypes.h"
#include "Cube.h"

typedef IO::StandardFile LidarFile; // Type for LiDAR octree and data files

struct LidarOctreeFileHeader // Header structure for LiDAR octree files
	{
	/* Elements: */
	public:
	Cube domain; // Domain of the octree's root node
	unsigned int maxNumPointsPerNode; // Maximum number of points stored in each octree node
	
	/* Constructors and destructors: */
	LidarOctreeFileHeader(const Cube& sDomain,unsigned int sMaxNumPointsPerNode) // Creates an octree file header
		:domain(sDomain),
		 maxNumPointsPerNode(sMaxNumPointsPerNode)
		{
		}
	LidarOctreeFileHeader(LidarFile& file) // Reads the header of a LiDAR octree file
		:domain(file),
		 maxNumPointsPerNode(file.read<unsigned int>())
		{
		}
	
	/* Methods: */
	static size_t getFileSize(void) // Returns the file size of the header
		{
		return Cube::getFileSize()+sizeof(unsigned int);
		};
	void write(LidarFile& file) const // Writes the header to file
		{
		/* Write the octree's domain: */
		domain.write(file);
		
		/* Write the maximum number of points per node: */
		file.write(maxNumPointsPerNode);
		};
	};

struct LidarOctreeFileNode // Structure for a node in a LiDAR octree file
	{
	/* Elements: */
	public:
	LidarFile::Offset childrenOffset; // Offset of the node's children in the octree file (or 0 if leaf node)
	Scalar detailSize; // Detail size of this node, for proper LOD computation
	unsigned int numPoints; // Number of points in the node
	LidarFile::Offset dataOffset; // Offset of the node's points and ancillary data in the points or data files, in units of record size
	
	/* Methods: */
	static size_t getFileSize(void) // Returns the file size of an octree node
		{
		return sizeof(LidarFile::Offset)+sizeof(Scalar)+sizeof(unsigned int)+sizeof(LidarFile::Offset);
		};
	void read(LidarFile& file) // Reads the octree node from file
		{
		file.read(childrenOffset);
		file.read(detailSize);
		file.read(numPoints);
		file.read(dataOffset);
		};
	void write(LidarFile& file) const // Writes the octree node to file
		{
		file.write(childrenOffset);
		file.write(detailSize);
		file.write(numPoints);
		file.write(dataOffset);
		};
	};

struct LidarDataFileHeader // Header structure for LiDAR point or ancillary data files
	{
	/* Elements: */
	public:
	unsigned int recordSize; // Size of a data element in the data file in bytes
	
	/* Constructors and destructors: */
	LidarDataFileHeader(unsigned int sRecordSize) // Creates a LiDAR data file header
		:recordSize(sRecordSize)
		{
		}
	LidarDataFileHeader(LidarFile& file) // Reads the header of a LiDAR data file
		:recordSize(file.read<unsigned int>())
		{
		}
	
	/* Methods: */
	static size_t getFileSize(void) // Returns the file size of the header
		{
		return sizeof(unsigned int);
		};
	void write(LidarFile& file) const // Writes the header to file
		{
		file.write(recordSize);
		}
	};

#endif
