/***********************************************************************
LidarOctreeFile - Structures defining octrees of LiDAR point data.
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

#ifndef LIDAROCTREEFILE_INCLUDED
#define LIDAROCTREEFILE_INCLUDED

#include <Misc/Endianness.h>
#include <Misc/LargeFile.h>
#include <Geometry/Point.h>
#include <Geometry/ValuedPoint.h>
#include <Geometry/Endianness.h>
#include <GL/gl.h>
#include <GL/GLColor.h>

typedef float Scalar; // Scalar type for 3D points
typedef Geometry::Point<Scalar,3> Point; // Type for 3D points
typedef GLColor<GLubyte,4> Color; // Type for RGB color values
typedef Geometry::ValuedPoint<Point,Color> LidarPoint; // Type for LiDAR measurement points (3D point with RGB color)

struct LidarOctreeFileHeader // Header structure for LiDAR octree files
	{
	/* Elements: */
	public:
	Point center; // Center of the cube containing all LiDAR points
	Scalar radius; // Radius (half side length) of cube containing all LiDAR points
	unsigned int maxNumPointsPerNode; // Maximum number of points stored in each octree node
	
	/* Methods: */
	static size_t getSize(void) // Returns the file size of the header
		{
		return sizeof(Point)+sizeof(Scalar)+sizeof(unsigned int);
		};
	void write(Misc::LargeFile& file) const // Writes the header to file
		{
		file.write(center.getComponents(),3);
		file.write(radius);
		file.write(maxNumPointsPerNode);
		};
	void read(Misc::LargeFile& file) // Reads the header from file
		{
		file.read(center.getComponents(),3);
		file.read(radius);
		file.read(maxNumPointsPerNode);
		};
	};

struct LidarOctreeFileNode // Structure for a node in a LiDAR octree file
	{
	/* Elements: */
	public:
	Misc::LargeFile::Offset childrenOffset; // Offset of the node's children in the octree file (or 0 if leaf node)
	Scalar detailSize; // Detail size of this node, for proper LOD computation
	unsigned int numPoints; // Number of points in the node
	Misc::LargeFile::Offset pointsOffset; // Offset of the node's points in the points file
	
	/* Methods: */
	static size_t getSize(void) // Returns the file size of an octree node
		{
		return sizeof(Misc::LargeFile::Offset)+sizeof(Scalar)+sizeof(Misc::LargeFile::Offset)+sizeof(unsigned int);
		// return sizeof(Misc::LargeFile::Offset)+sizeof(Misc::LargeFile::Offset)+sizeof(unsigned int);
		};
	void write(Misc::LargeFile& file) const // Writes the octree node to file
		{
		file.write(childrenOffset);
		file.write(detailSize);
		file.write(pointsOffset);
		file.write(numPoints);
		};
	void read(Misc::LargeFile& file) // Reads the octree node from file
		{
		file.read(childrenOffset);
		file.read(detailSize);
		file.read(pointsOffset);
		file.read(numPoints);
		};
	};

/******************************************
Helper functions for endianness conversion:
******************************************/

namespace Misc {

template <>
inline
void
swapEndianness(
	LidarOctreeFileHeader& ofh)
	{
	swapEndianness(ofh.center);
	swapEndianness(ofh.radius);
	swapEndianness(ofh.maxNumPointsPerNode);
	}

template <>
inline
void
swapEndianness(
	LidarOctreeFileNode& ofn)
	{
	swapEndianness(ofn.childrenOffset);
	swapEndianness(ofn.detailSize);
	swapEndianness(ofn.numPoints);
	swapEndianness(ofn.pointsOffset);
	}

template <>
class EndiannessSwapper<LidarPoint>
	{
	/* Methods: */
	public:
	static void swap(LidarPoint& lp)
		{
		swapEndianness(lp.getComponents(),3);
		}
	static void swap(LidarPoint* lps,size_t numLps)
		{
		for(size_t i=0;i<numLps;++i)
			swapEndianness(lps[i].getComponents(),3);
		}
	};

}

#endif
