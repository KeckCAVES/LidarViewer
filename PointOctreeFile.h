/***********************************************************************
PointOctreeFile - Structures defining the file structure of octree-based
multiresolution point sets.
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

#ifndef POINTOCTREEFILE_INCLUDED
#define POINTOCTREEFILE_INCLUDED

#include <Misc/LargeFile.h>
#include <Geometry/Point.h>
#include <Geometry/ValuedPoint.h>
#include <Geometry/Box.h>
#include <GL/gl.h>
#include <GL/GLColor.h>

typedef Geometry::Point<float,3> Point; // Type for 3D points
typedef GLColor<GLubyte,4> Color; // Type for RGB color values
typedef Geometry::ValuedPoint<Point,Color> OctreePoint; // Type for octree points (3D point with RGB color)
typedef Geometry::Box<float,3> Box; // Type for domain bounding boxes
typedef Misc::LargeFile PointOctreeFile; // Type representing point octree files
typedef PointOctreeFile::Offset FileOffset; // Type for offsets in point octree files

struct PointOctreeFileHeader // Header structure for point octree files
	{
	/* Elements: */
	public:
	Point center; // Center of the cube containing all points
	float size; // Half side length of cube containing all points
	unsigned int maxNumPointsPerNode; // Maximum number of points stored in each octree node
	
	/* Constructors and destructors: */
	PointOctreeFileHeader(void) // Dummy constructor
		{
		};
	PointOctreeFileHeader(PointOctreeFile& file) // Reads header from file
		{
		file.read(center.getComponents(),3);
		file.read(size);
		file.read(maxNumPointsPerNode);
		};
	
	/* Methods: */
	static size_t getSize(void) // Returns the file size of the header
		{
		return sizeof(Point::Scalar)*3+sizeof(float)+sizeof(unsigned int);
		};
	void write(PointOctreeFile& file) const // Writes the header to file
		{
		file.write(center.getComponents(),3);
		file.write(size);
		file.write(maxNumPointsPerNode);
		};
	};

struct PointOctreeFileNode // Structure for a node in a point octree file
	{
	/* Elements: */
	public:
	float detailSize; // The size of the smallest detail stored in this node (distance between two closest points)
	unsigned int numPoints; // Number of points in this octree node (point data immediately follows the octree node header)
	FileOffset childrenOffset; // File position of the first of this node's eight children (child nodes are stored consecutively); 0 if node is leaf
	
	/* Constructors and destructors: */
	PointOctreeFileNode(void) // Dummy constructor
		{
		};
	PointOctreeFileNode(PointOctreeFile& file) // Reads node from file
		{
		file.read(detailSize);
		file.read(numPoints);
		file.read(childrenOffset);
		};
	
	/* Methods: */
	static size_t getSize(void) // Returns the file size of an octree node
		{
		return sizeof(float)+sizeof(unsigned int)+sizeof(FileOffset);
		};
	void write(PointOctreeFile& file) const // Writes the octree node to file
		{
		file.write(detailSize);
		file.write(numPoints);
		file.write(childrenOffset);
		};
	};

#endif
