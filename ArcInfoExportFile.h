/***********************************************************************
ArcInfoExportFile - Class to represent GIS data parsed from ARC/INFO
export files in e00 format.
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

#ifndef ARCINFOEXPORTFILE_INCLUDED
#define ARCINFOEXPORTFILE_INCLUDED

#include <vector>
#include <Geometry/Point.h>
#include <Geometry/Box.h>
#include <Geometry/AffineTransformation.h>
#include <GL/gl.h>
#include <GL/GLVertex.h>
#include <GL/GLObject.h>

class ArcInfoExportFile:public GLObject
	{
	/* Embedded classes: */
	public:
	typedef float Scalar; // Scalar type for coordinates
	typedef Geometry::Point<Scalar,3> Point; // Type for coordinates
	typedef Geometry::Box<Scalar,3> Box; // Type for bounding boxes
	typedef Geometry::AffineTransformation<double,3> ReadTransform; // Type for coordinate transformations to apply during reading
	private:
	typedef GLVertex<void,0,void,0,void,Scalar,3> Vertex; // Type for vertices
	
	struct Polyline // Structure representing a polyline
		{
		/* Elements: */
		public:
		int index; // Coverage number
		int id; // Coverage ID
		int startNode; // ID of start node of polyline
		int endNode; // ID of end node of polyline
		int leftPolygonIndex; // Index of coverage polygon on left side of polyline
		int rightPolygonIndex; // Index of coverage polygon on right side of polyline
		unsigned int firstVertex; // Index of first polyline vertex in vertex array
		unsigned int numVertices; // Number of vertices in polyline
		};
	
	struct DataItem:public GLObject::DataItem
		{
		/* Elements: */
		public:
		bool hasVertexBufferObjectExtension; // Flag whether the local OpenGL supports vertex buffers
		GLuint polylineVerticesBufferObjectId; // ID of vertex buffer object containing polyline vertices
		
		/* Constructors and destructors: */
		DataItem(void);
		virtual ~DataItem(void);
		};
	
	/* Elements: */
	Box bbox; // Bounding box of all geometry elements
	std::vector<Point> polylineVertices; // Vector of vertices defining all polylines
	std::vector<Polyline> polylines; // Vector of polylines
	
	/* Private methods: */
	void calcBoundingBox(void); // Calculates the bounding box around all geometry parsed from the ARC/INFO export file
	
	/* Constructors and destructors: */
	public:
	ArcInfoExportFile(const char* fileName,const ReadTransform& readTransform =ReadTransform::identity); // Reads the ARC/INFO export file of the given name
	~ArcInfoExportFile(void);
	
	/* Methods: */
	virtual void initContext(GLContextData& contextData) const;
	const Box& getBoundingBox(void) const // Returns the bounding box around all geometry parsed from the ARC/INFO export file
		{
		return bbox;
		};
	template <class TransformParam>
	void transform(const TransformParam& transform) // Transforms all geometry parsed from the ARC/INFO export file
		{
		/* Transform all polyline vertices: */
		for(std::vector<Point>::iterator pIt=polylineVertices.begin();pIt!=polylineVertices.end();++pIt)
			*pIt=transform.transform(*pIt);
		
		/* Update the bounding box: */
		calcBoundingBox();
		}
	void glRenderAction(GLContextData& contextData) const; // Draws all geometry parsed from the ARC/INFO export file
	};

#endif
