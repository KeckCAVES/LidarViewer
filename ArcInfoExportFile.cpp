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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Misc/ThrowStdErr.h>
#include <Misc/File.h>
#include <GL/gl.h>
#include <GL/GLVertexArrayParts.h>
#include <GL/GLContextData.h>
#include <GL/GLExtensionManager.h>
#include <GL/Extensions/GLARBVertexBufferObject.h>
#include <GL/GLGeometryWrappers.h>

#include "ArcInfoExportFile.h"

namespace {

/****************
Helper functions:
****************/

inline ArcInfoExportFile::Point transformPoint(const ArcInfoExportFile::ReadTransform& transform,double px,double py)
	{
	return ArcInfoExportFile::Point(transform.transform(ArcInfoExportFile::ReadTransform::Point(px,py,0)));
	}

}

/********************************************
Methods of class ArcInfoExportFile::DataItem:
********************************************/

ArcInfoExportFile::DataItem::DataItem(void)
	:hasVertexBufferObjectExtension(GLARBVertexBufferObject::isSupported()),
	 polylineVerticesBufferObjectId(0)
	{
	/* Check if the vertex buffer object extension is supported: */
	if(hasVertexBufferObjectExtension)
		{
		/* Initialize the vertex buffer object extension: */
		GLARBVertexBufferObject::initExtension();
		
		/* Create a vertex buffer object: */
		glGenBuffersARB(1,&polylineVerticesBufferObjectId);
		}
	}

ArcInfoExportFile::DataItem::~DataItem(void)
	{
	/* Check if the vertex buffer object extension is supported: */
	if(hasVertexBufferObjectExtension)
		{
		/* Destroy the vertex buffer object: */
		glDeleteBuffersARB(1,&polylineVerticesBufferObjectId);
		}
	}


/**********************************
Methods of class ArcInfoExportFile:
**********************************/

void ArcInfoExportFile::calcBoundingBox(void)
	{
	/* Calculate the bounding box of all used polyline vertices: */
	bbox=Box::empty;
	for(std::vector<Polyline>::const_iterator plIt=polylines.begin();plIt!=polylines.end();++plIt)
		{
		const Point* vPtr=&polylineVertices[plIt->firstVertex];
		for(unsigned int i=0;i<plIt->numVertices;++i,++vPtr)
			bbox.addPoint(*vPtr);
		}
	}

ArcInfoExportFile::ArcInfoExportFile(const char* fileName,const ArcInfoExportFile::ReadTransform& readTransform)
	:bbox(Box::empty)
	{
	/* Open the input file: */
	Misc::File file(fileName,"rt");
	
	/* Check the file header: */
	char line[256];
	file.gets(line,sizeof(line));
	if(strncasecmp(line,"EXP ",4)!=0)
		Misc::throwStdErr("ArcInfoExportFile::ArcInfoExportFile: File %s is no valid ARC/INFO export file",fileName);
	if(atoi(line+4)!=0)
		Misc::throwStdErr("ArcInfoExportFile::ArcInfoExportFile: File %s is a compressed ARC/INFO export file",fileName);
	
	/* Read embedded ARC files until the end-of-file: */
	file.gets(line,sizeof(line));
	while(strncasecmp(line,"EOS",3)!=0)
		{
		int doubleFlag=atoi(line+3);
		if(strncasecmp(line,"ARC",3)==0)
			{
			/* Read an ARC file: */
			while(true)
				{
				/* Read the polyline header: */
				file.gets(line,sizeof(line));
				Polyline pl;
				sscanf(line,"%d %d %d %d %d %d %u",&pl.index,&pl.id,&pl.startNode,&pl.endNode,&pl.leftPolygonIndex,&pl.rightPolygonIndex,&pl.numVertices);
				if(pl.index==-1)
					break;
				
				/* Read the polyline vertices: */
				pl.firstVertex=polylineVertices.size();
				if(doubleFlag==2)
					{
					/* Single-precision points are two per line: */
					for(unsigned int i=0;i<pl.numVertices/2;++i)
						{
						file.gets(line,sizeof(line));
						double p1x,p1y,p2x,p2y;
						sscanf(line,"%lf %lf %lf %lf",&p1x,&p1y,&p2x,&p2y);
						polylineVertices.push_back(transformPoint(readTransform,p1x,p1y));
						polylineVertices.push_back(transformPoint(readTransform,p2x,p2y));
						}
					if(pl.numVertices%2==1)
						{
						file.gets(line,sizeof(line));
						double px,py;
						sscanf(line,"%lf %lf",&px,&py);
						polylineVertices.push_back(transformPoint(readTransform,px,py));
						}
					}
				else if(doubleFlag==3)
					{
					/* Double-precision points are one per line: */
					for(unsigned int i=0;i<pl.numVertices;++i)
						{
						file.gets(line,sizeof(line));
						double px,py;
						sscanf(line,"%lf %lf",&px,&py);
						polylineVertices.push_back(transformPoint(readTransform,px,py));
						}
					}
				
				/* Store the polyline: */
				polylines.push_back(pl);
				}
			}
		else if(strncasecmp(line,"SIN",3)==0)
			{
			/* Skip a SIN file: */
			do
				{
				file.gets(line,sizeof(line));
				}
			while(strncasecmp(line,"EOX",3)!=0);
			}
		else if(strncasecmp(line,"LOG",3)==0)
			{
			/* Skip a LOG file: */
			do
				{
				file.gets(line,sizeof(line));
				}
			while(strncasecmp(line,"EOL",3)!=0);
			}
		else if(strncasecmp(line,"PRJ",3)==0)
			{
			/* Skip a PRJ file: */
			do
				{
				file.gets(line,sizeof(line));
				}
			while(strncasecmp(line,"EOP",3)!=0);
			}
		else if(strncasecmp(line,"TX6",3)==0||strncasecmp(line,"TX7",3)==0||strncasecmp(line,"RXP",3)==0||strncasecmp(line,"RPL",3)==0)
			{
			/* Skip a text section: */
			do
				{
				file.gets(line,sizeof(line));
				}
			while(strncasecmp(line,"JABBERWOCKY",11)!=0); // Huh?
			}
		else if(strncasecmp(line,"MTD",3)==0)
			{
			/* Skip the Metadata section: */
			do
				{
				file.gets(line,sizeof(line));
				}
			while(strncasecmp(line,"EOD",3)!=0);
			}
		else if(strncasecmp(line,"IFO",3)==0)
			{
			/* Skip the INFO section: */
			do
				{
				file.gets(line,sizeof(line));
				}
			while(strncasecmp(line,"EOI",3)!=0);
			}
		else
			{
			/* Skip an unrecognized file: */
			while(true)
				{
				file.gets(line,sizeof(line));
				int in1,in2,in3,in4,in5,in6,in7;
				if(sscanf(line,"%d %d %d %d %d %d %d",&in1,&in2,&in3,&in4,&in5,&in6,&in7)==7&&in1==-1)
					break;
				}
			}
		
		/* Read the next file header: */
		file.gets(line,sizeof(line));
		}
	
	/* Compute the polylines' bounding box: */
	calcBoundingBox();
	}

ArcInfoExportFile::~ArcInfoExportFile(void)
	{
	}

void ArcInfoExportFile::initContext(GLContextData& contextData) const
	{
	/* Create a context data item and store it in the context: */
	DataItem* dataItem=new DataItem;
	contextData.addDataItem(this,dataItem);
	
	/* Check if the vertex buffer object extension is supported: */
	if(dataItem->hasVertexBufferObjectExtension)
		{
		/* Create a vertex buffer object to store the points' coordinates: */
		glBindBufferARB(GL_ARRAY_BUFFER_ARB,dataItem->polylineVerticesBufferObjectId);
		glBufferDataARB(GL_ARRAY_BUFFER_ARB,polylineVertices.size()*sizeof(Vertex),0,GL_STATIC_DRAW_ARB);
		
		/* Copy all points: */
		Vertex* vPtr=static_cast<Vertex*>(glMapBufferARB(GL_ARRAY_BUFFER_ARB,GL_WRITE_ONLY_ARB));
		for(std::vector<Point>::const_iterator pIt=polylineVertices.begin();pIt!=polylineVertices.end();++pIt,++vPtr)
			{
			for(int i=0;i<3;++i)
				vPtr->position[i]=(*pIt)[i];
			}
		glUnmapBufferARB(GL_ARRAY_BUFFER_ARB);
		
		/* Protect the vertex buffer object: */
		glBindBufferARB(GL_ARRAY_BUFFER_ARB,0);
		}
	}

void ArcInfoExportFile::glRenderAction(GLContextData& contextData) const
	{
	/* Get a pointer to the data item: */
	DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
	
	/* Save and set up OpenGL state: */
	GLVertexArrayParts::enable(Vertex::getPartsMask());
	
	/* Check if the vertex buffer object extension is supported: */
	if(dataItem->hasVertexBufferObjectExtension)
		{
		/* Bind the polylines' vertex buffer object: */
		glBindBufferARB(GL_ARRAY_BUFFER_ARB,dataItem->polylineVerticesBufferObjectId);
		
		/* Render the point set from the vertex buffer object: */
		glVertexPointer(static_cast<const Vertex*>(0));
		}
	else
		{
		/* Render the point set directly from the vertex vector: */
		glVertexPointer(0,&polylineVertices[0]);
		}
	
	/* Render all polylines: */
	for(std::vector<Polyline>::const_iterator plIt=polylines.begin();plIt!=polylines.end();++plIt)
		glDrawArrays(GL_LINE_STRIP,plIt->firstVertex,plIt->numVertices);
	
	if(dataItem->hasVertexBufferObjectExtension)
		{
		/* Protect the polylines' vertex buffer object: */
		glBindBufferARB(GL_ARRAY_BUFFER_ARB,0);
		}
	
	/* Restore OpenGL state: */
	GLVertexArrayParts::disable(Vertex::getPartsMask());
	}
