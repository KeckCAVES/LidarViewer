/***********************************************************************
SpherePrimitive - Class for spheres extracted from point clouds.
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

#define NONSTANDARD_TEMPLATES 1

#include <iostream>
#include <Misc/ThrowStdErr.h>
#include <Misc/File.h>
#include <Comm/MulticastPipe.h>
#include <Math/Math.h>
#include <Geometry/Vector.h>
#include <GL/gl.h>
#include <GL/GLColorTemplates.h>
#include <GL/GLContextData.h>
#include <GL/GLModels.h>
#include <GL/GLGeometryWrappers.h>

#include "LidarOctree.h"
#include "LidarSelectionExtractor.h"
#include "SphereFitter.h"
#include "LevenbergMarquardtMinimizer.h"

#include "SpherePrimitive.h"

/******************************************
Methods of class SpherePrimitive::DataItem:
******************************************/

SpherePrimitive::DataItem::DataItem(void)
	:displayListId(glGenLists(2))
	{
	}

SpherePrimitive::DataItem::~DataItem(void)
	{
	glDeleteLists(displayListId,2);
	}

/********************************
Methods of class SpherePrimitive:
********************************/

SpherePrimitive::SpherePrimitive(const LidarOctree* octree,Comm::MulticastPipe* pipe)
	{
	/* Extract all selected points from the octree: */
	LidarSelectionExtractor<SphereFitter::Point> lse;
	octree->processSelectedPoints(lse);
	
	if(lse.getPoints().size()>=4)
		{
		/* Create a sphere fitter: */
		SphereFitter sf(lse.getPoints());
		
		/* Minimize the target function: */
		Scalar f=LevenbergMarquardtMinimizer<SphereFitter>::minimize(sf);
		
		/* Store the number of points and the RMS residual: */
		numPoints=lse.getPoints().size();
		rms=Math::sqrt(f*Scalar(2)/Scalar(numPoints));
		
		/* Extract the sphere parameters: */
		point=sf.getCenter();
		radius=sf.getRadius();
		
		/* Print the sphere's equation: */
		std::cout<<"Sphere fitting "<<numPoints<<" points"<<std::endl;
		std::cout<<"Center point: ("<<point[0]<<", "<<point[1]<<", "<<point[2]<<")"<<std::endl;
		std::cout<<"Radius: "<<radius<<std::endl;
		std::cout<<"RMS approximation residual: "<<rms<<std::endl;
		
		if(pipe!=0)
			{
			/* Send the extracted primitive over the pipe: */
			pipe->write<int>(1);
			pipe->write<unsigned int>((unsigned int)(numPoints));
			pipe->write<Scalar>(rms);
			pipe->write<Scalar>(point.getComponents(),3);
			pipe->write<Scalar>(radius);
			pipe->finishMessage();
			}
		}
	else
		{
		if(pipe!=0)
			{
			pipe->write<int>(0);
			pipe->finishMessage();
			}
		Misc::throwStdErr("SpherePrimitive::SpherePrimitive: Not enough selected points");
		}
	}

SpherePrimitive::SpherePrimitive(Comm::MulticastPipe* pipe)
	{
	/* Read the status flag from the pipe: */
	if(!pipe->read<int>())
		Misc::throwStdErr("SpherePrimitive::SpherePrimitive: Not enough selected points");
	
	/* Read the number of points and the RMS residual: */
	numPoints=pipe->read<unsigned int>();
	rms=pipe->read<Scalar>();
	
	/* Read the sphere parameters: */
	pipe->read<Scalar>(point.getComponents(),3);
	radius=pipe->read<Scalar>();
	}

SpherePrimitive::SpherePrimitive(Misc::File& file,const Primitive::Vector& translation)
	{
	/* Read the number of points and the RMS residual: */
	numPoints=file.read<unsigned int>();
	rms=file.read<Scalar>();
	
	/* Read the sphere parameters: */
	file.read<Scalar>(point.getComponents(),3);
	point+=translation;
	radius=file.read<Scalar>();
	}

Primitive::Scalar SpherePrimitive::pick(const Primitive::Point& pickPoint,Primitive::Scalar maxPickDistance) const
	{
	/* Calculate the pick point's distance from the sphere's center: */
	Scalar centerDist=Geometry::dist(pickPoint,point);
	
	/* Return the pick point's distance from the sphere: */
	return Math::abs(centerDist-radius);
	}

void SpherePrimitive::initContext(GLContextData& contextData) const
	{
	/* Create a data item and store it in the context: */
	DataItem* dataItem=new DataItem;
	contextData.addDataItem(this,dataItem);
	
	/* Create the sphere rendering display lists: */
	glNewList(dataItem->displayListId,GL_COMPILE);
	
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	glDrawSphereIcosahedron(1.0,5);
	
	glEndList();
	
	glNewList(dataItem->displayListId+1,GL_COMPILE);
	
	glBlendFunc(GL_ONE,GL_ONE);
	glDisable(GL_CULL_FACE);
	glLineWidth(1.0f);
	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	glDrawSphereIcosahedron(1.0,5);
	
	glEndList();
	}

void SpherePrimitive::glRenderAction(GLContextData& contextData) const
	{
	/* Retrieve the data item: */
	DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
	
	glPushAttrib(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_ENABLE_BIT|GL_LINE_BIT|GL_POINT_BIT|GL_POLYGON_BIT);
	glDisable(GL_LIGHTING);
	
	/* Draw the sphere's center: */
	glPointSize(5.0f);
	glColor(surfaceColor);
	glBegin(GL_POINTS);
	glVertex(point);
	glEnd();
	
	glEnable(GL_BLEND);
	glDepthMask(GL_FALSE);
	
	glPushMatrix();
	glTranslate(point-Point::origin);
	glScale(radius);
	
	/* Draw the surface: */
	glColor(surfaceColor);
	glCallList(dataItem->displayListId);
	
	/* Draw the grid: */
	glColor(gridColor);
	glCallList(dataItem->displayListId+1);
	
	glPopMatrix();
	
	glPopAttrib();
	}

void SpherePrimitive::write(Misc::File& file,const Primitive::Vector& translation) const
	{
	/* Write the number of points and the RMS residual: */
	file.write<unsigned int>((unsigned int)(numPoints));
	file.write<Scalar>(rms);
	
	/* Write the sphere parameters: */
	file.write<Scalar>((point+translation).getComponents(),3);
	file.write<Scalar>(radius);
	}
