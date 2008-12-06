/***********************************************************************
CylinderPrimitive - Class for cylinders extracted from point clouds.
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
#include <GL/gl.h>
#include <GL/GLColorTemplates.h>
#include <GL/GLContextData.h>
#include <GL/GLGeometryWrappers.h>

#include "LidarOctree.h"
#include "LidarSelectionExtractor.h"
#include "CylinderFitter.h"
#include "LevenbergMarquardtMinimizer.h"

#include "CylinderPrimitive.h"

/********************************************
Methods of class CylinderPrimitive::DataItem:
********************************************/

CylinderPrimitive::DataItem::DataItem(void)
	:displayListId(glGenLists(2))
	{
	}

CylinderPrimitive::DataItem::~DataItem(void)
	{
	glDeleteLists(displayListId,2);
	}

/**********************************
Methods of class CylinderPrimitive:
**********************************/

CylinderPrimitive::CylinderPrimitive(const LidarOctree* octree,Comm::MulticastPipe* pipe)
	{
	/* Extract all selected points from the octree: */
	LidarSelectionExtractor<CylinderFitter::Point> lse;
	octree->processSelectedPoints(lse);
	
	if(lse.getPoints().size()>=6)
		{
		/* Try to fit a cylinder starting with all the primary axes: */
		Scalar minF=Math::Constants<Scalar>::max;
		for(int initialAxis=0;initialAxis<3;++initialAxis)
			{
			/* Create a cylinder fitter: */
			CylinderFitter cf(lse.getPoints(),initialAxis);
			
			/* Minimize the target function: */
			Scalar f=LevenbergMarquardtMinimizer<CylinderFitter>::minimize(cf);
			
			if(minF>f)
				{
				/* Store the target function: */
				minF=f;
				
				/* Extract the cylinder parameters: */
				center=cf.getCenter();
				axis=cf.getAxis();
				axis.normalize();
				radius=cf.getRadius();
				}
			}
		
		/* Store the number of points and the RMS residual: */
		numPoints=lse.getPoints().size();
		rms=Math::sqrt(minF*Scalar(2)/Scalar(numPoints));
		
		/* Calculate the point set's coverage along the cylinder axis: */
		Scalar min,max;
		std::vector<CylinderFitter::Point>::const_iterator pIt=lse.getPoints().begin();
		min=max=(*pIt-center)*axis;
		for(++pIt;pIt!=lse.getPoints().end();++pIt)
			{
			Scalar d=(*pIt-center)*axis;
			if(min>d)
				min=d;
			else if(max<d)
				max=d;
			}
		
		/* Set the cylinder's height and adjust the center point: */
		length=(max-min)*Scalar(1.1);
		center+=axis*Math::mid(min,max);
		
		/* Print the cylinder's equation: */
		std::cout<<"Cylinder fitting "<<numPoints<<" points"<<std::endl;
		std::cout<<"Center point: ("<<center[0]<<", "<<center[1]<<", "<<center[2]<<")"<<std::endl;
		std::cout<<"Axis direction: ("<<axis[0]<<", "<<axis[1]<<", "<<axis[2]<<")"<<std::endl;
		std::cout<<"Radius: "<<radius<<", height: "<<length<<std::endl;
		std::cout<<"RMS approximation residual: "<<rms<<std::endl;
		
		/* Compute an appropriate number of grid lines in x and y: */
		double aspect=(Scalar(2)*Math::Constants<Scalar>::pi*radius)/length;
		if(aspect>=1.0)
			{
			numX=10;
			numY=int(Math::floor(10.0/aspect+0.5));
			}
		else
			{
			numY=10;
			numX=int(Math::floor(10.0*aspect+0.5));
			}
		
		if(pipe!=0)
			{
			/* Send the extracted primitive over the pipe: */
			pipe->write<int>(1);
			pipe->write<unsigned int>((unsigned int)(numPoints));
			pipe->write<Scalar>(rms);
			pipe->write<Scalar>(center.getComponents(),3);
			pipe->write<Scalar>(axis.getComponents(),3);
			pipe->write<Scalar>(radius);
			pipe->write<Scalar>(length);
			pipe->write<int>(numX);
			pipe->write<int>(numY);
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
		Misc::throwStdErr("CylinderPrimitive::CylinderPrimitive: Not enough selected points");
		}
	}

CylinderPrimitive::CylinderPrimitive(Comm::MulticastPipe* pipe)
	{
	/* Read the status flag from the pipe: */
	if(!pipe->read<int>())
		Misc::throwStdErr("CylinderPrimitive::CylinderPrimitive: Not enough selected points");
	
	/* Read the number of points and the RMS residual: */
	numPoints=pipe->read<unsigned int>();
	rms=pipe->read<Scalar>();
	
	/* Read the cylinder parameters: */
	pipe->read<Scalar>(center.getComponents(),3);
	pipe->read<Scalar>(axis.getComponents(),3);
	radius=pipe->read<Scalar>();
	length=pipe->read<Scalar>();
	numX=pipe->read<int>();
	numY=pipe->read<int>();
	}

CylinderPrimitive::CylinderPrimitive(Misc::File& file,const Primitive::Vector& translation)
	{
	/* Read the number of points and the RMS residual: */
	numPoints=file.read<unsigned int>();
	rms=file.read<Scalar>();
	
	/* Read the cylinder parameters: */
	file.read<Scalar>(center.getComponents(),3);
	center+=translation;
	file.read<Scalar>(axis.getComponents(),3);
	radius=file.read<Scalar>();
	length=file.read<Scalar>();
	numX=file.read<int>();
	numY=file.read<int>();
	}

Primitive::Scalar CylinderPrimitive::pick(const Primitive::Point& pickPoint,Primitive::Scalar maxPickDistance) const
	{
	Scalar dist2=Scalar(0);
	
	/* Project the pick point onto the cylinder's axis: */
	Scalar axisParam=Math::abs((pickPoint-center)*axis)-Math::div2(length);
	
	/* Reject if the axis parameter is out of bounds: */
	if(axisParam>maxPickDistance)
		return axisParam;
	if(axisParam>Scalar(0))
		dist2+=Math::sqr(axisParam);
	
	/* Calculate the pick point's distance from the cylinder's axis: */
	Scalar axisDist=Geometry::mag(Geometry::cross(axis,pickPoint-center));
	dist2+=Math::sqr(axisDist-radius);
	
	/* Return the pick point's distance from the cylinder: */
	return Math::sqrt(dist2);
	}

void CylinderPrimitive::initContext(GLContextData& contextData) const
	{
	/* Create a data item and store it in the context: */
	DataItem* dataItem=new DataItem;
	contextData.addDataItem(this,dataItem);
	
	/* Create a coordinate system for the cylinder: */
	Vector cx=Geometry::normal(axis);
	cx.normalize();
	Vector cy=Geometry::cross(axis,cx);
	cy.normalize();
	Vector cz=axis*(length/Scalar(2));
	
	/* Create the cylinder rendering display lists: */
	glNewList(dataItem->displayListId,GL_COMPILE);
	
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_CULL_FACE);
	glBegin(GL_QUAD_STRIP);
	glNormal(1.0,0.0,0.0);
	glVertex(center+cx*radius+cz);
	glVertex(center+cx*radius-cz);
	for(int x=1;x<72;++x)
		{
		Scalar angle=Math::rad(Scalar(x)*Scalar(360/72));
		Vector d=cx*Math::cos(angle)+cy*Math::sin(angle);
		glNormal(d);
		d*=radius;
		glVertex(center+d+cz);
		glVertex(center+d-cz);
		}
	glNormal(1.0,0.0,0.0);
	glVertex(center+cx*radius+cz);
	glVertex(center+cx*radius-cz);
	glEnd();
	
	glEndList();
	
	glNewList(dataItem->displayListId+1,GL_COMPILE);
	
	glBlendFunc(GL_ONE,GL_ONE);
	glLineWidth(1.0f);
	glBegin(GL_LINES);
	for(int x=0;x<numX;++x)
		{
		Scalar angle=Math::rad(Scalar(x)*Scalar(360)/Scalar(numX));
		Vector d=cx*Math::cos(angle)+cy*Math::sin(angle);
		d*=radius;
		glVertex(center+d+cz);
		glVertex(center+d-cz);
		}
	glEnd();
	for(int y=0;y<=numY;++y)
		{
		Point center2=center+axis*((Scalar(y)/Scalar(numY)-Scalar(0.5))*length);
		glBegin(GL_LINE_LOOP);
		for(int x=0;x<72;++x)
			{
			Scalar angle=Math::rad(Scalar(x)*Scalar(360/72));
			Vector d=cx*Math::cos(angle)+cy*Math::sin(angle);
			d*=radius;
			glVertex(center2+d);
			}
		glEnd();
		}
	
	glEndList();
	}

void CylinderPrimitive::glRenderAction(GLContextData& contextData) const
	{
	/* Retrieve the data item: */
	DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
	
	glPushAttrib(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_ENABLE_BIT|GL_LINE_BIT|GL_POLYGON_BIT);
	glDisable(GL_LIGHTING);
	
	/* Draw the cylinder's central axis: */
	glLineWidth(3.0f);
	glColor(surfaceColor);
	glBegin(GL_LINES);
	Vector z=axis*Math::div2(length);
	glVertex(center-z);
	glVertex(center+z);
	glEnd();
	
	glEnable(GL_BLEND);
	glDepthMask(GL_FALSE);
	
	/* Draw the surface: */
	glColor(surfaceColor);
	glCallList(dataItem->displayListId);
	
	/* Draw the grid: */
	glColor(gridColor);
	glCallList(dataItem->displayListId+1);
	
	glPopAttrib();
	}

void CylinderPrimitive::write(Misc::File& file,const Primitive::Vector& translation) const
	{
	/* Write the number of points and the RMS residual: */
	file.write<unsigned int>((unsigned int)(numPoints));
	file.write<Scalar>(rms);
	
	/* Write the cylinder parameters: */
	file.write<Scalar>((center+translation).getComponents(),3);
	file.write<Scalar>(axis.getComponents(),3);
	file.write<Scalar>(radius);
	file.write<Scalar>(length);
	file.write<int>(numX);
	file.write<int>(numY);
	}
