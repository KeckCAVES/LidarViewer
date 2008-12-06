/***********************************************************************
PlanePrimitive - Class for planes extracted from point clouds.
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

#include <iostream>
#include <Misc/Utility.h>
#include <Misc/ThrowStdErr.h>
#include <Misc/File.h>
#include <Comm/MulticastPipe.h>
#include <Math/Math.h>
#include <Geometry/Vector.h>
#include <GL/gl.h>
#include <GL/GLColorTemplates.h>
#include <GL/GLContextData.h>
#include <GL/GLGeometryWrappers.h>

#include "LidarOctree.h"
#include "LidarPlaneExtractor.h"
#include "LidarPlaneFitter.h"

#include "PlanePrimitive.h"

/*****************************************
Methods of class PlanePrimitive::DataItem:
*****************************************/

PlanePrimitive::DataItem::DataItem(void)
	:displayListId(glGenLists(2))
	{
	}

PlanePrimitive::DataItem::~DataItem(void)
	{
	glDeleteLists(displayListId,2);
	}

/*******************************
Methods of class PlanePrimitive:
*******************************/

PlanePrimitive::PlanePrimitive(const LidarOctree* octree,Comm::MulticastPipe* pipe)
	{
	/* Create a LiDAR plane extractor: */
	LidarPlaneExtractor lpe;
	
	/* Process all selected points: */
	octree->processSelectedPoints(lpe);
	
	if(lpe.getNumPoints()>=3)
		{
		/* Extract the plane's coordinate frame: */
		LidarPlaneExtractor::Point centroid;
		LidarPlaneExtractor::Vector planeFrame[3];
		double lengths[3];
		lpe.calcPlane(centroid,planeFrame,lengths);
		
		/* Ensure that (planeFrame, planeNormal) is a right-handed system: */
		if(Geometry::cross(planeFrame[0],planeFrame[1])*planeFrame[2]<Scalar(0))
			planeFrame[2]=-planeFrame[2];
		plane=Plane(planeFrame[2],centroid);
		
		/* Calculate the bounding box of the selected points in plane coordinates: */
		LidarPlaneFitter lpf(centroid,planeFrame);
		octree->processSelectedPoints(lpf);
		
		/* Store the number of points and the RMS residual: */
		numPoints=lpe.getNumPoints();
		rms=lpf.getRMS();
		
		/* Calculate the extracted plane's rectangle: */
		double min[2],max[2];
		for(int i=0;i<2;++i)
			{
			min[i]=lpf.getMin(i);
			max[i]=lpf.getMax(i);
			}
		double size=max[0]-min[0];
		if(size<max[1]-min[1])
			size=max[1]-min[1];
		for(int i=0;i<2;++i)
			{
			min[i]-=0.1*size;
			max[i]+=0.1*size;
			}
		for(int i=0;i<4;++i)
			{
			points[i]=centroid;
			points[i]+=planeFrame[0]*(i&0x1?max[0]:min[0]);
			points[i]+=planeFrame[1]*(i&0x2?max[1]:min[1]);
			}
		
		/* Compute an appropriate number of grid lines in x and y: */
		double aspect=(max[0]-min[0])/(max[1]-min[1]);
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
		
		/* Print the plane equation: */
		std::cout<<"Plane fitting "<<numPoints<<" points"<<std::endl;
		std::cout<<"Centroid: ("<<centroid[0]<<", "<<centroid[1]<<", "<<centroid[2]<<")"<<std::endl;
		LidarPlaneExtractor::Vector normal=planeFrame[2];
		normal.normalize();
		std::cout<<"Normal vector: ("<<normal[0]<<", "<<normal[1]<<", "<<normal[2]<<")"<<std::endl;
		std::cout<<"RMS approximation residual: "<<rms<<std::endl;
		
		if(pipe!=0)
			{
			/* Send the extracted primitive over the pipe: */
			pipe->write<int>(1);
			pipe->write<unsigned int>((unsigned int)(numPoints));
			pipe->write<Scalar>(rms);
			pipe->write<LidarPlaneExtractor::Point::Scalar>(centroid.getComponents(),3);
			pipe->write<LidarPlaneExtractor::Vector::Scalar>(planeFrame[2].getComponents(),3);
			for(int i=0;i<4;++i)
				pipe->write<Point::Scalar>(points[i].getComponents(),3);
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
		Misc::throwStdErr("PlanePrimitive::PlanePrimitive: Not enough selected points");
		}
	}

PlanePrimitive::PlanePrimitive(Comm::MulticastPipe* pipe)
	{
	/* Read the status flag from the pipe: */
	if(!pipe->read<int>())
		Misc::throwStdErr("PlanePrimitive::PlanePrimitive: Not enough selected points");
	
	/* Read the number of points and the RMS residual: */
	numPoints=pipe->read<unsigned int>();
	rms=pipe->read<Scalar>();
	
	/* Read the plane equation: */
	LidarPlaneExtractor::Point centroid;
	LidarPlaneExtractor::Vector normal;
	pipe->read<LidarPlaneExtractor::Point::Scalar>(centroid.getComponents(),3);
	pipe->read<LidarPlaneExtractor::Vector::Scalar>(normal.getComponents(),3);
	plane=Plane(normal,centroid);
	for(int i=0;i<4;++i)
		pipe->read<Point::Scalar>(points[i].getComponents(),3);
	numX=pipe->read<int>();
	numY=pipe->read<int>();
	}

PlanePrimitive::PlanePrimitive(Misc::File& file,const Primitive::Vector& translation)
	{
	/* Read the number of points and the RMS residual: */
	numPoints=file.read<unsigned int>();
	rms=file.read<Scalar>();
	
	/* Read the plane equation: */
	Plane::Vector normal;
	Scalar offset;
	file.read<Scalar>(normal.getComponents(),3);
	offset=file.read<Scalar>();
	offset+=normal*translation;
	plane=Plane(normal,offset);
	for(int i=0;i<4;++i)
		{
		file.read<Point::Scalar>(points[i].getComponents(),3);
		points[i]+=translation;
		}
	numX=file.read<int>();
	numY=file.read<int>();
	}

Primitive::Scalar PlanePrimitive::pick(const Primitive::Point& pickPoint,Primitive::Scalar maxPickDistance) const
	{
	Scalar dist2=Scalar(0);
	
	/* Calculate the pick point's distance from the plane: */
	Scalar planeDist=Math::abs(plane.calcDistance(pickPoint));
	
	/* Reject if too far away: */
	if(planeDist>=maxPickDistance)
		return planeDist;
	dist2+=Math::sqr(planeDist);
	
	/* Check if the pick point is within the bounds of the plane primitive's rectangle: */
	Scalar dist;
	if((dist=(pickPoint-points[0])*Geometry::cross(points[1]-points[0],plane.getNormal()))>Scalar(0))
		dist2+=Math::sqr(dist);
	if((dist=(pickPoint-points[1])*Geometry::cross(points[3]-points[1],plane.getNormal()))>Scalar(0))
		dist2+=Math::sqr(dist);
	if((dist=(pickPoint-points[3])*Geometry::cross(points[2]-points[3],plane.getNormal()))>Scalar(0))
		dist2+=Math::sqr(dist);
	if((dist=(pickPoint-points[2])*Geometry::cross(points[0]-points[2],plane.getNormal()))>Scalar(0))
		dist2+=Math::sqr(dist);
	
	return Math::sqrt(dist2);
	}

void PlanePrimitive::initContext(GLContextData& contextData) const
	{
	/* Create a data item and store it in the context: */
	DataItem* dataItem=new DataItem;
	contextData.addDataItem(this,dataItem);
	
	/* Create the plane rendering display lists: */
	glNewList(dataItem->displayListId,GL_COMPILE);
	
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_CULL_FACE);
	glBegin(GL_QUADS);
	glVertex(points[0]);
	glVertex(points[1]);
	glVertex(points[3]);
	glVertex(points[2]);
	glEnd();
	
	glEndList();
	
	glNewList(dataItem->displayListId+1,GL_COMPILE);
	
	glBlendFunc(GL_ONE,GL_ONE);
	glLineWidth(1.0f);
	glBegin(GL_LINES);
	for(int x=0;x<=numX;++x)
		{
		glVertex(Geometry::affineCombination(points[0],points[1],Scalar(x)/Scalar(numX)));
		glVertex(Geometry::affineCombination(points[2],points[3],Scalar(x)/Scalar(numX)));
		}
	for(int y=0;y<=numY;++y)
		{
		glVertex(Geometry::affineCombination(points[0],points[2],Scalar(y)/Scalar(numY)));
		glVertex(Geometry::affineCombination(points[1],points[3],Scalar(y)/Scalar(numY)));
		}
	glEnd();
	
	glEndList();
	}

void PlanePrimitive::glRenderAction(GLContextData& contextData) const
	{
	/* Retrieve the data item: */
	DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
	
	glPushAttrib(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_ENABLE_BIT|GL_LINE_BIT|GL_POLYGON_BIT);
	glDisable(GL_LIGHTING);
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

void PlanePrimitive::write(Misc::File& file,const Primitive::Vector& translation) const
	{
	/* Write the number of points and the RMS residual: */
	file.write<unsigned int>((unsigned int)(numPoints));
	file.write<Scalar>(rms);
	
	/* Write the plane equation: */
	file.write<Scalar>(plane.getNormal().getComponents(),3);
	file.write<Scalar>(plane.getOffset()+plane.getNormal()*translation);
	for(int i=0;i<4;++i)
		file.write<Point::Scalar>((points[i]+translation).getComponents(),3);
	file.write<int>(numX);
	file.write<int>(numY);
	}
