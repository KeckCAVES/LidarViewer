/***********************************************************************
LinePrimitive - Class for lines extracted from point clouds by
intersecting two plane primitives.
Copyright (c) 2008 Oliver Kreylos

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
#include <Math/Constants.h>
#include <GL/gl.h>
#include <GL/GLColorTemplates.h>
#include <GL/GLGeometryWrappers.h>

#include "PlanePrimitive.h"

#include "LinePrimitive.h"

/******************************
Methods of class LinePrimitive:
******************************/

LinePrimitive::LinePrimitive(const PlanePrimitive* p1,const PlanePrimitive* p2,Comm::MulticastPipe* pipe)
	{
	/* Get the two planes' plane equations: */
	const PlanePrimitive* ps[2];
	ps[0]=p1;
	ps[1]=p2;
	PlanePrimitive::Plane planes[2];
	for(int i=0;i<2;++i)
		planes[i]=ps[i]->getPlane();
	
	/* Create the underdetermined linear system: */
	double a[3][4];
	for(int i=0;i<2;++i)
		{
		for(int j=0;j<3;++j)
			a[i][j]=double(planes[i].getNormal()[j]);
		a[i][3]=double(planes[i].getOffset());
		}
	for(int j=0;j<4;++j)
		a[2][j]=0.0;
	
	/* Find the null space of the underdetermined system: */
	int rowIndices[3];
	for(int i=0;i<3;++i)
		rowIndices[i]=i;
	for(int step=0;step<3-1;++step)
		{
		/* Find the full pivot: */
		double pivot=Math::abs(a[step][step]);
		int pivotRow=step;
		int pivotCol=step;
		for(int i=step;i<3;++i)
			for(int j=step;j<3;++j)
				{
				double val=Math::abs(a[i][j]);
				if(pivot<val)
					{
					pivot=val;
					pivotRow=i;
					pivotCol=j;
					}
				}
		
		/* Swap current and pivot rows if necessary: */
		if(pivotRow!=step)
			{
			/* Swap rows step and pivotRow: */
			for(int j=0;j<4;++j)
				Misc::swap(a[step][j],a[pivotRow][j]);
			}
		
		/* Swap current and pivot columns if necessary: */
		if(pivotCol!=step)
			{
			/* Swap columns step and pivotCol: */
			for(int i=0;i<3;++i)
				Misc::swap(a[i][step],a[i][pivotCol]);
			Misc::swap(rowIndices[step],rowIndices[pivotCol]);
			}
		
		/* Combine all rows with the current row: */
		for(int i=step+1;i<3;++i)
			{
			/* Combine rows i and step: */
			double factor=-a[i][step]/a[step][step];
			for(int j=step+1;j<4;++j)
				a[i][j]+=a[step][j]*factor;
			}
		}
	
	/* Calculate the swizzled result using backsubstitution: */
	double x[3],y[3];
	x[3-1]=1.0;
	y[3-1]=0.0;
	for(int i=3-2;i>=0;--i)
		{
		x[i]=0.0;
		y[i]=a[i][3];
		for(int j=i+1;j<3;++j)
			{
			x[i]-=a[i][j]*x[j];
			y[i]-=a[i][j]*y[j];
			}
		x[i]/=a[i][i];
		y[i]/=a[i][i];
		}
	
	/* Unswizzle the result: */
	for(int i=0;i<3;++i)
		{
		axis[rowIndices[i]]=Scalar(x[i]);
		center[rowIndices[i]]=Scalar(y[i]);
		}
	
	Scalar axisLen=Geometry::mag(axis);
	if(axisLen>Scalar(0))
		{
		/* Normalize the line direction: */
		axis/=axisLen;
		
		/* Find the extents of both planes' rectangles on the line: */
		double min=Math::Constants<double>::max;
		double max=Math::Constants<double>::min;
		for(int plane=0;plane<2;++plane)
			for(int i=0;i<4;++i)
				{
				double param=(ps[plane]->getPoint(i)-center)*axis;
				if(min>param)
					min=param;
				if(max<param)
					max=param;
				}
		
		/* Calculate the line's length, and adjust the center point: */
		length=(max-min)*Scalar(1.1);
		center+=axis*Math::mid(min,max);
		
		/* Print the line's equation: */
		std::cout<<"Line intersecting two planes"<<std::endl;
		std::cout<<"Center point: ("<<center[0]<<", "<<center[1]<<", "<<center[2]<<")"<<std::endl;
		std::cout<<"Axis direction: ("<<axis[0]<<", "<<axis[1]<<", "<<axis[2]<<")"<<std::endl;
		std::cout<<"Length: "<<length<<std::endl;
		
		if(pipe!=0)
			{
			/* Send the extracted primitive over the pipe: */
			pipe->write<int>(1);
			pipe->write<Scalar>(center.getComponents(),3);
			pipe->write<Scalar>(axis.getComponents(),3);
			pipe->write<Scalar>(length);
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
		Misc::throwStdErr("LinePrimitive::LinePrimitive: Given planes do not intersect");
		}
	}

LinePrimitive::LinePrimitive(Comm::MulticastPipe* pipe)
	{
	/* Read the status flag from the pipe: */
	if(!pipe->read<int>())
		Misc::throwStdErr("LinePrimitive::LinePrimitive: Given planes do not intersect");
	
	/* Read the line parameters: */
	pipe->read<Scalar>(center.getComponents(),3);
	pipe->read<Scalar>(axis.getComponents(),3);
	length=pipe->read<Scalar>();
	}

LinePrimitive::LinePrimitive(Misc::File& file,const Primitive::Vector& translation)
	{
	/* Read the line parameters: */
	file.read<Scalar>(center.getComponents(),3);
	center+=translation;
	file.read<Scalar>(axis.getComponents(),3);
	length=file.read<Scalar>();
	}

Primitive::Scalar LinePrimitive::pick(const Primitive::Point& pickPoint,Primitive::Scalar maxPickDistance) const
	{
	Scalar dist2=Scalar(0);
	
	/* Project the pick point onto the line's axis: */
	Scalar axisParam=Math::abs((pickPoint-center)*axis)-Math::div2(length);
	
	/* Reject if the axis parameter is out of bounds: */
	if(axisParam>maxPickDistance)
		return axisParam;
	if(axisParam>Scalar(0))
		dist2+=Math::sqr(axisParam);
	
	/* Calculate the pick point's distance from the line's axis: */
	Scalar axisDist2=Geometry::sqr(Geometry::cross(axis,pickPoint-center));
	dist2+=axisDist2;
	
	/* Return the pick point's distance from the line: */
	return Math::sqrt(dist2);
	}

void LinePrimitive::glRenderAction(GLContextData& contextData) const
	{
	glPushAttrib(GL_COLOR_BUFFER_BIT|GL_ENABLE_BIT|GL_LINE_BIT);
	glDisable(GL_LIGHTING);
	
	/* Draw the cylinder's central axis: */
	glLineWidth(3.0f);
	glColor(surfaceColor);
	glBegin(GL_LINES);
	Vector z=axis*Math::div2(length);
	glVertex(center-z);
	glVertex(center+z);
	glEnd();
	
	glPopAttrib();
	}

void LinePrimitive::write(Misc::File& file,const Primitive::Vector& translation) const
	{
	/* Write the line parameters: */
	file.write<Scalar>((center+translation).getComponents(),3);
	file.write<Scalar>(axis.getComponents(),3);
	file.write<Scalar>(length);
	}
