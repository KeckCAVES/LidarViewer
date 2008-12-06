/***********************************************************************
PointPrimitive - Class for points extracted from point clouds by
intersecting three plane primitives or one line primitive and one plane
primitive.
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
#include <GL/gl.h>
#include <GL/GLColorTemplates.h>
#include <GL/GLGeometryWrappers.h>

#include "PlanePrimitive.h"
#include "LinePrimitive.h"

#include "PointPrimitive.h"

/*******************************
Methods of class PointPrimitive:
*******************************/

PointPrimitive::PointPrimitive(const PlanePrimitive* p1,const PlanePrimitive* p2,const PlanePrimitive* p3,Comm::MulticastPipe* pipe)
	{
	/* Get the three planes' plane equations: */
	const PlanePrimitive* ps[3];
	ps[0]=p1;
	ps[1]=p2;
	ps[2]=p3;
	PlanePrimitive::Plane planes[3];
	for(int i=0;i<3;++i)
		planes[i]=ps[i]->getPlane();
	
	/* Create the linear system: */
	double a[3][4];
	for(int i=0;i<3;++i)
		{
		for(int j=0;j<3;++j)
			a[i][j]=double(planes[i].getNormal()[j]);
		a[i][3]=double(planes[i].getOffset());
		}
	
	/* Solve the linear system: */
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
	double x[3];
	for(int i=3-1;i>=0;--i)
		{
		x[i]=a[i][3];
		for(int j=i+1;j<3;++j)
			x[i]-=a[i][j]*x[j];
		x[i]/=a[i][i];
		}
	
	/* Unswizzle the result: */
	for(int i=0;i<3;++i)
		point[rowIndices[i]]=Scalar(x[i]);
	
	/* Print the point's equation: */
	std::cout<<"Point intersecting three planes"<<std::endl;
	std::cout<<"Point: ("<<point[0]<<", "<<point[1]<<", "<<point[2]<<")"<<std::endl;
	
	if(pipe!=0)
		{
		/* Send the extracted primitive over the pipe: */
		pipe->write<int>(1);
		pipe->write<Scalar>(point.getComponents(),3);
		pipe->finishMessage();
		}
	}

PointPrimitive::PointPrimitive(const PlanePrimitive* p,const LinePrimitive* l,Comm::MulticastPipe* pipe)
	{
	/* Get the plane's plane equation: */
	const PlanePrimitive::Plane& plane=p->getPlane();
	
	/* Get the line's line equation: */
	const Point& center=l->getCenter();
	const Vector& axis=l->getAxis();
	
	/* Intersect the plane and the line: */
	Scalar denominator=axis*plane.getNormal();
	if(denominator!=Scalar(0))
		{
		/* Calculate the intersection point: */
		Scalar lambda=(plane.getOffset()-center*plane.getNormal())/denominator;
		point=center+axis*lambda;
		
		/* Print the point's equation: */
		std::cout<<"Point intersecting one plane and one line"<<std::endl;
		std::cout<<"Point: ("<<point[0]<<", "<<point[1]<<", "<<point[2]<<")"<<std::endl;
		
		if(pipe!=0)
			{
			/* Send the extracted primitive over the pipe: */
			pipe->write<int>(1);
			pipe->write<Scalar>(point.getComponents(),3);
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
		Misc::throwStdErr("PointPrimitive::PointPrimitive: Plane and line do not intersect");
		}
	}

PointPrimitive::PointPrimitive(Comm::MulticastPipe* pipe)
	{
	/* Read the status flag from the pipe: */
	if(!pipe->read<int>())
		Misc::throwStdErr("PointPrimitive::PointPrimitive: No valid intersection point");
	
	/* Read the point parameters: */
	pipe->read<Scalar>(point.getComponents(),3);
	}

PointPrimitive::PointPrimitive(Misc::File& file,const Primitive::Vector& translation)
	{
	/* Read the point parameters: */
	file.read<Scalar>(point.getComponents(),3);
	point+=translation;
	}

Primitive::Scalar PointPrimitive::pick(const Primitive::Point& pickPoint,Primitive::Scalar maxPickDistance) const
	{
	/* Return the distance from the pick point to the point: */
	return Geometry::dist(pickPoint,point);
	}

void PointPrimitive::glRenderAction(GLContextData& contextData) const
	{
	glPushAttrib(GL_COLOR_BUFFER_BIT|GL_ENABLE_BIT|GL_POINT_BIT);
	glDisable(GL_LIGHTING);
	
	/* Draw the point: */
	glPointSize(5.0f);
	glColor(surfaceColor);
	glBegin(GL_POINTS);
	glVertex(point);
	glEnd();
	
	glPopAttrib();
	}

void PointPrimitive::write(Misc::File& file,const Primitive::Vector& translation) const
	{
	/* Write the point parameters: */
	file.write<Scalar>((point+translation).getComponents(),3);
	}
