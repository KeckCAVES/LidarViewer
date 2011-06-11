/***********************************************************************
LinePrimitive - Class for lines extracted from point clouds by
intersecting two plane primitives.
Copyright (c) 2008-2010 Oliver Kreylos

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
#include <IO/File.h>
#include <Comm/MulticastPipe.h>
#include <Math/Math.h>
#include <Math/Constants.h>
#include <Geometry/Box.h>
#include <Geometry/PCACalculator.h>
#include <GL/gl.h>
#include <GL/GLColorTemplates.h>
#include <GL/GLGeometryWrappers.h>

#include "LidarOctree.h"
#include "PlanePrimitive.h"

#include "LinePrimitive.h"

class LidarLineExtractor
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Point<double,3> Point; // Type for points
	typedef Geometry::Vector<double,3> Vector; // Type for vectors
	typedef Geometry::Box<double,3> Box; // Type for bounding boxes
	
	/* Elements: */
	private:
	Box bb; // Bounding box of all processed points
	Geometry::PCACalculator<3> pca; // Helper object to accumulate the points' covariance matrix and calculate their PCA
	
	/* Constructors and destructors: */
	public:
	LidarLineExtractor(void)
		:bb(Box::empty)
		{
		};
	
	/* Methods: */
	void operator()(const LidarPoint& lp) // Process the given LiDAR point
		{
		/* Add the node point to the bounding box: */
		bb.addPoint(lp);
		
		/* Add the point to the PCA calculator: */
		pca.accumulatePoint(lp);
		};
	size_t getNumPoints(void) const // Returns the number of processed points
		{
		return pca.getNumPoints();
		}
	const Box& getBB(void) const // Returns the processed points' bounding box
		{
		return bb;
		};
	void calcLine(Point& centroid,Vector& axis) // Returns the least-squares line
		{
		/* Calculate the processed points' centroid: */
		centroid=pca.calcCentroid();
		
		/* Calculate the point set's covariance matrix: */
		pca.calcCovariance();
		
		/* Calculate the covariance matrix' eigenvalues: */
		double evs[3];
		pca.calcEigenvalues(evs);
		
		/* Get the "longest" eigenvector: */
		axis=pca.calcEigenvector(evs[0]);
		};
	};

class LidarLineFitter
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Point<double,3> Point; // Type for points
	typedef Geometry::Vector<double,3> Vector; // Type for vectors
	
	/* Elements: */
	private:
	Point centroid; // Line's centroid
	Vector axis; // Line's normalized axis
	double min,max; // Bounding interval of all points in line's coordinate system
	size_t numPoints; // Number of accumulated points
	double ms; // Accumulated RMS distance from points to line
	
	/* Constructors and destructors: */
	public:
	LidarLineFitter(const Point& sCentroid,const Vector& sAxis)
		:centroid(sCentroid),axis(sAxis),
		 min(Math::Constants<double>::max),
		 max(Math::Constants<double>::min),
		 numPoints(0),ms(0.0)
		{
		/* Normalize the axis vector: */
		axis.normalize();
		};
	
	/* Methods: */
	void operator()(const LidarPoint& lp) // Process the given LiDAR point
		{
		/* Transform the point to line coordinates: */
		Vector lpc=Point(lp)-centroid;
		double x=lpc*axis;
		
		/* Add the point to the bounding interval: */
		if(min>x)
			min=x;
		if(max<x)
			max=x;
		
		/* Add the point to the RMS distance: */
		++numPoints;
		ms+=Geometry::sqr(lpc-axis*x);
		};
	double getMin(void) const
		{
		return min;
		};
	double getMax(void) const
		{
		return max;
		};
	double getRMS(void) const
		{
		return Math::sqrt(ms/double(numPoints));
		};
	};

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
			pipe->write<unsigned int>((unsigned int)(numPoints));
			pipe->write<Scalar>(rms);
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

LinePrimitive::LinePrimitive(const LidarOctree* octree,Comm::MulticastPipe* pipe)
	{
	/* Create a LiDAR line extractor: */
	LidarLineExtractor lle;
	
	/* Process all selected points: */
	octree->processSelectedPoints(lle);
	
	if(lle.getNumPoints()>=2)
		{
		/* Extract the line's coordinate frame: */
		LidarLineExtractor::Point centroid;
		LidarLineExtractor::Vector laxis;
		lle.calcLine(centroid,laxis);
		
		/* Calculate the bounding interval of the selected points in line coordinates: */
		LidarLineFitter llf(centroid,laxis);
		octree->processSelectedPoints(llf);
		double min=llf.getMin();
		double max=llf.getMax();
		double size=max-min;
		min-=0.1*size;
		max+=0.1*size;
		
		/* Store the number of points and the RMS residual: */
		numPoints=lle.getNumPoints();
		rms=llf.getRMS();
		
		/* Store the line: */
		laxis.normalize();
		center=Point(centroid+laxis*Math::mid(min,max));
		axis=Vector(laxis);
		length=Scalar(max-min);
		
		/* Print the line's equation: */
		std::cout<<"Line fitting "<<numPoints<<" points"<<std::endl;
		std::cout<<"Center point: ("<<center[0]<<", "<<center[1]<<", "<<center[2]<<")"<<std::endl;
		std::cout<<"Axis direction: ("<<axis[0]<<", "<<axis[1]<<", "<<axis[2]<<")"<<std::endl;
		std::cout<<"Length: "<<length<<std::endl;
		std::cout<<"RMS approximation residual: "<<rms<<std::endl;
		
		if(pipe!=0)
			{
			/* Send the extracted primitive over the pipe: */
			pipe->write<int>(1);
			pipe->write<unsigned int>((unsigned int)(numPoints));
			pipe->write<Scalar>(rms);
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
		Misc::throwStdErr("LinePrimitive::LinePrimitive: Not enough selected points");
		}
	}

LinePrimitive::LinePrimitive(Comm::MulticastPipe* pipe)
	{
	/* Read the status flag from the pipe: */
	if(!pipe->read<int>())
		Misc::throwStdErr("LinePrimitive::LinePrimitive: Undefined line equation");
	
	/* Read the number of points and the RMS residual: */
	numPoints=pipe->read<unsigned int>();
	rms=pipe->read<Scalar>();
	
	/* Read the line parameters: */
	pipe->read<Scalar>(center.getComponents(),3);
	pipe->read<Scalar>(axis.getComponents(),3);
	length=pipe->read<Scalar>();
	}

LinePrimitive::LinePrimitive(IO::File& file,const Primitive::Vector& translation)
	{
	/* Read the number of points and the RMS residual: */
	numPoints=file.read<unsigned int>();
	rms=file.read<Scalar>();
	
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

void LinePrimitive::write(IO::File& file,const Primitive::Vector& translation) const
	{
	/* Write the number of points and the RMS residual: */
	file.write<unsigned int>((unsigned int)(numPoints));
	file.write<Scalar>(rms);
	
	/* Write the line parameters: */
	file.write<Scalar>((center+translation).getComponents(),3);
	file.write<Scalar>(axis.getComponents(),3);
	file.write<Scalar>(length);
	}
