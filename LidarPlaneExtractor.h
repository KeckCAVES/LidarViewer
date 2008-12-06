/***********************************************************************
LidarPlaneExtractor - Point processor functor class to extract least-
squares planes from sets of selected LiDAR points.
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

#ifndef LIDARPLANEEXTRACTOR_INCLUDED
#define LIDARPLANEEXTRACTOR_INCLUDED

#include <Misc/Utility.h>
#include <Math/Math.h>
#include <Math/Constants.h>
#include <Geometry/ComponentArray.h>
#include <Geometry/Point.h>
#include <Geometry/Vector.h>
#include <Geometry/Matrix.h>
#include <Geometry/Plane.h>
#include <Geometry/Box.h>

class LidarPlaneExtractor
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Point<double,3> Point; // Type for points
	typedef Geometry::Vector<double,3> Vector; // Type for vectors
	typedef Geometry::Box<float,3> Box; // Type for bounding boxes
	private:
	typedef Geometry::ComponentArray<double,3> CA;
	typedef Geometry::Matrix<double,3,3> Matrix;
	
	/* Elements: */
	private:
	Box bb; // Bounding box of all processed points
	double pxpxs,pxpys,pxpzs,pypys,pypzs,pzpzs,pxs,pys,pzs; // Accumulated components of covariance matrix
	size_t numPoints; // Number of accumulated points
	
	/* Private methods: */
	Vector calcEigenvector(const Matrix& cov,double eigenvalue) const // Returns the eigenvector of the given matrix for the given eigenvalue
		{
		/* Create the modified covariance matrix: */
		Matrix c=cov;
		for(int i=0;i<3;++i)
			c(i,i)-=eigenvalue;
		
		/* Find the null space of the modified covariance matrix: */
		int rowIndices[3];
		for(int i=0;i<3;++i)
			rowIndices[i]=i;
		for(int step=0;step<3-1;++step)
			{
			/* Find the full pivot: */
			double pivot=Math::abs(c(step,step));
			int pivotRow=step;
			int pivotCol=step;
			for(int i=step;i<3;++i)
				for(int j=step;j<3;++j)
					{
					double val=Math::abs(c(i,j));
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
				for(int j=0;j<3;++j)
					Misc::swap(c(step,j),c(pivotRow,j));
				}
			
			/* Swap current and pivot columns if necessary: */
			if(pivotCol!=step)
				{
				/* Swap columns step and pivotCol: */
				for(int i=0;i<3;++i)
					Misc::swap(c(i,step),c(i,pivotCol));
				Misc::swap(rowIndices[step],rowIndices[pivotCol]);
				}
			
			/* Combine all rows with the current row: */
			for(int i=step+1;i<3;++i)
				{
				/* Combine rows i and step: */
				double factor=-c(i,step)/c(step,step);
				for(int j=step+1;j<3;++j)
					c(i,j)+=c(step,j)*factor;
				}
			}
		
		/* Calculate the swizzled result using backsubstitution: */
		CA x;
		x[3-1]=1.0;
		for(int i=3-2;i>=0;--i)
			{
			x[i]=0.0;
			for(int j=i+1;j<3;++j)
				x[i]-=c(i,j)*x[j];
			x[i]/=c(i,i);
			}
		
		/* Unswizzle and normalize the result: */
		Vector result;
		for(int i=0;i<3;++i)
			result[rowIndices[i]]=x[i];
		result.normalize();
		return result;
		};
	
	/* Constructors and destructors: */
	public:
	LidarPlaneExtractor(void)
		:bb(Box::empty),
		 pxpxs(0.0),pxpys(0.0),pxpzs(0.0),pypys(0.0),pypzs(0.0),pzpzs(0.0),pxs(0.0),pys(0.0),pzs(0.0),
		 numPoints(0)
		{
		};
	
	/* Methods: */
	void operator()(const LidarPoint& lp) // Process the given LiDAR point
		{
		/* Add the node point to the bounding box: */
		bb.addPoint(lp);
		
		/* Accumulate the node point: */
		pxpxs+=double(lp[0])*double(lp[0]);
		pxpys+=double(lp[0])*double(lp[1]);
		pxpzs+=double(lp[0])*double(lp[2]);
		pypys+=double(lp[1])*double(lp[1]);
		pypzs+=double(lp[1])*double(lp[2]);
		pzpzs+=double(lp[2])*double(lp[2]);
		pxs+=double(lp[0]);
		pys+=double(lp[1]);
		pzs+=double(lp[2]);
		
		++numPoints;
		};
	size_t getNumPoints(void) const // Returns the number of processed points
		{
		return numPoints;
		};
	const Box& getBB(void) const // Returns the processed points' bounding box
		{
		return bb;
		};
	void calcPlane(Point& centroid,Vector plane[3],double lengths[3]) // Returns the least-squares plane and its aligned normalized coordinate frame and the lengths of the eigenvectors
		{
		if(numPoints<3)
			return;
		double np=double(numPoints);
		
		/* Calculate the processed points' covariance matrix: */
		Matrix c;
		c(0,0)=(pxpxs-pxs*pxs/np)/np;
		c(0,1)=(pxpys-pxs*pys/np)/np;
		c(0,2)=(pxpzs-pxs*pzs/np)/np;
		c(1,0)=c(0,1);
		c(1,1)=(pypys-pys*pys/np)/np;
		c(1,2)=(pypzs-pys*pzs/np)/np;
		c(2,0)=c(0,2);
		c(2,1)=c(1,2);
		c(2,2)=(pzpzs-pzs*pzs/np)/np;
		
		/* Calculate the coefficients of the covariance matrix' characteristic polynomial: */
		double cp[3];
		cp[0]=-c(0,0)-c(1,1)-c(2,2);
		cp[1]=c(0,0)*c(1,1)+c(0,0)*c(2,2)+c(1,1)*c(2,2)-c(0,1)*c(1,0)-c(0,2)*c(2,0)-c(1,2)*c(2,1);
		cp[2]=-c(0,0)*(c(1,1)*c(2,2)-c(1,2)*c(2,1))+c(0,1)*(c(1,0)*c(2,2)-c(1,2)*c(2,0))-c(0,2)*(c(1,0)*c(2,1)-c(1,1)*c(2,0));
		
		/* Find all roots of the characteristic polynomial: */
		double roots[3];
		double q=(Math::sqr(cp[0])-3.0*cp[1])/9.0;
		double q3=Math::sqr(q)*q;
		double r=((2.0*Math::sqr(cp[0])-9.0*cp[1])*cp[0]+27.0*cp[2])/54.0;
		if(Math::sqr(r)<q3)
			{
			/* There are three real roots: */
			double theta=Math::acos(r/Math::sqrt(q3));
			roots[0]=-2.0*Math::sqrt(q)*Math::cos(theta/3.0)-cp[0]/3.0;
			roots[1]=-2.0*Math::sqrt(q)*Math::cos((theta+2.0*Math::Constants<double>::pi)/3.0)-cp[0]/3.0;
			roots[2]=-2.0*Math::sqrt(q)*Math::cos((theta-2.0*Math::Constants<double>::pi)/3.0)-cp[0]/3.0;
			}
		else
			{
			/* There is only one real root: */
			double a=Math::pow(Math::abs(r)+Math::sqrt(Math::sqr(r)-q3),1.0/3.0);
			if(r>0.0)
				a=-a;
			double b=a==0.0?0.0:q/a;
			roots[0]=a+b-cp[0]/3.0;
			roots[1]=roots[0];
			roots[2]=roots[0];
			}
		
		/* Use Newton iteration to clean up the roots: */
		for(int i=0;i<3;++i)
			for(int j=0;j<5;++j)
				{
				double f=((roots[i]+cp[0])*roots[i]+cp[1])*roots[i]+cp[2];
				double fp=(3.0*roots[i]+2.0*cp[0])*roots[i]+cp[1];
				double s=f/fp;
				roots[i]-=s;
				}
		
		/* Sort the eigenvalues by descending absolute value: */
		if(Math::abs(roots[0])<Math::abs(roots[1]))
			Misc::swap(roots[0],roots[1]);
		if(Math::abs(roots[1])<Math::abs(roots[2]))
			Misc::swap(roots[1],roots[2]);
		if(Math::abs(roots[0])<Math::abs(roots[1]))
			Misc::swap(roots[0],roots[1]);
		
		/* Return all eigenvectors and eigenvalues: */
		for(int i=0;i<3;++i)
			{
			plane[i]=calcEigenvector(c,roots[i]);
			lengths[i]=roots[i];
			}
		
		/* Calculate the processed points' centroid: */
		centroid=Point(pxs/np,pys/np,pzs/np);
		};
	};

#endif
