/***********************************************************************
NormalCalculator - Functor class to calculate a normal vector for a
point in a LiDAR data set.
Copyright (c) 2008-2014 Oliver Kreylos

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
#include <iomanip>
#include <Misc/Utility.h>
#include <Misc/ThrowStdErr.h>
#include <Math/Math.h>

#include "NormalCalculator.h"

/*********************************
Methods of class NormalCalculator:
*********************************/

NormalCalculator::Plane::Vector NormalCalculator::calcEigenvector(const NormalCalculator::Matrix& cov,double eigenvalue)
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
	Plane::Vector result;
	for(int i=0;i<3;++i)
		result[rowIndices[i]]=x[i];
	result.normalize();
	return result;
	}

NormalCalculator::Plane::Vector NormalCalculator::calcNormal(const NormalCalculator::Matrix& cov)
	{
	/* Calculate the coefficients of the covariance matrix' characteristic polynomial: */
	double cp[3];
	cp[0]=-cov(0,0)-cov(1,1)-cov(2,2);
	cp[1]=cov(0,0)*cov(1,1)+cov(0,0)*cov(2,2)+cov(1,1)*cov(2,2)-cov(0,1)*cov(1,0)-cov(0,2)*cov(2,0)-cov(1,2)*cov(2,1);
	cp[2]=-cov(0,0)*(cov(1,1)*cov(2,2)-cov(1,2)*cov(2,1))+cov(0,1)*(cov(1,0)*cov(2,2)-cov(1,2)*cov(2,0))-cov(0,2)*(cov(1,0)*cov(2,1)-cov(1,1)*cov(2,0));
	
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
	
	/* Calculate the smallest eigenvector: */
	return calcEigenvector(cov,roots[2]);
	}

/***************************************
Methods of class RadiusNormalCalculator:
***************************************/

RadiusNormalCalculator::RadiusNormalCalculator(Scalar sRadius)
	:radius2(Math::sqr(sRadius))
	{
	}

void RadiusNormalCalculator::prepare(const Point& newQueryPoint)
	{
	/* Copy the query point: */
	queryPoint=newQueryPoint;
	
	/* Reset the PCA accumulator: */
	pxpxs=0.0;
	pxpys=0.0;
	pxpzs=0.0;
	pypys=0.0;
	pypzs=0.0;
	pzpzs=0.0;
	pxs=0.0;
	pys=0.0;
	pzs=0.0;
	numPoints=0;
	
	/* Reset the closest distance: */
	closestDist2=radius2;
	}

NormalCalculator::Plane RadiusNormalCalculator::calcPlane(void) const
	{
	if(numPoints<3)
		Misc::throwStdErr("RadiusNormalCalculator::calcPlane: Too few processed points, have %u instead of 3",numPoints);
	
	/* Calculate the processed points' covariance matrix: */
	double np=double(numPoints);
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
	
	/* Return the plane equation: */
	return Plane(calcNormal(c),Plane::Point(pxs/np,pys/np,pzs/np));
	}

/*********************************************
Methods of class NumberRadiusNormalCalculator:
*********************************************/

NumberRadiusNormalCalculator::NumberRadiusNormalCalculator(unsigned int sMaxNumNeighbors)
	:maxNumNeighbors(sMaxNumNeighbors),maxDist2(Math::Constants<Scalar>::max),
	 neighbors(new Neighbor[maxNumNeighbors])
	{
	}

NumberRadiusNormalCalculator::NumberRadiusNormalCalculator(unsigned int sMaxNumNeighbors,Scalar sMaxDist)
	:maxNumNeighbors(sMaxNumNeighbors),maxDist2(Math::sqr(sMaxDist)),
	 neighbors(new Neighbor[maxNumNeighbors])
	{
	}

NumberRadiusNormalCalculator::NumberRadiusNormalCalculator(const NumberRadiusNormalCalculator& source)
	:maxNumNeighbors(source.maxNumNeighbors),maxDist2(source.maxDist2),
	 neighbors(new Neighbor[maxNumNeighbors]),
	 currentNumNeighbors(source.currentNumNeighbors),currentMaxDist2(source.currentMaxDist2)
	{
	/* Copy the source's neighbor heap: */
	for(unsigned int i=0;i<maxNumNeighbors;++i)
		neighbors[i]=source.neighbors[i];
	}

NumberRadiusNormalCalculator& NumberRadiusNormalCalculator::operator=(const NumberRadiusNormalCalculator& source)
	{
	if(this!=&source)
		{
		if(maxNumNeighbors!=source.maxNumNeighbors)
			{
			maxNumNeighbors=source.maxNumNeighbors;
			delete[] neighbors;
			neighbors=new Neighbor[maxNumNeighbors];
			}
		maxDist2=source.maxDist2;
		currentNumNeighbors=source.currentNumNeighbors;
		currentMaxDist2=source.currentMaxDist2;
		
		/* Copy the source's neighbor heap: */
		for(unsigned int i=0;i<maxNumNeighbors;++i)
			neighbors[i]=source.neighbors[i];
		}
	return *this;
	}

NumberRadiusNormalCalculator::~NumberRadiusNormalCalculator(void)
	{
	delete[] neighbors;
	}

void NumberRadiusNormalCalculator::operator()(const LidarPoint& point)
	{
	Scalar dist2=Geometry::sqrDist(point,queryPoint);
	if(dist2<currentMaxDist2)
		{
		if(currentNumNeighbors<maxNumNeighbors)
			{
			/* Insert the new point into the heap: */
			unsigned int insertionPos=currentNumNeighbors;
			while(insertionPos>0)
				{
				unsigned int parent=(insertionPos-1)>>1;
				if(neighbors[parent].dist2>=dist2)
					break;
				neighbors[insertionPos]=neighbors[parent];
				insertionPos=parent;
				}
			neighbors[insertionPos].point=point;
			neighbors[insertionPos].dist2=dist2;
			
			/* Increment the current neighborhood size and check if it became full: */
			++currentNumNeighbors;
			if(currentNumNeighbors==maxNumNeighbors)
				currentMaxDist2=neighbors[0].dist2;
			}
		else
			{
			/* Replace the currently farthest-away neighbor in the heap: */
			unsigned int insertionPos=0;
			while(true)
				{
				unsigned int biggestIndex=insertionPos;
				Scalar biggest=dist2;
				unsigned int child=(insertionPos<<1);
				for(int i=0;i<2;++i)
					{
					++child;
					if(child<maxNumNeighbors&&neighbors[child].dist2>biggest)
						{
						biggestIndex=child;
						biggest=neighbors[child].dist2;
						}
					}
				if(biggestIndex==insertionPos)
					break;
				neighbors[insertionPos]=neighbors[biggestIndex];
				insertionPos=biggestIndex;
				}
			neighbors[insertionPos].point=point;
			neighbors[insertionPos].dist2=dist2;
			
			/* Update the current neighborhood radius: */
			currentMaxDist2=neighbors[0].dist2;
			}
		}
	}

void NumberRadiusNormalCalculator::prepare(const Point& newQueryPoint)
	{
	/* Copy the query point: */
	queryPoint=newQueryPoint;
	
	/* Reset the neighborhood collector: */
	currentNumNeighbors=0;
	currentMaxDist2=maxDist2;
	}

NormalCalculator::Plane NumberRadiusNormalCalculator::calcPlane(void) const
	{
	if(currentNumNeighbors<3)
		Misc::throwStdErr("NumberRadiusNormalCalculator::calcPlane: Too few processed points, have %u instead of 3",currentNumNeighbors);
	
	/* Calculate the processed points' covariance matrix: */
	double pxpxs=0.0;
	double pxpys=0.0;
	double pxpzs=0.0;
	double pypys=0.0;
	double pypzs=0.0;
	double pzpzs=0.0;
	double pxs=0.0;
	double pys=0.0;
	double pzs=0.0;
	for(unsigned int i=0;i<currentNumNeighbors;++i)
		{
		/* Accumulate the neighbor: */
		const Point& lp=neighbors[i].point;
		pxpxs+=double(lp[0])*double(lp[0]);
		pxpys+=double(lp[0])*double(lp[1]);
		pxpzs+=double(lp[0])*double(lp[2]);
		pypys+=double(lp[1])*double(lp[1]);
		pypzs+=double(lp[1])*double(lp[2]);
		pzpzs+=double(lp[2])*double(lp[2]);
		pxs+=double(lp[0]);
		pys+=double(lp[1]);
		pzs+=double(lp[2]);
		}
	
	double np=double(currentNumNeighbors);
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
	
	/* Return the plane equation: */
	return Plane(calcNormal(c),Plane::Point(pxs/np,pys/np,pzs/np));
	}

Scalar NumberRadiusNormalCalculator::getClosestDist(void) const
	{
	/* Find the closest non-identical neighbor: */
	Scalar result2=maxDist2;
	for(unsigned int i=0;i<currentNumNeighbors;++i)
		{
		if(neighbors[i].dist2>Scalar(0)&&result2>neighbors[i].dist2)
			result2=neighbors[i].dist2;
		}
	
	return Math::sqrt(result2);
	}
