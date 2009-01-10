/***********************************************************************
NormalCalculator - Functor class to calculate a normal vector for a
point in a LiDAR data set.
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
#include <iomanip>
#include <Misc/Utility.h>
#include <Misc/ThrowStdErr.h>
#include <Math/Math.h>

#include "NormalCalculator.h"

/*********************************
Methods of class NormalCalculator:
*********************************/

NormalCalculator::Plane::Vector NormalCalculator::calcEigenvector(const NormalCalculator::Matrix& cov,double eigenvalue) const
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

NormalCalculator::NormalCalculator(const Point& sQueryPoint,Scalar sRadius2)
	:queryPoint(sQueryPoint),
	 radius2(sRadius2),
	 pxpxs(0.0),pxpys(0.0),pxpzs(0.0),pypys(0.0),pypzs(0.0),pzpzs(0.0),pxs(0.0),pys(0.0),pzs(0.0),
	 numPoints(0)
	{
	}

NormalCalculator::Plane NormalCalculator::calcPlane(void) const
	{
	if(numPoints<3)
		Misc::throwStdErr("PlaneFitter::calcPlane: Too few processed points, have %u instead of 3",numPoints);
	
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
	
	/* Calculate the smallest eigenvector: */
	Plane::Vector normal=calcEigenvector(c,roots[2]);
	
	/* Calculate the processed points' centroid: */
	Plane::Point centroid=Plane::Point(pxs/np,pys/np,pzs/np);
	
	/* Return the plane equation: */
	return Plane(normal,centroid);
	}

namespace {

/**************
Helper classes:
**************/

class FindPoint // Class to find a point inside an octree node
	{
	/* Elements: */
	private:
	Point queryPoint; // The position of the point to find
	const LidarPoint* foundPoint; // The found LiDAR point
	
	/* Constructors and destructors: */
	public:
	FindPoint(const Point& sQueryPoint)
		:queryPoint(sQueryPoint),
		 foundPoint(0)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& lp)
		{
		if(Geometry::sqrDist(lp,queryPoint)==Scalar(0))
			foundPoint=&lp;
		}
	const Point& getQueryPoint(void) const
		{
		return queryPoint;
		}
	Scalar getQueryRadius2(void) const
		{
		return Scalar(0);
		}
	const LidarPoint* getFoundPoint(void) const
		{
		return foundPoint;
		}
	};

class NormalAverager // Class to average normal vectors of collapsed points during subsampling
	{
	/* Elements: */
	private:
	Point queryPoint; // The LiDAR point whose neighbors to find
	Scalar radius2; // Squared radius of search sphere around query point
	const Vector& pointNormal; // Normal vector onto which to project near normals before averaging
	const LidarPoint* pointBase; // Base pointer of processed node's point array
	const Vector* normalBase; // Base pointer of processed node's normal vector array
	Vector normal; // The averaged normal vector
	
	/* Constructors and destructors: */
	public:
	NormalAverager(const Point& sQueryPoint,Scalar sRadius2,const Vector& sPointNormal)
		:queryPoint(sQueryPoint),
		 radius2(sRadius2),
		 pointNormal(sPointNormal),
		 pointBase(0),normalBase(0),
		 normal(Vector::zero)
		{
		}
	
	/* Methods: */
	void setArrays(const LidarPoint* sPointBase,const Vector* sNormalBase) // Sets the point and normal vector arrays for the next process
		{
		pointBase=sPointBase;
		normalBase=sNormalBase;
		}
	void operator()(const LidarPoint& lp)
		{
		/* Check if the point is inside the search radius: */
		if(Geometry::sqrDist(lp,queryPoint)<=radius2)
			{
			/* Get the point's normal vector by mapping its index into the normal array: */
			const Vector& nearNormal=normalBase[&lp-pointBase];
			
			/* Accumulate the result normal: */
			if(nearNormal*pointNormal>=Scalar(0))
				normal+=nearNormal;
			else
				normal-=nearNormal;
			}
		}
	const Point& getQueryPoint(void) const
		{
		return queryPoint;
		}
	Scalar getQueryRadius2(void) const
		{
		return radius2;
		}
	const Vector& getNormal(void) const // Returns the averaged normal vector
		{
		return normal;
		}
	};

}

/*********************************
Methods of class NormalCalculator:
*********************************/

NodeNormalCalculator::NodeNormalCalculator(LidarProcessOctree& sLpo,Scalar sRadius,const char* normalFileName)
	:lpo(sLpo),
	 radius2(Math::sqr(sRadius)),
	 normalBuffer(new Vector[lpo.getMaxNumPointsPerNode()]),
	 normalDataSize(sizeof(Scalar)*3),
	 normalFile(normalFileName,"w+b",LidarFile::LittleEndian),
	 numProcessedNodes(0)
	{
	/* Create the child normal buffers: */
	for(int i=0;i<8;++i)
		childNormalBuffers[i]=new Vector[lpo.getMaxNumPointsPerNode()];
	
	/* Write the normal file's header: */
	LidarDataFileHeader dfh((unsigned int)(normalDataSize));
	dfh.write(normalFile);
	}

NodeNormalCalculator::~NodeNormalCalculator(void)
	{
	/* Delete all buffers: */
	delete[] normalBuffer;
	for(int i=0;i<8;++i)
		delete[] childNormalBuffers[i];
	}

void NodeNormalCalculator::operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel)
	{
	/* Check if this node is a leaf or an interior node: */
	if(node.isLeaf())
		{
		/* Calculate a normal vector for each LiDAR point in this node: */
		for(unsigned int i=0;i<node.getNumPoints();++i)
			{
			/* Create a plane fitter for the point: */
			NormalCalculator normalCalculator(node[i],radius2);
			
			/* Process the point's neighborhood: */
			lpo.processPointsDirected(normalCalculator);
			
			/* Get the point's normal vector: */
			if(normalCalculator.getNumPoints()>=3)
				{
				/* Get the fitted plane's normal vector: */
				normalBuffer[i]=normalCalculator.calcPlane().getNormal();
				
				#if 0
				/* Flip the normal vector so it always points "upwards" (hacky): */
				if(normalBuffer[i][2]<Scalar(0))
					normalBuffer[i]=-normalBuffer[i];
				#endif
				}
			else
				{
				/* Solitary point; assign dummy normal vector: */
				normalBuffer[i]=Vector::zero;
				}
			}
		}
	else
		{
		/*****************************
		Subsample the node's children:
		*****************************/
		
		/* Get pointers to the node's children and load their normal vector arrays: */
		LidarProcessOctree::Node* children[8];
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			children[childIndex]=lpo.getChild(&node,childIndex);
			if(children[childIndex]->getNumPoints()>0)
				{
				normalFile.seekSet(LidarDataFileHeader::getFileSize()+normalDataSize*children[childIndex]->getDataOffset());
				normalFile.read(childNormalBuffers[childIndex],children[childIndex]->getNumPoints());
				}
			}
		
		/* Find the points that were collapsed onto each node point and average their normal vectors: */
		Scalar averageRadius2=Math::sqr(node.getDetailSize()*Scalar(1.5));
		for(unsigned int i=0;i<node.getNumPoints();++i)
			{
			/*****************************************************************
			Find the exact normal vector of this point's "ancestor" to
			properly average normals:
			*****************************************************************/
			
			/* Find the child node containing this point's ancestor: */
			int pointChildIndex=node.getDomain().findChild(node[i]);
			LidarProcessOctree::Node* pointChild=children[pointChildIndex];
			
			/* Find the point's ancestor: */
			FindPoint fp(node[i]);
			lpo.processNodePointsDirected(pointChild,fp);
			if(fp.getFoundPoint()==0)
				Misc::throwStdErr("Things are fucked up!");
			
			/* Retrieve the ancestor's normal vector: */
			const Vector& pointNormal=childNormalBuffers[pointChildIndex][fp.getFoundPoint()-pointChild->getPoints()];
			
			/* Create a functor to average normal vectors from the point's neighborhood: */
			NormalAverager normalAverager(node[i],averageRadius2,pointNormal);
			for(int childIndex=0;childIndex<8;++childIndex)
				{
				/* Check if the child node's domain overlaps the search sphere: */
				if(children[childIndex]->getDomain().sqrDist(node[i])<=averageRadius2)
					{
					/* Search for neighbors in this child node: */
					normalAverager.setArrays(children[childIndex]->getPoints(),childNormalBuffers[childIndex]);
					lpo.processNodePointsDirected(children[childIndex],normalAverager);
					}
				}
			
			/* Store the averaged normal vector: */
			normalBuffer[i]=normalAverager.getNormal();
			if(normalBuffer[i]!=Vector::zero)
				normalBuffer[i].normalize();
			}
		}
	
	/* Write the node's normal vectors to the normal file: */
	normalFile.seekSet(LidarDataFileHeader::getFileSize()+normalDataSize*node.getDataOffset());
	normalFile.write(normalBuffer,node.getNumPoints());
	
	/* Update the progress counter: */
	++numProcessedNodes;
	int percent=int((numProcessedNodes*100+lpo.getNumNodes()/2)/lpo.getNumNodes());
	std::cout<<"\b\b\b\b"<<std::setw(3)<<percent<<"%"<<std::flush;
	}
