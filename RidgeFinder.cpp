/***********************************************************************
RidgeFinder - Functor classes to find ridges in LiDAR data sets.
Copyright (c) 2009 Oliver Kreylos

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

#include "RidgeFinder.h"

/********************************
Methods of class PlaneCalculator:
********************************/

PlaneCalculator::calcPlane(void) const
	{
	/* Calculate the covariance matrix of the point neighborhood: */
	pca.calcCovariance();
	
	/* Calculate the eigenvalues of the covariance matrix: */
	double evs[3];
	pca.calcEigenvalues(evs);
	
	/* Calculate the eigenvector of the smallest eigenvalue: */
	Plane::Vector normal=pca.calcEigenvector(evs[2]);
	
	/* Calculate the centroid of the point neighborhood: */
	Plane::Point centroid=pca.calcCentroid();
	
	/* Return the resulting plane equation: */
	return Plane(normal,centroid);
	}

/********************************
Methods of class PointClassifier:
********************************/

PointClassifier::PointClassifier(const Point& sQueryPoint,Scalar sRadius2,const PlaneCalculator::Plane& sPlane)
	:queryPoint(sQueryPoint),radius2(sRadius2),plane(sPlane)
	{
	/* Calculate the transformation into the plane's local coordinate system: */
	Vector z=plane.getNormal();
	z.normalize();
	Vector x=Geometry::normal(z);
	x.normalize();
	Vector y=Geometry::cross(z,x);
	z.normalize();
	Point o=plane.project(node[i]);
	ATransform::Matrix m;
	for(int i=0;i<3;++i)
		{
		m(i,0)=x[i];
		m(i,1)=y[i];
		m(i,2)=z[i];
		m(i,3)=o[i];
		}
	toPlane=ATransform(m);
	toPlane.doInvert();
	}

bool PointClassifier::isRidge(void) const
	{
	if(pcas[0].getNumPoints()>=2&&pcas[1].getNumPoints()>=2)
		{
		double radius=Math::sqrt(radius2);
		
		/* Calculate eigenvalues of both PCAs: */
		double evs[2][2];
		unsigned int numEvs[2];
		bool isRidges[2];
		for(int i=0;i<2;++i)
			{
			pcas[i].calcCovariance();
			numEvs[i]=pcas[i].calcEigenvalues(evs[i]);
			isRidges[i]=numEvs[i]==2&&Math::abs(evs[i][0])>=radius*0.9&&Math::abs(evs[i][1])<=radius*0.333;
			}
		
		if(isRidges[0]&&!isRidges[1])
			{
			ATransform::Point qp=toPlane.transform(queryPoint);
			return pcas[0].calcEigenvector(evs[0][1])*(PCACalculator<2>::Point(qp[0],qp[1])-pcas[0].getCentroid())<=1.0;
			}
		else if(!isRidges[0]&&isRidges[1])
			{
			ATransform::Point qp=toPlane.transform(queryPoint);
			return pcas[1].calcEigenvector(evs[1][1])*(PCACalculator<2>::Point(qp[0],qp[1])-pcas[1].getCentroid())<=1.0;
			}
		else
			return false;
		}
	else
		return false;
	}

/****************************
Methods of class RidgeFinder:
****************************/

RidgeFinder::RidgeFinder(LidarProcessOctree& sLpo,Scalar sRadius,const char* colorFileName)
	:lpo(sLpo),
	 radius2(Math::sqr(sRadius)),
	 colorBuffer(new Color[lpo.getMaxNumPointsPerNode()]),
	 colorDataSize(sizeof(Color)),
	 colorFile(colorFileName,"w+b",LidarFile::LittleEndian),
	 numProcessedNodes(0)
	{
	/* Write the color file's header: */
	LidarDataFileHeader dfh((unsigned int)(colorDataSize));
	dfh.write(colorFile);
	}

RidgeFinder::~RidgeFinder(void)
	{
	delete[] colorBuffer;
	}

void RidgeFinder::operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel)
	{
	/* Check if this node is a leaf or an interior node: */
	if(node.isLeaf())
		{
		/* Process each point in the node: */
		for(unsigned int i=0;i<node.getNumPoints();++i)
			{
			/* Calculate a best-fitting plane equation for the point's neighborhood: */
			PlaneCalculator planeCalculator(node[i],radius2);
			lpo.processPointsDirected(planeCalculator);
			
			/* Check if the point has neighbors: */
			if(planeCalculator.getNumPoints()>=3)
				{
				/* Run a ridge classifier on the point: */
				PointClassifier pc(queryPoint,radius2,planeCalculator.getPlane());
				lpo.processPointDirected(pc);
				
				if(pc.isRidge())
					colorBuffer[i]=Color(255,255,255);
				else
					colorBuffer[i]=Color(0,0,0);
				}
			else
				{
				/* Declare the point a non-ridge: */
				colorBuffer[i]=Color(0,0,0);
				}
			}
		}
	else
		{
		}
	
	/* Write the node's colors to the color file: */
	colorFile.seekSet(LidarDataFileHeader::getFileSize()+colorDataSize*node.getDataOffset());
	colorFile.write(colorBuffer,node.getNumPoints());
	
	/* Update the progress counter: */
	++numProcessedNodes;
	int percent=int((numProcessedNodes*100+lpo.getNumNodes()/2)/lpo.getNumNodes());
	std::cout<<"\b\b\b\b"<<std::setw(3)<<percent<<"%"<<std::flush;
	}
