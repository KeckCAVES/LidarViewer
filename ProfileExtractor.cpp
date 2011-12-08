/***********************************************************************
ProfileExtractor - Algorithm to extract straight-line profiles from 2.5D
LiDAR data.
Copyright (c) 2010-2011 Oliver Kreylos

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

#include "ProfileExtractor.h"

#include <vector>
#include <Cluster/MulticastPipe.h>
#include <SceneGraph/GroupNode.h>
#include <SceneGraph/ColorNode.h>
#include <SceneGraph/CoordinateNode.h>
#include <SceneGraph/IndexedLineSetNode.h>
#include <SceneGraph/ShapeNode.h>

#include "LidarOctree.h"
#include "SceneGraph.h"

class Sampler // Class to sample the implicit surface defined by a LiDAR point cloud at an (x, y) position
	{
	/* Elements: */
	private:
	double pos[2]; // Sample point's (x, y) position
	double dx[2],dy[2]; // Sample frame orientations in x and y
	double step[2]; // Sample step size along dx and dy orientations
	double numLobes; // Number of lobes in Lanczos reconstruction filter
	double accumulator; // Accumulated convolution between point cloud and filter
	double weightSum; // Sum of weights for all accumulated points
	double absWeightSum; // Sum of absolute weights for all accumulated points
	
	/* Constructors and destructors: */
	public:
	Sampler(const double sDx[2],const double sDy[2],const double sStep[2],int sNumLobes) // Creates an empty sampler
		:numLobes(double(sNumLobes))
		{
		for(int i=0;i<2;++i)
			{
			dx[i]=sDx[i];
			dy[i]=sDy[i];
			step[i]=sStep[i];
			}
		}
	
	/* Methods: */
	void setPos(const Point& p)
		{
		/* Copy the sample position: */
		for(int i=0;i<2;++i)
			pos[i]=double(p[i]);
		
		/* Reset the filter accumulator: */
		accumulator=0.0;
		weightSum=0.0;
		absWeightSum=0.0;
		}
	void operator()(const LidarPoint& p)
		{
		/* Calculate the sample coordinates of the LiDAR point: */
		double lp[2];
		lp[0]=double(p[0]-pos[0])*dx[0]+double(p[1]-pos[1])*dx[1];
		lp[1]=double(p[0]-pos[0])*dy[0]+double(p[1]-pos[1])*dy[1];
		
		/* Calculate the Lanczos filter weight for the LiDAR point: */
		double weight=1.0;
		for(int i=0;i<2;++i)
			{
			double x=Math::Constants<double>::pi*lp[i]/step[i];
			if(Math::abs(x)>=numLobes)
				weight=0.0;
			else if(x!=0.0)
				{
				weight*=Math::sin(x)/x;
				x/=numLobes;
				weight*=Math::sin(x)/x;
				}
			}
		
		/* Accumulate the weighted point: */
		accumulator+=double(p[2])*weight;
		weightSum+=weight;
		absWeightSum+=Math::abs(weight);
		}
	double getWeightSum(void) const // Returns the sum of point weights
		{
		return weightSum;
		}
	double getAbsWeightSum(void) const // Returns the sum of absolute point weights
		{
		return absWeightSum;
		}
	double getValue(void) const // Returns the convolution result
		{
		/* Return the weighted average: */
		return accumulator/weightSum;
		}
	};

void extractProfile(const LidarOctree* octree,const Point& p0,const Point& p1,double segmentLength,int oversampling,double filterWidth,int numLobes,Cluster::MulticastPipe* pipe)
	{
	std::vector<Point> profilePoints;
	
	if(pipe==0||pipe->isMaster())
		{
		/* Calculate the distribution of sample points along the profile: */
		double pLen=Geometry::dist(p0,p1);
		int numSegments=int(pLen/segmentLength+0.5);
		segmentLength=pLen/double(numSegments);
		numSegments*=oversampling;
		
		/* Create the profile point list: */
		profilePoints.reserve(numSegments+1);
		profilePoints.push_back(p0);
		for(int i=1;i<numSegments;++i)
			{
			float lambda=float(double(i)/double(numSegments));
			profilePoints.push_back(Geometry::affineCombination(p0,p1,lambda));
			}
		profilePoints.push_back(p1);
		
		/* Set up the sampling filter frame: */
		double dx[2],dy[2],step[2];
		dx[0]=double(p1[0])-double(p0[0]);
		dx[1]=double(p1[1])-double(p0[1]);
		double dxMag=Math::sqrt(Math::sqr(dx[0])+Math::sqr(dx[1]));
		dx[0]/=dxMag;
		dx[1]/=dxMag;
		dy[0]=dx[1];
		dy[1]=-dx[0];
		step[0]=segmentLength;
		step[1]=filterWidth;
		Sampler sampler(dx,dy,step,numLobes);
		
		/* Sample all profile points: */
		for(std::vector<Point>::iterator ppIt=profilePoints.begin();ppIt!=profilePoints.end();++ppIt)
			{
			/* Calculate the bounding box of the sampling filter's support: */
			Vector bx=Vector(dx[0],dx[1],0.0f)*(float(segmentLength)*float(numLobes));
			Vector by=Vector(dy[0],dy[1],0.0f)*(float(filterWidth)*float(numLobes));
			Box frame=Box::empty;
			frame.min[2]=Math::Constants<float>::min;
			frame.max[2]=Math::Constants<float>::max;
			frame.addPoint(*ppIt-bx-by);
			frame.addPoint(*ppIt-bx+by);
			frame.addPoint(*ppIt+bx-by);
			frame.addPoint(*ppIt+bx+by);
			
			/* Sample the point: */
			sampler.setPos(*ppIt);
			octree->processPointsInBox(frame,sampler);
			
			/* Set the profile point's elevation: */
			(*ppIt)[2]=sampler.getValue();
			}
		
		if(pipe!=0)
			{
			/* Send the extracted profile points to the slaves: */
			pipe->write<unsigned int>(profilePoints.size());
			for(std::vector<Point>::iterator ppIt=profilePoints.begin();ppIt!=profilePoints.end();++ppIt)
				pipe->write<Point::Scalar>(ppIt->getComponents(),3);
			pipe->flush();
			}
		}
	else
		{
		/* Retrieve the extracted profile points from the master: */
		unsigned int numProfilePoints=pipe->read<unsigned int>();
		profilePoints.reserve(numProfilePoints);
		for(unsigned int i=0;i<numProfilePoints;++i)
			{
			Point pp;
			pipe->read<Point::Scalar>(pp.getComponents(),3);
			profilePoints.push_back(pp);
			}
		}
	
	/* Add the profile to the scene graph root: */
	SceneGraph::GroupNode* root=new SceneGraph::GroupNode;
	getSceneGraphRoot().children.appendValue(root);
	getSceneGraphRoot().update();
	
	SceneGraph::ShapeNode* s=new SceneGraph::ShapeNode;
	root->children.appendValue(s);
	{
	SceneGraph::IndexedLineSetNode* ils=new SceneGraph::IndexedLineSetNode;
	s->geometry.setValue(ils);
	
	SceneGraph::ColorNode* color=new SceneGraph::ColorNode;
	ils->color.setValue(color);
	color->color.appendValue(SceneGraph::Color(0.0f,0.5f,0.0f));
	color->update();
	
	SceneGraph::CoordinateNode* coord=new SceneGraph::CoordinateNode;
	ils->coord.setValue(coord);
	for(std::vector<Point>::iterator ppIt=profilePoints.begin();ppIt!=profilePoints.end();++ppIt)
		coord->point.appendValue(*ppIt);
	coord->update();
	
	for(unsigned int i=0;i<profilePoints.size();++i)
		ils->coordIndex.appendValue(i);
	
	ils->colorPerVertex.setValue(false);
	ils->lineWidth.setValue(3.0f);
	ils->update();
	}
	s->update();
	
	root->update();
	}
