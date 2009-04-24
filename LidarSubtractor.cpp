/***********************************************************************
LidarSubtractor - Post-processing filter to subtract a (relatively
small) point set from a LiDAR data set and create a new LiDAR data set.
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

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <Misc/File.h>
#include <Math/Math.h>
#include <Geometry/ArrayKdTree.h>

#include "LidarTypes.h"
#include "LidarProcessOctree.h"
#include "PointAccumulator.h"
#include "LidarOctreeCreator.h"

class MatchingPointFinder
	{
	/* Elements: */
	private:
	Point p; // Searched point
	Scalar epsilon,epsilon2; // Maximum search distance
	bool found; // Flag whether a matching point was found
	
	/* Constructors and destructors: */
	public:
	MatchingPointFinder(const Point& sP,Scalar sEpsilon)
		:p(sP),epsilon(sEpsilon),epsilon2(Math::sqr(epsilon)),
		 found(false)
		{
		}
	
	/* Methods: */
	const Point& getQueryPosition(void) const
		{
		return p;
		}
	bool operator()(const Point& node,int splitDimension)
		{
		/* Compare node's point to current closest point: */
		if(Geometry::sqrDist(node,p)<epsilon2)
			found=true;
		
		/* Stop traversal if split plane is farther away than epsilon: */
		return epsilon>Math::abs(node[splitDimension]-p[splitDimension]);
		};
	bool isFound(void) const
		{
		return found;
		}
	};

class NodePointSubtractor
	{
	/* Elements: */
	private:
	LidarProcessOctree& basePoints;
	Geometry::ArrayKdTree<Point>& subtractPoints;
	Scalar epsilon;
	PointAccumulator& pa;
	
	/* Constructors and destructors: */
	public:
	NodePointSubtractor(LidarProcessOctree& sBasePoints,Geometry::ArrayKdTree<Point>& sSubtractPoints,Scalar sEpsilon,PointAccumulator& sPa)
		:basePoints(sBasePoints),subtractPoints(sSubtractPoints),
		 epsilon(sEpsilon),
		 pa(sPa)
		{
		}
	
	/* Methods: */
	void operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel)
		{
		/* Check if this node is a leaf: */
		if(node.isLeaf())
			{
			/* Look for each node's point in the subtract set: */
			for(unsigned int i=0;i<node.getNumPoints();++i)
				{
				MatchingPointFinder mpf(node[i],epsilon);
				subtractPoints.traverseTreeDirected(mpf);
				if(!mpf.isFound())
					{
					/* Point has no match in subtract set; add it to point accumulator: */
					pa.addPoint(node[i]);
					}
				}
			}
		}
	};

int main(int argc,char* argv[])
	{
	/* Parse the command line and load the input files: */
	LidarProcessOctree* basePoints=0; // Octree containing the base point data
	const char* subtractFileName=0; // Name of ASCII file containing the point set to subtract
	int asciiColumnIndices[3]; // Column indices of x, y, z point components in ASCII file
	Scalar epsilon=Scalar(1.0e-7); // Maximum match point distance
	const char* outputFileName=0; // Name of resulting LiDAR octree file
	unsigned int maxPointsPerNode=1024; // Node size of resulting LiDAR octree file
	PointAccumulator pa; // Point accumulator holding the subtracted point set
	std::string tempPointFileNameTemplate="/tmp/LidarPreprocessorTempPoints";
	
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"o")==0)
				{
				++i;
				if(i<argc)
					outputFileName=argv[i];
				else
					std::cerr<<"Dangling -o flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"np")==0)
				{
				++i;
				if(i<argc)
					maxPointsPerNode=(unsigned int)(atoi(argv[i]));
				else
					std::cerr<<"Dangling -np flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"ooc")==0)
				{
				++i;
				if(i<argc)
					pa.setMemorySize(atoi(argv[i]),pa.getMaxNumPointsPerNode());
				else
					std::cerr<<"Dangling -ooc flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"to")==0)
				{
				++i;
				if(i<argc)
					{
					std::string tempOctreeFileNameTemplate=argv[i];
					pa.setTempOctreeFileNameTemplate(tempOctreeFileNameTemplate+"XXXXXX");
					}
				else
					std::cerr<<"Dangling -to flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"tp")==0)
				{
				++i;
				if(i<argc)
					tempPointFileNameTemplate=argv[i];
				else
					std::cerr<<"Dangling -tp flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"eps")==0)
				{
				++i;
				if(i<argc)
					epsilon=Scalar(atof(argv[i]));
				else
					std::cerr<<"Dangling -eps flag on command line"<<std::endl;
				}
			}
		else if(basePoints==0)
			{
			/* Create a processing octree: */
			basePoints=new LidarProcessOctree(argv[i],size_t(64)*size_t(1024*1024));
			}
		else if(subtractFileName==0)
			{
			/* Store the subtraction file name: */
			subtractFileName=argv[i];
			for(int i=0;i<3;++i)
				asciiColumnIndices[i]=i;
			}
		else
			std::cerr<<"Ignoring command line argument "<<argv[i]<<std::endl;
		}
	
	/* Load the subtraction point set into a kd-tree: */
	std::cout<<"Loading subtraction points from "<<subtractFileName<<"..."<<std::flush;
	std::vector<Point> subtractPoints;
	Misc::File subtractFile(subtractFileName,"rt");
	while(!subtractFile.eof())
		{
		char line[256];
		subtractFile.gets(line,sizeof(line));
		float p[3];
		if(sscanf(line,"%f %f %f",&p[0],&p[1],&p[2])==3)
			subtractPoints.push_back(Point(p));
		}
	std::cout<<" done"<<std::endl;
	
	std::cout<<"Creating kd-tree of "<<subtractPoints.size()<<" subtraction points..."<<std::flush;
	Geometry::ArrayKdTree<Point>* subtractPointTree=new Geometry::ArrayKdTree<Point>(subtractPoints.size());
	Point* points=subtractPointTree->accessPoints();
	for(size_t i=0;i<subtractPoints.size();++i)
		points[i]=subtractPoints[i];
	subtractPointTree->releasePoints(4);
	std::cout<<" done"<<std::endl;
	
	/* Process the base point set: */
	{
	std::cout<<"Subtracting points..."<<std::flush;
	NodePointSubtractor nps(*basePoints,*subtractPointTree,epsilon,pa);
	basePoints->processNodesPostfix(nps);
	pa.finishReading();
	std::cout<<" done"<<std::endl;
	}
	
	/* Clear input data structures: */
	delete basePoints;
	delete subtractPointTree;
	
	/* Construct an octree with less than maxPointsPerNode points per leaf: */
	LidarOctreeCreator tree(pa.getMaxNumCacheablePoints(),maxPointsPerNode,pa.getTempOctrees(),tempPointFileNameTemplate+"XXXXXX");
	
	/* Delete the temporary point octrees: */
	pa.deleteTempOctrees();
	
	/* Write the octree structure and data to the destination LiDAR file: */
	tree.write(outputFileName);
	
	return 0;
	}
