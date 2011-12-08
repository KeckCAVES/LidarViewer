/***********************************************************************
PointSetSimilarity - Tool to calculate a similarity measure between two
(relatively small) subsets of the same original point set.
Copyright (c) 2009-2011 Oliver Kreylos

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
#include <iostream>
#include <IO/OpenFile.h>
#include <IO/ValueSource.h>
#include <Math/Math.h>
#include <Geometry/ArrayKdTree.h>

#include "LidarTypes.h"

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

typedef Geometry::ArrayKdTree<Point> PointTree;

PointTree* loadPointSet(const char* pointFileName)
	{
	std::vector<Point> points;
	{
	std::cout<<"Loading points from "<<pointFileName<<"..."<<std::flush;
	IO::ValueSource pointSource(IO::openFile(pointFileName));
	pointSource.setWhitespace(',',true);
	pointSource.setPunctuation('\n',true);
	pointSource.skipWs();
	while(!pointSource.eof())
		{
		/* Read the next point: */
		Point p;
		for(int i=0;i<3;++i)
			p[i]=Scalar(pointSource.readNumber());
		points.push_back(p);
		
		/* Skip the rest of the line: */
		pointSource.skipLine();
		pointSource.skipWs();
		}
	std::cout<<" done"<<std::endl;
	}
	
	std::cout<<"Creating kd-tree of "<<points.size()<<" points..."<<std::flush;
	PointTree* pointTree=new PointTree(points.size());
	Point* ptps=pointTree->accessPoints();
	for(size_t i=0;i<points.size();++i)
		ptps[i]=points[i];
	pointTree->releasePoints(4);
	std::cout<<" done"<<std::endl;
	
	return pointTree;
	}

int main(int argc,char* argv[])
	{
	/* Load the two point sets: */
	PointTree* trees[2];
	for(int i=0;i<2;++i)
		trees[i]=loadPointSet(argv[1+i]);
	
	/* Get the epsilon value: */
	Scalar epsilon=Scalar(atof(argv[3]));
	
	/* Compare the point sets: */
	int numMatches=0;
	for(int tree=0;tree<2;++tree)
		{
		PointTree& tree1=*trees[tree];
		PointTree& tree2=*trees[1-tree];
		for(int i=0;i<tree1.getNumNodes();++i)
			{
			/* Find a match for the point in tree1 in tree2: */
			MatchingPointFinder mpf(tree1.getNode(i),epsilon);
			tree2.traverseTreeDirected(mpf);
			if(mpf.isFound())
				++numMatches;
			}
		}
	
	/* Print the result: */
	std::cout<<"Similarity percentage between the two point sets: "<<double(numMatches)*100.0/double(trees[0]->getNumNodes()+trees[1]->getNumNodes())<<"%"<<std::endl;
	
	for(int i=0;i<2;++i)
		delete trees[i];
	
	return 0;
	}
