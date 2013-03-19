/***********************************************************************
SubtractorHelper - Helper functions and classes for point subtraction
utilities.
Copyright (c) 2009-2013 Oliver Kreylos

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

#include "SubtractorHelper.h"

#include <stdexcept>
#include <vector>
#include <iostream>
#include <Misc/FileNameExtensions.h>
#include <IO/File.h>
#include <IO/OpenFile.h>
#include <IO/ReadAheadFilter.h>
#include <IO/ValueSource.h>
#include <Geometry/Vector.h>

Geometry::ArrayKdTree<Point>* loadSubtractSet(const char* fileName,const Geometry::Vector<double,3>& offset,int numThreads)
	{
	Geometry::ArrayKdTree<Point>* subtractPointTree=0;
	if(Misc::hasCaseExtension(fileName,".bin"))
		{
		std::cout<<"Loading subtraction points from binary file "<<fileName<<"..."<<std::flush;
		
		try
			{
			/* Load a point set from a binary file in octree coordinates: */
			IO::FilePtr pointFile(IO::openFile(fileName));
			pointFile->setEndianness(Misc::LittleEndian);
			unsigned int numPoints=pointFile->read<unsigned int>();
			
			subtractPointTree=new Geometry::ArrayKdTree<Point>(numPoints);
			Point* pPtr=subtractPointTree->accessPoints();
			if(offset!=Geometry::Vector<double,3>::zero)
				{
				Geometry::Vector<float,3> o=offset;
				for(unsigned int i=0;i<numPoints;++i,++pPtr)
					{
					/* Read a LiDAR point, apply an offset, and store a regular point: */
					LidarPoint lp=pointFile->read<LidarPoint>();
					lp-=o;
					*pPtr=lp;
					}
				}
			else
				{
				for(unsigned int i=0;i<numPoints;++i,++pPtr)
					{
					/* Read a LiDAR point and store a regular point: */
					LidarPoint lp=pointFile->read<LidarPoint>();
					*pPtr=lp;
					}
				}
			std::cout<<" done"<<std::endl;
			
			/* Create the kd-tree: */
			std::cout<<"Creating kd-tree of "<<numPoints<<" subtraction points..."<<std::flush;
			if(numThreads>1)
				subtractPointTree->releasePoints(numThreads);
			else
				subtractPointTree->releasePoints();
			std::cout<<" done"<<std::endl;
			}
		catch(std::runtime_error err)
			{
			std::cout<<" failed"<<std::endl;
			std::cerr<<"Caught exception "<<err.what()<<"while reading binary subtraction file "<<fileName<<"; terminating"<<std::endl;
			delete subtractPointTree;
			return 0;
			}
		}
	else
		{
		std::cout<<"Loading subtraction points from ASCII file "<<fileName<<"..."<<std::flush;
		
		try
			{
			IO::ValueSource subtractSource(new IO::ReadAheadFilter(IO::openFile(fileName)));
			subtractSource.setWhitespace(',',true);
			subtractSource.setPunctuation('\n',true);
			subtractSource.skipWs();
			std::vector<Point> subtractPoints;
			while(!subtractSource.eof())
				{
				/* Read the next point and subtract the base file's point offset: */
				Point p;
				for(int i=0;i<3;++i)
					p[i]=Scalar(subtractSource.readNumber()-offset[i]);
				subtractPoints.push_back(p);
				
				/* Skip the rest of the line: */
				subtractSource.skipLine();
				subtractSource.skipWs();
				}
			std::cout<<" done"<<std::endl;
			
			std::cout<<"Creating kd-tree of "<<subtractPoints.size()<<" subtraction points..."<<std::flush;
			subtractPointTree=new Geometry::ArrayKdTree<Point>(subtractPoints.size());
			Point* points=subtractPointTree->accessPoints();
			for(size_t i=0;i<subtractPoints.size();++i)
				points[i]=subtractPoints[i];
			if(numThreads>1)
				subtractPointTree->releasePoints(numThreads);
			else
				subtractPointTree->releasePoints();
			std::cout<<" done"<<std::endl;
			}
		catch(std::runtime_error err)
			{
			std::cout<<" failed"<<std::endl;
			std::cerr<<"Caught exception "<<err.what()<<"while reading ASCII subtraction file "<<fileName<<"; terminating"<<std::endl;
			delete subtractPointTree;
			return 0;
			}
		}
	
	return subtractPointTree;
	}
