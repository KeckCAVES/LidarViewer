/***********************************************************************
LidarSpotRemover - Utility to extract isolated points (spot noise) from
LiDAR data sets.
Copyright (c) 2013 Oliver Kreylos

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
#include <string.h>
#include <iostream>
#include <IO/SeekableFile.h>
#include <IO/OpenFile.h>

#include "LidarTypes.h"
#include "LidarProcessOctree.h"

/**************
Helper classes:
**************/

class NeighborCounter
	{
	/* Elements: */
	private:
	Point center;
	float radius2;
	unsigned int minNumNeighbors;
	unsigned int numNeighbors;
	
	/* Constructors and destructors: */
	public:
	NeighborCounter(const Point& sCenter,float sRadius2,unsigned int sMinNumNeighbors)
		:center(sCenter),radius2(sRadius2),
		 minNumNeighbors(sMinNumNeighbors),
		 numNeighbors(0)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& point)
		{
		/* Check if the point is inside the neighborhood sphere; if so, increase neighbor count: */
		if(Geometry::sqrDist(center,point)<=radius2)
			{
			/* Stop traversing the tree if the required number of neighbors have been found: */
			if(++numNeighbors==minNumNeighbors)
				radius2=0.0f;
			}
		}
	const Point& getQueryPoint(void) const
		{
		return center;
		}
	Scalar getQueryRadius2(void) const
		{
		return radius2;
		}
	unsigned int getNumNeighbors(void) const
		{
		return numNeighbors;
		}
	};

class OutlierSaver
	{
	/* Elements: */
	private:
	LidarProcessOctree& lpo;
	float radius,radius2;
	unsigned int minNumNeighbors;
	IO::FilePtr outlierFile;
	size_t numOutliers;
	
	/* Constructors and destructors: */
	public:
	OutlierSaver(LidarProcessOctree& sLpo,float sRadius,unsigned int sMinNumNeighbors,IO::FilePtr sOutlierFile)
		:lpo(sLpo),
		 radius(sRadius),radius2(radius*radius),
		 minNumNeighbors(sMinNumNeighbors),
		 outlierFile(sOutlierFile),
		 numOutliers(0)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& point)
		{
		/* Count the number of points in the point's neighborhood: */
		NeighborCounter nc(point,radius2,minNumNeighbors);
		lpo.processPointsDirected(nc);
		
		if(nc.getNumNeighbors()<minNumNeighbors)
			{
			/* Point is an outlier; save it to the file: */
			outlierFile->write(point);
			++numOutliers;
			}
		}
	size_t getNumOutliers(void) const
		{
		return numOutliers;
		}
	};

int main(int argc,char* argv[])
	{
	/* Parse the command line: */
	size_t memoryCache=128;
	const char* lidarFileName=0;
	float neighborhoodRadius=1.0f;
	unsigned int minNumNeighbors=3;
	const char* outlierFileName="Outliers.bin";
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"cache")==0)
				{
				++i;
				memoryCache=size_t(atoi(argv[i]));
				}
			else if(strcasecmp(argv[i]+1,"radius")==0)
				{
				++i;
				neighborhoodRadius=float(atof(argv[i]));
				}
			else if(strcasecmp(argv[i]+1,"minNeighbors")==0)
				{
				++i;
				minNumNeighbors=(unsigned int)atoi(argv[i]);
				}
			}
		else if(lidarFileName==0)
			lidarFileName=argv[i];
		else if(outlierFileName==0)
			outlierFileName=argv[i];
		}
	if(lidarFileName==0)
		{
		std::cerr<<"No LiDAR file name provided"<<std::endl;
		return 1;
		}
	
	/* Open the LiDAR octree file: */
	LidarProcessOctree lpo(lidarFileName,memoryCache*size_t(1024*1024));
	
	/* Open the outlier file: */
	IO::SeekableFilePtr outlierFile(IO::openSeekableFile(outlierFileName,IO::File::WriteOnly));
	outlierFile->setEndianness(Misc::LittleEndian);
	
	/* Write a bogus file header: */
	outlierFile->write<unsigned int>(0);
	
	/* Process the LiDAR file: */
	OutlierSaver os(lpo,neighborhoodRadius,minNumNeighbors,outlierFile);
	lpo.processPoints(os);
	
	/* Finalize the outlier file: */
	outlierFile->setWritePosAbs(0);
	outlierFile->write<unsigned int>(os.getNumOutliers());
	std::cout<<os.getNumOutliers()<<" outlier points written to "<<outlierFileName<<std::endl;
	
	return 0;
	}
