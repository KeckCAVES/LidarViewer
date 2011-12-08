/***********************************************************************
LidarIlluminator - Post-processing filter to calculate normal vectors
for each point in a LiDAR data set.
Copyright (c) 2008-2011 Oliver Kreylos

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
#include <string>
#include <iostream>

#include "LidarTypes.h"
#include "LidarProcessOctree.h"
#include "NormalCalculator.h"

int main(int argc,char* argv[])
	{
	const char* fileName=0;
	unsigned int maxNumNeighbors=0;
	Scalar radius=Scalar(0);
	int cacheSize=512;
	const char* outlierFileName=0;
	unsigned int numThreads=1;
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"neighbors")==0)
				{
				++i;
				maxNumNeighbors=atoi(argv[i]);
				}
			else if(strcasecmp(argv[i]+1,"radius")==0)
				{
				++i;
				radius=Scalar(atof(argv[i]));
				}
			else if(strcasecmp(argv[i]+1,"cache")==0)
				{
				++i;
				cacheSize=atoi(argv[i]);
				}
			else if(strcasecmp(argv[i]+1,"threads")==0)
				{
				++i;
				numThreads=atoi(argv[i]);
				}
			else if(strcasecmp(argv[i]+1,"saveOutliers")==0)
				{
				++i;
				outlierFileName=argv[i];
				}
			}
		else if(fileName==0)
			fileName=argv[i];
		}
	if(fileName==0)
		{
		std::cerr<<"No file name provided"<<std::endl;
		return 1;
		}
	
	/* Create a processing octree: */
	LidarProcessOctree lpo(fileName,size_t(cacheSize)*size_t(1024*1024));
	
	/* Create the output file name: */
	std::string normalFileName=fileName;
	normalFileName.push_back('/');
	normalFileName.append("Normals");
	
	/* Check which version of the normal calculator to use: */
	if(maxNumNeighbors>0)
		{
		/* Create a fixed-neighborhood normal vector calculation functor: */
		if(radius>Scalar(0))
			{
			NumberRadiusNormalCalculator nc(maxNumNeighbors,radius);
			NodeNormalCalculator<NumberRadiusNormalCalculator> nodeNormalCalculator(lpo,nc,normalFileName.c_str(),numThreads);
			
			/* Enable saving of outliers if requested: */
			if(outlierFileName!=0)
				nodeNormalCalculator.saveOutlierPoints(outlierFileName);
			
			/* Calculate normal vectors for all points in the octree: */
			std::cout<<"Calculating normal vectors...   0%"<<std::flush;
			lpo.processNodesPostfix(nodeNormalCalculator);
			std::cout<<std::endl;
			}
		else
			{
			NumberRadiusNormalCalculator nc(maxNumNeighbors);
			NodeNormalCalculator<NumberRadiusNormalCalculator> nodeNormalCalculator(lpo,nc,normalFileName.c_str(),numThreads);
			
			/* Enable saving of outliers if requested: */
			if(outlierFileName!=0)
				nodeNormalCalculator.saveOutlierPoints(outlierFileName);
			
			/* Calculate normal vectors for all points in the octree: */
			std::cout<<"Calculating normal vectors...   0%"<<std::flush;
			lpo.processNodesPostfix(nodeNormalCalculator);
			std::cout<<std::endl;
			}
		}
	else
		{
		/* Create a fixed-radius normal vector calculation functor: */
		RadiusNormalCalculator nc(radius);
		NodeNormalCalculator<RadiusNormalCalculator> nodeNormalCalculator(lpo,nc,normalFileName.c_str(),numThreads);
		
		/* Enable saving of outliers if requested: */
		if(outlierFileName!=0)
			nodeNormalCalculator.saveOutlierPoints(outlierFileName);
		
		/* Calculate normal vectors for all points in the octree: */
		std::cout<<"Calculating normal vectors...   0%"<<std::flush;
		lpo.processNodesPostfix(nodeNormalCalculator);
		std::cout<<std::endl;
		}
	
	return 0;
	}
