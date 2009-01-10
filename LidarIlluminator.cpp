/***********************************************************************
LidarIlluminator - Post-processing filter to calculate normal vectors
for each point in a LiDAR data set.
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

#include <stdlib.h>
#include <string>
#include <iostream>

#include "LidarTypes.h"
#include "LidarProcessOctree.h"
#include "NormalCalculator.h"

int main(int argc,char* argv[])
	{
	const char* fileName=0;
	Scalar radius=Scalar(0);
	int cacheSize=512;
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"radius")==0)
				{
				++i;
				radius=Scalar(atof(argv[i]));
				}
			else if(strcasecmp(argv[i]+1,"cache")==0)
				{
				++i;
				cacheSize=atoi(argv[i]);
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
	
	/* Calculate normal vectors for all points in the octree: */
	std::string normalFileName=fileName;
	normalFileName.push_back('/');
	normalFileName.append("Normals");
	NodeNormalCalculator nodeNormalCalculator(lpo,radius,normalFileName.c_str());
	std::cout<<"Calculating normal vectors...   0%"<<std::flush;
	lpo.processNodesPostfix(nodeNormalCalculator);
	std::cout<<std::endl;
	
	return 0;
	}
