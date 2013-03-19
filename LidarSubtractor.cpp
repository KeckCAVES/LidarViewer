/***********************************************************************
LidarSubtractor - Post-processing filter to subtract a (relatively
small) point set from a LiDAR data set and create a new LiDAR data set.
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

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <Misc/StandardValueCoders.h>
#include <Misc/ConfigurationFile.h>
#include <IO/File.h>
#include <IO/OpenFile.h>
#include <Geometry/ArrayKdTree.h>

#include "LidarTypes.h"
#include "LidarProcessOctree.h"
#include "PointAccumulator.h"
#include "LidarOctreeCreator.h"
#include "SubtractorHelper.h"

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
					pa.addPoint(PointAccumulator::Point(node[i].getComponents()),PointAccumulator::Color(node[i].value.getRgba()));
					}
				}
			}
		}
	};

int main(int argc,char* argv[])
	{
	/* Set default values for all parameters: */
	unsigned int memoryCacheSize=512;
	unsigned int tempOctreeMaxNumPointsPerNode=4096;
	std::string tempOctreeFileNameTemplate="/tmp/LidarPreprocessorTempOctree";
	unsigned int maxNumPointsPerNode=4096;
	int numThreads=1;
	std::string tempPointFileNameTemplate="/tmp/LidarPreprocessorTempPoints";
	
	try
		{
		/* Open LidarViewer's configuration file: */
		Misc::ConfigurationFile configFile(LIDARVIEWER_CONFIGFILENAME);
		Misc::ConfigurationFileSection cfg=configFile.getSection("/LidarPreprocessor");
		
		/* Override program settings from configuration file: */
		memoryCacheSize=cfg.retrieveValue<unsigned int>("./memoryCacheSize",memoryCacheSize);
		tempOctreeMaxNumPointsPerNode=cfg.retrieveValue<unsigned int>("./tempOctreeMaxNumPointsPerNode",tempOctreeMaxNumPointsPerNode);
		tempOctreeFileNameTemplate=cfg.retrieveValue<std::string>("./tempOctreeFileNameTemplate",tempOctreeFileNameTemplate);
		maxNumPointsPerNode=cfg.retrieveValue<unsigned int>("./maxNumPointsPerNode",maxNumPointsPerNode);
		numThreads=cfg.retrieveValue<int>("./numThreads",numThreads);
		tempPointFileNameTemplate=cfg.retrieveValue<std::string>("./tempPointFileNameTemplate",tempPointFileNameTemplate);
		}
	catch(std::runtime_error err)
		{
		/* Just ignore the error */
		}
	
	/* Parse the command line and load the input files: */
	const char* baseFileName=0; // Name of the base LiDAR file
	unsigned baseMemoryCacheSize=64; // Memory cache size for base octree in MB
	const char* subtractFileName=0; // Name of ASCII or binary file containing the point set to subtract
	int asciiColumnIndices[3]={0,1,2}; // Column indices of x, y, z point components in ASCII file
	Scalar epsilon=Scalar(1.0e-7); // Maximum match point distance
	PointAccumulator::Vector pointOffset=PointAccumulator::Vector::zero; // Offset vector added to points during octree creation
	const char* outputFileName=0; // Name of resulting LiDAR octree file
	bool haveOffset=false; // Flag whether an explicit point offset was specified on the command line
	bool offsetSubtractPoints=true; // Flag whether to offset subtraction points to native octree coordinates
	
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
					maxNumPointsPerNode=(unsigned int)(atoi(argv[i]));
				else
					std::cerr<<"Dangling -np flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"nt")==0)
				{
				++i;
				if(i<argc)
					numThreads=atoi(argv[i]);
				else
					std::cerr<<"Dangling -nt flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"ooc")==0)
				{
				++i;
				if(i<argc)
					memoryCacheSize=(unsigned int)(atoi(argv[i]));
				else
					std::cerr<<"Dangling -ooc flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"to")==0)
				{
				++i;
				if(i<argc)
					tempOctreeFileNameTemplate=argv[i];
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
			else if(strcasecmp(argv[i]+1,"lasOffset")==0)
				{
				if(i+3<argc)
					{
					for(int j=0;j<3;++j)
						{
						++i;
						pointOffset[j]=atof(argv[i]);
						}
					haveOffset=true;
					}
				else
					std::cerr<<"Dangling -lasOffset flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"lasOffsetFile")==0)
				{
				++i;
				if(i<argc)
					{
					/* Read the point offset from a binary file: */
					try
						{
						IO::FilePtr offsetFile=IO::openFile(argv[i]);
						offsetFile->setEndianness(Misc::LittleEndian);
						offsetFile->read(pointOffset.getComponents(),3);
						haveOffset=true;
						}
					catch(std::runtime_error err)
						{
						/* Print a warning and carry on: */
						std::cerr<<"Ignoring lasOffsetFile argument due to error "<<err.what()<<" when reading file "<<argv[i]<<std::endl;
						}
					}
				else
					std::cerr<<"Dangling -lasOffsetFile flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"noOffset")==0)
				offsetSubtractPoints=false;
			else if(strcasecmp(argv[i]+1,"eps")==0)
				{
				++i;
				if(i<argc)
					epsilon=Scalar(atof(argv[i]));
				else
					std::cerr<<"Dangling -eps flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"booc")==0)
				{
				++i;
				if(i<argc)
					baseMemoryCacheSize=(unsigned int)(atoi(argv[i]));
				else
					std::cerr<<"Dangling -booc flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"columns")==0)
				{
				if(i+3<argc)
					{
					for(int j=0;j<3;++j)
						{
						++i;
						asciiColumnIndices[j]=atoi(argv[i]);
						}
					}
				else
					std::cerr<<"Dangling -columns flag on command line"<<std::endl;
				}
			}
		else if(baseFileName==0)
			{
			/* Store the LiDAR file name: */
			baseFileName=argv[i];
			}
		else if(subtractFileName==0)
			{
			/* Store the subtraction file name: */
			subtractFileName=argv[i];
			}
		else
			std::cerr<<"Ignoring command line argument "<<argv[i]<<std::endl;
		}
	
	/* Open the base LiDAR file: */
	LidarProcessOctree* basePoints;
	try
		{
		/* Create a processing octree: */
		basePoints=new LidarProcessOctree(baseFileName,size_t(baseMemoryCacheSize)*size_t(1024*1024));
		}
	catch(std::runtime_error err)
		{
		std::cerr<<"Cannot open LiDAR file "<<baseFileName<<" due to exception "<<err.what()<<"; terminating"<<std::endl;
		return 1;
		}
	
	/* Create a point accumulator: */
	PointAccumulator pa;
	pa.setMemorySize(memoryCacheSize,tempOctreeMaxNumPointsPerNode);
	pa.setTempOctreeFileNameTemplate(tempOctreeFileNameTemplate+"XXXXXX");
	if(!haveOffset)
		{
		/* Use the base point set's point offset for the resulting file: */
		pointOffset=-basePoints->getOffset();
		}
	else
		{
		PointAccumulator::Vector totalOffset=pointOffset-PointAccumulator::Vector(basePoints->getOffset());
		if(totalOffset!=PointAccumulator::Vector::zero)
			pa.setPointOffset(totalOffset);
		}
	
	/* Load the subtraction point set into a kd-tree: */
	PointKdTree* subtractPointTree=loadSubtractSet(subtractFileName,offsetSubtractPoints?Geometry::Vector<double,3>(basePoints->getOffset()):Geometry::Vector<double,3>::zero);
	if(subtractPointTree==0)
		return 1;
	
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
	LidarOctreeCreator tree(pa.getMaxNumCacheablePoints(),maxNumPointsPerNode,numThreads,pa.getTempOctrees(),tempPointFileNameTemplate+"XXXXXX");
	
	/* Delete the temporary point octrees: */
	pa.deleteTempOctrees();
	
	/* Write the octree structure and data to the destination LiDAR file: */
	tree.write(size_t(memoryCacheSize)*size_t(1024*1024),outputFileName);
	
	/* Check if a point offset was defined: */
	if(pointOffset!=PointAccumulator::Vector::zero)
		{
		/* Write the point offsets to an offset file: */
		std::string offsetFileName=outputFileName;
		offsetFileName.append("/Offset");
		IO::FilePtr offsetFile(IO::openFile(offsetFileName.c_str(),IO::File::WriteOnly));
		offsetFile->setEndianness(Misc::LittleEndian);
		offsetFile->write(pointOffset.getComponents(),3);
		}
	
	return 0;
	}
