/***********************************************************************
LidarExporter - Utility to export points from LiDAR files to ASCII files
in xyzrgb format.
Copyright (c) 2010-2012 Oliver Kreylos

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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <Misc/File.h>
#include <IO/SeekableFile.h>
#include <IO/OpenFile.h>

#include "LidarTypes.h"
#include "LidarProcessOctree.h"

class PointSaver
	{
	/* Elements: */
	private:
	LidarProcessOctree::OffsetVector offset;
	Misc::File resultFile;
	size_t numPoints;
	
	/* Constructors and destructors: */
	public:
	PointSaver(const LidarProcessOctree::OffsetVector& sOffset,const char* resultFileName)
		:offset(sOffset),
		 resultFile(resultFileName,"wb"),
		 numPoints(0)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& point)
		{
		double pos[3];
		for(int i=0;i<3;++i)
			pos[i]=double(point[i])+double(offset[i]);
		fprintf(resultFile.getFilePtr(),"%.12g %.12g %.12g %u %u %u\n",pos[0],pos[1],pos[2],point.value[0],point.value[1],point.value[2]);
		++numPoints;
		}
	size_t getNumPoints(void) const
		{
		return numPoints;
		}
	};

class LasPointSaver
	{
	/* Elements: */
	private:
	IO::SeekableFilePtr lasFile;
	double scale[3],offset[3];
	double min[3],max[3];
	size_t numPoints;
	
	/* Constructors and destructors: */
	public:
	LasPointSaver(const char* lasFileName,const double sScale[3],const double sOffset[3])
		:lasFile(IO::openSeekableFile(lasFileName,IO::File::WriteOnly)),
		 numPoints(0)
		{
		for(int i=0;i<3;++i)
			{
			scale[i]=sScale[i];
			offset[i]=sOffset[i];
			min[i]=Math::Constants<double>::max;
			max[i]=Math::Constants<double>::min;
			}
		
		/* Create the initial LAS file header: */
		char signature[5]="LASF";
		lasFile->write<char>(signature,4);
		lasFile->write<unsigned short>(0);
		lasFile->write<unsigned short>(0);
		lasFile->write<unsigned int>(0);
		lasFile->write<unsigned short>(0);
		lasFile->write<unsigned short>(0);
		char dummy[32]="";
		lasFile->write<char>(dummy,8);
		lasFile->write<unsigned char>(1);
		lasFile->write<unsigned char>(2);
		lasFile->write<char>(dummy,32);
		lasFile->write<char>(dummy,32);
		lasFile->write<unsigned short>(1);
		lasFile->write<unsigned short>(2011);
		lasFile->write<unsigned short>(227);
		lasFile->write<unsigned int>(227);
		lasFile->write<unsigned int>(0);
		lasFile->write<unsigned char>(0);
		lasFile->write<unsigned short>(20);
		lasFile->write<unsigned int>(0);
		lasFile->write<unsigned int>(0);
		lasFile->write<unsigned int>(0);
		lasFile->write<unsigned int>(0);
		lasFile->write<unsigned int>(0);
		lasFile->write<unsigned int>(0);
		lasFile->write<double>(scale,3);
		lasFile->write<double>(offset,3);
		for(int i=0;i<3;++i)
			{
			lasFile->write<double>(max[i]);
			lasFile->write<double>(min[i]);
			}
		}
	~LasPointSaver(void)
		{
		/* Write the final LAS header: */
		lasFile->setWritePosAbs(107);
		lasFile->write<unsigned int>(numPoints);
		lasFile->write<unsigned int>(numPoints);
		lasFile->setWritePosAbs(179);
		for(int i=0;i<3;++i)
			{
			lasFile->write<double>(max[i]);
			lasFile->write<double>(min[i]);
			}
		}
	
	/* Methods: */
	void operator()(const LidarPoint& point)
		{
		/* Quantize the point position: */
		int p[3];
		for(int i=0;i<3;++i)
			p[i]=int(Math::floor((double(point[i])-offset[i])/scale[i]+0.5));
		
		/* Write the point record: */
		lasFile->write<int>(p,3);
		lasFile->write<unsigned short>(0);
		lasFile->write<char>(0);
		lasFile->write<char>(0);
		lasFile->write<unsigned char>(0);
		lasFile->write<unsigned char>(0);
		lasFile->write<unsigned short>(0);
		
		/* Update LAS header: */
		for(int i=0;i<3;++i)
			{
			if(min[i]>point[i])
				min[i]=point[i];
			if(max[i]<point[i])
				max[i]=point[i];
			}
		++numPoints;
		}
	size_t getNumPoints(void) const
		{
		return numPoints;
		}
	};

int main(int argc,char* argv[])
	{
	/* Process command line: */
	const char* lidarFile=0;
	bool writeLas=false;
	const char* outputFile=0;
	int cacheSize=512;
	bool haveBox=false;
	double box[6];
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"cache")==0)
				{
				++i;
				cacheSize=atoi(argv[i]);
				}
			else if(strcasecmp(argv[i]+1,"box")==0)
				{
				haveBox=true;
				for(int j=0;j<6;++j)
					{
					++i;
					box[j]=atof(argv[i]);
					}
				}
			else if(strcasecmp(argv[i]+1,"las")==0)
				writeLas=true;
			else
				std::cerr<<"Ignoring command line option "<<argv[i]<<std::endl;
			}
		else if(lidarFile==0)
			lidarFile=argv[i];
		else if(outputFile==0)
			outputFile=argv[i];
		else
			std::cerr<<"Ignoring command line argument "<<argv[i]<<std::endl;
		}
	if(lidarFile==0||outputFile==0)
		{
		std::cerr<<"Usage: "<<argv[0]<<" [-cache <cache size>] [-box <box spec>] <LiDAR file name> [-las] <output file name>"<<std::endl;
		std::cerr<<"  -cache <cache size> sets the size of the LiDAR memory cache in MB (default: 512)"<<std::endl;
		std::cerr<<"  -box <box spec> specifies a box in source coordinates from which to export points (default: export all points)"<<std::endl;
		std::cerr<<"     box specification: <min_x> <min_y> <min_z> <max_x> <max_y> <max_z>"<<std::endl;
		std::cerr<<"  -las requests to write exported points into a LAS-like file (default: write into ASCII file)"<<std::endl;
		return 1;
		}
	
	/* Open the input and output files: */
	LidarProcessOctree lpo(lidarFile,size_t(cacheSize)*1024*1024);
	Box lbox;
	if(haveBox)
		{
		/* Transform the box to LiDAR coordinates: */
		for(int i=0;i<3;++i)
			{
			lbox.min[i]=Box::Scalar(box[i]-lpo.getOffset()[i]);
			lbox.max[i]=Box::Scalar(box[3+i]-lpo.getOffset()[i]);
			}
		}
	
	if(writeLas)
		{
		double scale[3]={0.001,0.001,0.001};
		double offset[3];
		for(int i=0;i<3;++i)
			offset[i]=lpo.getDomain().getCenter(i);
		LasPointSaver ps(outputFile,scale,offset);
		
		/* Process the LiDAR file: */
		if(haveBox)
			{
			/* Exports points from inside the box: */
			lpo.processPointsInBox(lbox,ps);
			}
		else
			{
			/* Export all points: */
			lpo.processPoints(ps);
			}
		
		/* Print statistics: */
		std::cout<<ps.getNumPoints()<<" points saved"<<std::endl;
		}
	else
		{
		PointSaver ps(lpo.getOffset(),outputFile);
		
		/* Process the LiDAR file: */
		if(haveBox)
			{
			/* Exports points from inside the box: */
			lpo.processPointsInBox(lbox,ps);
			}
		else
			{
			/* Export all points: */
			lpo.processPoints(ps);
			}
		
		/* Print statistics: */
		std::cout<<ps.getNumPoints()<<" points saved"<<std::endl;
		}
	
	return 0;
	}
