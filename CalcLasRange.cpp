/***********************************************************************
CalcLasRange - Utility program to calculate the spatial extent and range
of point return intensities or colors stored in a set of .LAS files.
Copyright (c) 2006-2013 Oliver Kreylos

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
#include <iostream>
#include <fstream>
#include <IO/File.h>
#include <Comm/OpenFile.h>
#include <Math/Math.h>
#include <Math/Constants.h>

size_t checkLasFileRange(const char* fileName,double bbox[6],bool getColorRange,float intensityRange[2],float rgbRange[3][2])
	{
	/* Open the LAS input file: */
	IO::FilePtr file(Comm::openFile(fileName));
	file->resizeReadBuffer(1048576);
	file->setEndianness(Misc::LittleEndian);
	
	/* Process the LAS file header: */
	unsigned char pointDataFormat=0;
	unsigned int numPointRecords=0;
	unsigned int pointDataOffset=0;
	try
		{
		/* Read the LAS file header: */
		char signature[4];
		file->read(signature,4);
		if(memcmp(signature,"LASF",4)!=0)
			{
			std::cout<<"File "<<fileName<<" is not a LAS file"<<std::endl;
			return 0;
			}
		
		file->read<unsigned short>(); // Ignore file source ID
		file->read<unsigned short>(); // Ignore reserved field
		file->read<unsigned int>(); // Ignore project ID
		file->read<unsigned short>(); // Ignore project ID
		file->read<unsigned short>(); // Ignore project ID
		file->skip<char>(8); // Ignore project ID
		file->skip<char>(2); // Ignore file version number
		file->skip<char>(32); // Ignore system identifier
		file->skip<char>(32); // Ignore generating software
		file->read<unsigned short>(); // Ignore file creation day of year
		file->read<unsigned short>(); // Ignore file creation year
		file->read<unsigned short>(); // Ignore header size
		pointDataOffset=file->read<unsigned int>();
		file->read<unsigned int>(); // Ignore number of variable-length records
		pointDataFormat=file->read<unsigned char>();
		unsigned short pointDataRecordLength=file->read<unsigned short>();
		numPointRecords=file->read<unsigned int>();
		unsigned int numPointsByReturn[5];
		file->read(numPointsByReturn,5);
		double scale[3];
		file->read(scale,3);
		double offset[3];
		file->read(offset,3);
		double min[3],max[3];
		for(int i=0;i<3;++i)
			{
			max[i]=file->read<double>();
			min[i]=file->read<double>();
			if(bbox[i]>min[i])
				bbox[i]=min[i];
			if(bbox[3+i]<max[i])
				bbox[3+i]=max[i];
			}
		
		#if 0
		std::cout<<"Input file contains "<<numPointRecords<<" points."<<std::endl;
		std::cout<<"Point record size: "<<pointDataRecordLength<<" bytes"<<std::endl;
		std::cout<<"Point transformation scale and offset: ("<<scale[0]<<", "<<scale[1]<<", "<<scale[2]<<"), ("<<offset[0]<<", "<<offset[1]<<", "<<offset[2]<<")"<<std::endl;
		std::cout<<"Point set bounds: ["<<min[0]<<", "<<max[0]<<"] x ["<<min[1]<<", "<<max[1]<<"] x ["<<min[2]<<", "<<max[2]<<"]"<<std::endl;
		#endif
		}
	catch(std::runtime_error err)
		{
		std::cout<<"File "<<fileName<<" caused exception "<<err.what()<<" while reading LAS file header"<<std::endl;
		return 0;
		}
	
	if(!getColorRange)
		{
		std::cout<<" done."<<std::endl;
		return numPointRecords;
		}
	
	/* Read all points: */
	std::cout<<"Reading input points..."<<std::flush;
	unsigned int pointIndex=0;
	try
		{
		/* Skip to the beginning of the point data records: */
		if(pointDataOffset<227)
			{
			std::cout<<"File "<<fileName<<" has invalid LAS file header"<<std::endl;
			return 0;
			}
		file->skip<unsigned char>(pointDataOffset-227);
		
		for(pointIndex=0;pointIndex<numPointRecords;++pointIndex)
			{
			/* Read the point position: */
			int pos[3];
			file->read(pos,3);
			
			/* Read the point intensity: */
			float intensity=float(file->read<unsigned short>());
			if(intensityRange[0]>intensity)
				intensityRange[0]=intensity;
			if(intensityRange[1]<intensity)
				intensityRange[1]=intensity;
			
			/* Skip irrelevant information: */
			file->skip<char>(4);
			file->read<unsigned short>();
			if(pointDataFormat&0x1)
				file->read<double>();
			if(pointDataFormat>=2)
				{
				/* Read the point color: */
				unsigned short rgb[3];
				file->read(rgb,3);
				for(int j=0;j<3;++j)
					{
					float value=float(rgb[j]);
					if(rgbRange[j][0]>value)
						rgbRange[j][0]=value;
					if(rgbRange[j][1]<value)
						rgbRange[j][1]=value;
					}
				}
			}
		std::cout<<" done."<<std::endl;
		}
	catch(std::runtime_error err)
		{
		std::cout<<" terminated after "<<pointIndex<<" point records due to exception "<<err.what()<<std::endl;
		}
	
	return pointIndex;
	}

int main(int argc,char* argv[])
	{
	double bbox[6];
	for(int i=0;i<3;++i)
		{
		bbox[i]=Math::Constants<double>::max;
		bbox[3+i]=Math::Constants<double>::min;
		}
	float intensityRange[2];
	intensityRange[0]=Math::Constants<float>::max;
	intensityRange[1]=Math::Constants<float>::min;
	float rgbRange[3][2];
	for(int i=0;i<3;++i)
		{
		rgbRange[i][0]=Math::Constants<float>::max;
		rgbRange[i][1]=Math::Constants<float>::min;
		}
	
	/* Parse the command line: */
	bool getColorRange=true;
	const char* sceneGraphName=0;
	int argi;
	for(argi=1;argi<argc&&argv[argi][0]=='-';++argi)
		{
		if(strcasecmp(argv[argi]+1,"boxOnly")==0)
			{
			/* Don't read all points and calculate intensity and color ranges: */
			getColorRange=false;
			}
		else if(strcasecmp(argv[argi]+1,"showCoverage")==0)
			{
			/* Create a scene graph file showing the extents and names of all LAS files: */
			++argi;
			sceneGraphName=argv[argi];
			}
		}
	
	std::ofstream sceneGraph;
	if(sceneGraphName!=0)
		{
		/* Create a scene graph file: */
		sceneGraph.open(sceneGraphName);
		sceneGraph<<"#VRML V2.0 utf8"<<std::endl;
		sceneGraph.precision(12);
		}
	
	/* Process all LAS files from the command line: */
	int firstLasFile=argi;
	int numLasFiles=argc-firstLasFile;
	size_t totalNumPoints=0;
	for(;argi<argc;++argi)
		{
		std::cout<<"Checking file "<<argi+1-firstLasFile<<" of "<<numLasFiles<<":"<<argv[argi]<<"..."<<std::flush;
		if(sceneGraphName!=0)
			{
			/* Create a new bounding box for this LAS file: */
			double fileBbox[6];
			for(int i=0;i<3;++i)
				{
				fileBbox[i]=Math::Constants<double>::max;
				fileBbox[3+i]=Math::Constants<double>::min;
				}
			
			/* Get the LAS file's extents: */
			totalNumPoints+=checkLasFileRange(argv[argi],fileBbox,getColorRange,intensityRange,rgbRange);
			
			/* Create a scene graph node for this LAS file: */
			sceneGraph<<"\n\
				Shape\n\
					{\n";
			
			if(argi==firstLasFile)
				{
				sceneGraph<<"\
					appearance DEF Red Appearance\n\
						{\n\
						material Material\n\
							{\n\
							diffuseColor 0.0 0.0 0.0\n\
							specularColor 0.0 0.0 0.0\n\
							emissiveColor 1.0 0.0 0.0\n\
							}\n\
						}\n";
				}
			else
				{
				sceneGraph<<"\
					appearance USE Red\n";
				}
			
			sceneGraph<<"\
					geometry IndexedLineSet\n\
						{\n\
						coord Coordinate\n\
							{\n\
							point\n\
								[\n";
			
			for(int i=0;i<8;++i)
				sceneGraph<<"\t\t\t\t\t\t\t\t"<<fileBbox[(i&0x1)?3:0]<<' '<<fileBbox[(i&0x2)?4:1]<<' '<<fileBbox[(i&0x4)?5:2]<<std::endl;
			sceneGraph<<"\
								]\n\
							}\n\
						coordIndex\n\
							[\n\
							0,1,3,2,0,4,5,7,6,4,-1\n\
							1,5,-1\n\
							3,7,-1\n\
							2,6,-1\n\
							]\n\
						lineWidth 3.0\n\
						}\n\
					}\n";
			
			sceneGraph<<"\n\
				Transform\n\
					{\n";
			sceneGraph<<"\t\t\t\t\ttranslation "<<Math::mid(fileBbox[0],fileBbox[3])<<" "<<Math::mid(fileBbox[1],fileBbox[4])<<" "<<fileBbox[5]<<std::endl;
			sceneGraph<<"\
					children\n\
						[\n\
						Shape\n\
							{\n\
							appearance USE Red\n\
							geometry Text\n\
								{\n";
			
			/* Show the LAS file's base name: */
			const char* basePtr=argv[argi];
			for(const char* fnPtr=argv[argi];*fnPtr!='\0';++fnPtr)
				if(*fnPtr=='/')
					basePtr=fnPtr+1;
			sceneGraph<<"\t\t\t\t\t\t\t\tstring \""<<basePtr<<"\""<<std::endl;
			
			if(argi==firstLasFile)
				{
				sceneGraph<<"\
								fontStyle DEF Sans FontStyle\n\
									{\n\
									family \"SANS\"\n\
									style \"PLAIN\"\n\
									size 75.0\n\
									justify [\"MIDDLE\" \"MIDDLE\"]\n\
									}\n";
				}
			else
				{
				sceneGraph<<"\
								fontStyle USE Sans\n";
				}
			
			sceneGraph<<"\t\t\t\t\t\t\t\tmaxExtent "<<(fileBbox[3]-fileBbox[0])*0.9<<std::endl;
			
			sceneGraph<<"\
								}\n\
							}\n\
						]\n\
					}\n";
			
			/* Merge the LAS file's bounding box with the overall box: */
			for(int i=0;i<3;++i)
				{
				if(bbox[i]>fileBbox[i])
					bbox[i]=fileBbox[i];
				if(bbox[3+i]<fileBbox[3+i])
					bbox[3+i]=fileBbox[3+i];
				}
			}
		else
			totalNumPoints+=checkLasFileRange(argv[argi],bbox,getColorRange,intensityRange,rgbRange);
		}
	std::cout<<"Total number of points: "<<totalNumPoints<<std::endl;
	std::cout<<"Overall point bounding box: ["<<bbox[0]<<", "<<bbox[3]<<"] x ["<<bbox[1]<<", "<<bbox[4]<<"] x ["<<bbox[2]<<", "<<bbox[5]<<"]"<<std::endl;
	if(getColorRange)
		{
		std::cout<<"Overall intensity range: ["<<intensityRange[0]<<", "<<intensityRange[1]<<"]"<<std::endl;
		std::cout<<"Overall RGB range: ["<<rgbRange[0][0]<<", "<<rgbRange[0][1]<<"], ["<<rgbRange[1][0]<<", "<<rgbRange[1][1]<<"], ["<<rgbRange[2][0]<<", "<<rgbRange[2][1]<<"]"<<std::endl;
		}
	
	if(sceneGraphName!=0)
		{
		/* Close the scene graph file: */
		sceneGraph.close();
		}
	
	std::cout<<"Recommended offset vector for LiDAR pre-processor: -lasOffset";
	for(int i=0;i<3;++i)
		std::cout<<' '<<-Math::mid(bbox[i],bbox[3+i]);
	std::cout<<std::endl;
	
	return 0;
	}
