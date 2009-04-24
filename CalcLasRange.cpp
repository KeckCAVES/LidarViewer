/***********************************************************************
CalcLasRange - Utility program to calculate the range of point return
intensities stored in a .LAS file.
Copyright (c) 2006-2008 Oliver Kreylos

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
#include <Misc/LargeFile.h>
#include <Math/Constants.h>

void checkLasFileRange(const char* fileName,float range[2])
	{
	/* Open the LAS input file: */
	Misc::LargeFile file(fileName,"rb",Misc::LargeFile::LittleEndian);
	
	/* Read the LAS file header: */
	char signature[4];
	file.read(signature,4);
	if(memcmp(signature,"LASF",4)!=0)
		return;
	
	char dummy[32];
	file.read<unsigned short>(); // Ignore file source ID
	file.read<unsigned short>(); // Ignore reserved field
	file.read<unsigned int>(); // Ignore project ID
	file.read<unsigned short>(); // Ignore project ID
	file.read<unsigned short>(); // Ignore project ID
	file.read(dummy,8); // Ignore project ID
	file.read(dummy,2); // Ignore file version number
	file.read(dummy,32); // Ignore system identifier
	file.read(dummy,32); // Ignore generating software
	file.read<unsigned short>(); // Ignore file creation day of year
	file.read<unsigned short>(); // Ignore file creation year
	file.read<unsigned short>(); // Ignore header size
	Misc::LargeFile::Offset pointDataOffset=Misc::LargeFile::Offset(file.read<unsigned int>());
	file.read<unsigned int>(); // Ignore number of variable-length records
	unsigned char pointDataFormat=file.read<unsigned char>();
	unsigned short pointDataRecordLength=file.read<unsigned short>();
	unsigned int numPointRecords=file.read<unsigned int>();
	unsigned int numPointsByReturn[5];
	file.read(numPointsByReturn,5);
	double scale[3];
	file.read(scale,3);
	double offset[3];
	file.read(offset,3);
	double min[3],max[3];
	for(int i=0;i<3;++i)
		{
		max[i]=file.read<double>();
		min[i]=file.read<double>();
		}
	
	std::cout<<"Input file contains "<<numPointRecords<<" points."<<std::endl;
	std::cout<<"Point transformation scale and offset: ("<<scale[0]<<", "<<scale[1]<<", "<<scale[2]<<"), ("<<offset[0]<<", "<<offset[1]<<", "<<offset[2]<<")"<<std::endl;
	std::cout<<"Point set bounds: ["<<min[0]<<", "<<max[0]<<"] x ["<<min[1]<<", "<<max[1]<<"] x ["<<min[2]<<", "<<max[2]<<"]"<<std::endl;
	
	/* Read all points: */
	std::cout<<"Reading input points..."<<std::flush;
	file.seekSet(pointDataOffset);
	for(unsigned int i=0;i<numPointRecords;++i)
		{
		/* Read the point position: */
		int pos[3];
		file.read(pos,3);
		
		/* Read the point intensity: */
		float intensity=float(file.read<unsigned short>());
		if(range[0]>intensity)
			range[0]=intensity;
		if(range[1]<intensity)
			range[1]=intensity;
		
		/* Ignore the rest of the point record: */
		char dummy[4];
		file.read(dummy,4);
		file.read<unsigned short>();
		if(pointDataFormat==1)
			file.read<double>();
		}
	std::cout<<" done."<<std::endl;
	}

int main(int argc,char* argv[])
	{
	float range[2];
	range[0]=Math::Constants<float>::max;
	range[1]=Math::Constants<float>::min;
	for(int i=1;i<argc;++i)
		{
		std::cout<<"Checking file "<<argv[i]<<std::endl;
		checkLasFileRange(argv[i],range);
		}
	
	std::cout<<"Overall value range: ["<<range[0]<<", "<<range[1]<<"]"<<std::endl;
	
	return 0;
	}
