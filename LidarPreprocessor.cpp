/***********************************************************************
New version of LiDAR data preprocessor.
Copyright (c) 2005-2008 Oliver Kreylos

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
#include <math.h>
#include <iostream>
#include <vector>
#include <Misc/ThrowStdErr.h>
#include <Misc/Endianness.h>
#include <Misc/File.h>
#include <Misc/LargeFile.h>
#include <Misc/Timer.h>
#include <Math/Math.h>
#include <Math/Constants.h>

#include "TempPointOctree.h"
#include "PointOctreeCreator.h"
#include "LidarOctreeCreator.h"

namespace {

/**************
Helper classes:
**************/

class OctreePointLoader // Class to load point sets from input file and manage out-of-core storage
	{
	/* Elements: */
	private:
	unsigned int maxNumCachablePoints; // Maximum number of points allowed in memory at a time
	std::vector<OctreePoint>* points; // Pointer to vector holding current in-memory point set
	unsigned int maxNumPointsPerNode; // Maximum number of points per node in the temporary octrees
	std::vector<TempPointOctree*> octrees; // List of octrees holding out-of-memory point sets
	std::string tempPointOctreeFileNameTemplate; // File name template for temporary octrees
	
	/* Private methods: */
	void savePoints(void) // Saves the current in-memory point set to a temporary octree file
		{
		/* Create a temporary octree for the current in-memory point set: */
		std::cout<<std::endl<<"Storing "<<points->size()<<" points as temporary octree..."<<std::flush;
		char tpofnt[1024];
		strcpy(tpofnt,tempPointOctreeFileNameTemplate.c_str());
		TempPointOctree* tpo=new TempPointOctree(tpofnt,maxNumPointsPerNode,&(*points)[0],points->size());
		std::cout<<" done"<<std::endl;
		octrees.push_back(tpo);
		
		/* Clear the point set: */
		points->clear();
		};
	
	/* Constructors and destructors: */
	public:
	OctreePointLoader(void) // Creates an in-core point loader
		:maxNumCachablePoints(~0x0U),
		 points(new std::vector<OctreePoint>),
		 tempPointOctreeFileNameTemplate("/tmp/LidarPreprocessorTempOctreeXXXXXX")
		{
		};
	~OctreePointLoader(void) // Destroys all point sets
		{
		/* Delete all point octrees: */
		for(std::vector<TempPointOctree*>::iterator poIt=octrees.begin();poIt!=octrees.end();++poIt)
			delete *poIt;
		
		/* Delete the in-memory point set: */
		delete points;
		};
	
	/* Methods: */
	void setTempFileNameTemplate(const char* newTempFileNameTemplate) // Sets the template for temporary octree file names
		{
		tempPointOctreeFileNameTemplate=newTempFileNameTemplate;
		};
	void setMemorySize(unsigned int memorySize,unsigned int newMaxNumPointsPerNode) // Limits the point loader to the given amount of memory in megabytes
		{
		/* Set the memory limit: */
		maxNumCachablePoints=(memorySize*1024U*1024U+sizeof(OctreePoint)-1)/sizeof(OctreePoint);
		maxNumPointsPerNode=newMaxNumPointsPerNode;
		
		/* Check if the current point set is already too large: */
		if(points->size()>maxNumCachablePoints)
			{
			/* Save the current point set: */
			savePoints();
			}
		
		/* Allocate the maximum amount of memory for the in-memory point set: */
		points->reserve(maxNumCachablePoints);
		};
	void addOctree(const char* octreeFileNameStem) // Adds an existing LiDAR octree to the current point set
		{
		/* Load the existing octree as a temporary octree: */
		TempPointOctree* tpo=new TempPointOctree(octreeFileNameStem,true);
		octrees.push_back(tpo);
		};
	void addPoint(const OctreePoint& op) // Pushes an octree point into the current point set
		{
		/* Check if the current in-memory point set is too big: */
		if(points->size()==maxNumCachablePoints)
			{
			/* Save the current point set: */
			savePoints();
			}
		
		/* Store the new point: */
		points->push_back(op);
		};
	bool finishReading(void) // Finishes reading points from source files; returns true if out-of-core processing is necessary
		{
		#if 1
		if(!points->empty())
			{
			/* Write the leftover in-memory points into another temporary octree: */
			savePoints();
			}
		
		/* Delete the in-memory point set: */
		delete points;
		points=0;
		
		/* Signal out-of-core processing: */
		return true;
		#else
		if(!octrees.empty())
			{
			if(!points->empty())
				{
				/* Write the leftover in-memory points into another temporary octree: */
				savePoints();
				}
			
			/* Delete the in-memory point set: */
			delete points;
			points=0;
			
			/* Signal out-of-core processing: */
			return true;
			}
		else
			return false;
		#endif
		};
	unsigned int getNumPoints(void) const // Returns the number of points in the in-memory point set
		{
		return points->size();
		};
	OctreePoint* getPoints(void) const // Returns the in-memory point set
		{
		return &(*points)[0];
		};
	std::vector<TempPointOctree*>& getTempOctrees(void) // Returns the list of temporary octrees
		{
		return octrees;
		};
	unsigned int getMaxNumCachablePoints(void) const // Returns the maximum number of points to be held in memory
		{
		return maxNumCachablePoints;
		};
	void deleteTempOctrees(void) // Deletes the temporary point octrees
		{
		/* Delete all point octrees: */
		for(std::vector<TempPointOctree*>::iterator poIt=octrees.begin();poIt!=octrees.end();++poIt)
			delete *poIt;
		octrees.clear();
		};
	};

}

void loadPointFileBin(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the binary input file: */
	Misc::LargeFile file(fileName,"rb",Misc::LargeFile::LittleEndian);
	
	/* Read the number of points in the file: */
	unsigned int numPoints=file.read<unsigned int>();
	
	/* Read all points: */
	for(unsigned int i=0;i<numPoints;++i)
		{
		/* Read the point position and intensity from the input file: */
		float rp[4];
		file.read(rp,4);
		
		/* Store the point position: */
		OctreePoint p;
		for(int j=0;j<3;++j)
			p[j]=rp[j];
		
		/* Convert the point intensity to and RGB color according to the color mask: */
		for(int j=0;j<3;++j)
			p.value[j]=GLubyte(Math::floor(rp[3]*colorMask[j]+0.5f));
		p.value[3]=GLubyte(255);
		
		/* Store the point: */
		opl.addPoint(p);
		}
	}

void loadPointFileBinRgb(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the binary input file: */
	Misc::LargeFile file(fileName,"rb",Misc::LargeFile::LittleEndian);
	
	/* Read the number of points in the file: */
	unsigned int numPoints=file.read<unsigned int>();
	
	/* Read all points: */
	for(unsigned int i=0;i<numPoints;++i)
		{
		/* Read the point position and color from the input file: */
		float rp[3];
		file.read(rp,3);
		GLubyte rcol[4];
		file.read(rcol,4);
		
		/* Store the point position: */
		OctreePoint p;
		for(int j=0;j<3;++j)
			p[j]=rp[j];
		
		/* Modify the RGB color according to the color mask: */
		for(int j=0;j<3;++j)
			p.value[j]=GLubyte(Math::floor(float(rcol[j])*colorMask[j]+0.5f));
		p.value[3]=GLubyte(255);
		
		/* Store the point: */
		opl.addPoint(p);
		}
	}

void loadPointFileTxt(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the ASCII input file: */
	Misc::File file(fileName,"rt");
	
	/* Read all lines from the input file: */
	while(!file.eof())
		{
		/* Read the next line from the file: */
		char line[256];
		file.gets(line,sizeof(line));
		
		/* Parse the point coordinates from the line: */
		OctreePoint p;
		if(sscanf(line,"%f %f %f",&p[0],&p[1],&p[2])==3)
			{
			for(int i=0;i<4;++i)
				p.value[i]=GLubyte(255);
			
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	}

void loadPointFileXyzi(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the ASCII input file: */
	Misc::File file(fileName,"rt");
	
	/* Read all lines from the input file: */
	while(!file.eof())
		{
		/* Read the next line from the file: */
		char line[256];
		file.gets(line,sizeof(line));
		
		/* Parse the point coordinates and intensity from the line: */
		OctreePoint p;
		float intensity;
		if(sscanf(line,"%f %f %f %f",&p[0],&p[1],&p[2],&intensity)==4&&intensity!=0.0f)
			{
			/* Convert the read intensity to an RGB color according to the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=GLubyte(Math::floor(intensity*colorMask[i]+0.5f));
			p.value[3]=GLubyte(255);
			
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	}

void loadPointFileGenericCsv(OctreePointLoader& opl,const char* fileName,const int columnIndices[4],const float colorMask[3])
	{
	/* Create the mapping from column indices to point components: */
	int maxColumnIndex=columnIndices[0];
	for(int i=1;i<4;++i)
		if(maxColumnIndex<columnIndices[i])
			maxColumnIndex=columnIndices[i];
	int* componentColumnIndices=new int[maxColumnIndex+1];
	for(int i=0;i<=maxColumnIndex;++i)
		componentColumnIndices[i]=-1;
	int numComponents=0;
	for(int i=0;i<4;++i)
		if(columnIndices[i]>=0)
			{
			componentColumnIndices[columnIndices[i]]=i;
			++numComponents;
			}
	
	/* Open the ASCII input file: */
	Misc::File file(fileName,"rt");
	
	/* Read all lines from the input file: */
	float intensityMin=Math::Constants<float>::max;
	float intensityMax=Math::Constants<float>::min;
	char line[256];
	while(!file.eof())
		{
		/* Read the next line from the file: */
		file.gets(line,sizeof(line));
		
		/* Separate and parse the columns: */
		int numParsedColumns=0;
		float componentValues[4];
		componentValues[3]=1.0f;
		char* columnStart=line;
		for(int columnIndex=0;columnIndex<=maxColumnIndex;++columnIndex)
			{
			/* Find the end of the current column: */
			char* columnEnd;
			for(columnEnd=columnStart;*columnEnd!=','&&*columnEnd!='\n';++columnEnd)
				;
			char terminator=*columnEnd;
			
			if(componentColumnIndices[columnIndex]>=0)
				{
				/* Skip whitespace and quotes at the beginning and the end: */
				char* valueStart=columnStart;
				while(valueStart!=columnEnd&&(isspace(*valueStart)||*valueStart=='\"'))
					++valueStart;
				char* valueEnd=columnEnd;
				while(valueEnd!=valueStart&&(isspace(valueEnd[-1])||valueEnd[-1]=='\"'))
					--valueEnd;
				
				if(valueStart!=valueEnd)
					{
					/* Parse the column's value: */
					*valueEnd='\0';
					char* conversionEnd;
					componentValues[componentColumnIndices[columnIndex]]=strtof(valueStart,&conversionEnd);
					if(conversionEnd!=valueStart)
						++numParsedColumns;
					}
				}
			
			/* Go to the next column: */
			if(terminator=='\n')
				break;
			columnStart=columnEnd+1;
			}
		
		if(numParsedColumns==numComponents)
			{
			/* Store the point position: */
			OctreePoint p;
			for(int i=0;i<3;++i)
				p[i]=componentValues[i];
			
			/* Convert the read intensity to an RGB color according to the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=GLubyte(Math::floor(componentValues[3]*colorMask[i]+0.5f));
			p.value[3]=GLubyte(255);
			
			if(intensityMin>componentValues[3])
				intensityMin=componentValues[3];
			if(intensityMax<componentValues[3])
				intensityMax=componentValues[3];
			
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	std::cout<<"Point data intensity range: "<<intensityMin<<", "<<intensityMax<<std::endl;
	
	/* Clean up: */
	delete[] componentColumnIndices;
	}

void loadPointFileGenericCsvRgb(OctreePointLoader& opl,const char* fileName,const int columnIndices[6],const float colorMask[3])
	{
	/* Create the mapping from column indices to point components: */
	int maxColumnIndex=columnIndices[0];
	for(int i=1;i<6;++i)
		if(maxColumnIndex<columnIndices[i])
			maxColumnIndex=columnIndices[i];
	int* componentColumnIndices=new int[maxColumnIndex+1];
	for(int i=0;i<=maxColumnIndex;++i)
		componentColumnIndices[i]=-1;
	int numComponents=0;
	for(int i=0;i<6;++i)
		if(columnIndices[i]>=0)
			{
			componentColumnIndices[columnIndices[i]]=i;
			++numComponents;
			}
	
	/* Open the ASCII input file: */
	Misc::File file(fileName,"rt");
	
	/* Read all lines from the input file: */
	char line[256];
	while(!file.eof())
		{
		/* Read the next line from the file: */
		file.gets(line,sizeof(line));
		
		/* Separate and parse the columns: */
		int numParsedColumns=0;
		float componentValues[6];
		for(int i=3;i<6;++i)
			componentValues[i]=1.0f;
		char* columnStart=line;
		for(int columnIndex=0;columnIndex<=maxColumnIndex;++columnIndex)
			{
			/* Find the end of the current column: */
			char* columnEnd;
			for(columnEnd=columnStart;*columnEnd!=','&&*columnEnd!='\n';++columnEnd)
				;
			char terminator=*columnEnd;
			
			if(componentColumnIndices[columnIndex]>=0)
				{
				/* Skip whitespace and quotes at the beginning and the end: */
				char* valueStart=columnStart;
				while(valueStart!=columnEnd&&(isspace(*valueStart)||*valueStart=='\"'))
					++valueStart;
				char* valueEnd=columnEnd;
				while(valueEnd!=valueStart&&(isspace(valueEnd[-1])||valueEnd[-1]=='\"'))
					--valueEnd;
				
				if(valueStart!=valueEnd)
					{
					/* Parse the column's value: */
					*valueEnd='\0';
					char* conversionEnd;
					componentValues[componentColumnIndices[columnIndex]]=strtof(valueStart,&conversionEnd);
					if(conversionEnd!=valueStart)
						++numParsedColumns;
					}
				}
			
			/* Go to the next column: */
			if(terminator=='\n')
				break;
			columnStart=columnEnd+1;
			}
		
		if(numParsedColumns==numComponents)
			{
			/* Store the point position: */
			OctreePoint p;
			for(int i=0;i<3;++i)
				p[i]=componentValues[i];
			
			/* Filter the read RGB color by the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=GLubyte(Math::floor(componentValues[3+i]*colorMask[i]+0.5f));
			p.value[3]=GLubyte(255);
			
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	
	/* Clean up: */
	delete[] componentColumnIndices;
	}

void loadPointFileCsv(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the ASCII input file: */
	Misc::File file(fileName,"rt");
	
	/* Read all lines from the input file: */
	char line[256];
	while(!file.eof())
		{
		/* Read the next line from the file: */
		file.gets(line,sizeof(line));
		
		/* Parse the point coordinates and intensity from the line: */
		OctreePoint p;
		float intensity;
		if(sscanf(line,"%f,%f,%f,%f",&p[0],&p[1],&p[2],&intensity)==4&&intensity!=0.0f)
			{
			/* Convert the read intensity to an RGB color according to the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=GLubyte(Math::floor(intensity*colorMask[i]+0.5f));
			p.value[3]=GLubyte(255);
			
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	}

void loadPointFileCsvRgb(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the ASCII input file: */
	Misc::File file(fileName,"rt");
	
	/* Read all lines from the input file: */
	char line[256];
	while(!file.eof())
		{
		/* Read the next line from the file: */
		file.gets(line,sizeof(line));
		
		/* Parse the point coordinates and intensity from the line: */
		OctreePoint p;
		float rgb[3];
		if(sscanf(line,"%f,%f,%f,%f,%f,%f",&p[0],&p[1],&p[2],&rgb[0],&rgb[1],&rgb[2])==6)
			{
			/* Convert the read RGB values to an RGB color according to the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=GLubyte(Math::floor(rgb[i]*colorMask[i]+0.5f));
			p.value[3]=GLubyte(255);
			
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	}

void loadPointFileAll(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the ASCII input file: */
	Misc::File file(fileName,"rt");
	
	/* Read all lines from the input file: */
	while(!file.eof())
		{
		/* Read the next line from the file: */
		char line[256];
		file.gets(line,sizeof(line));
		
		/* Parse the point coordinates and intensity from the line: */
		float dummy;
		OctreePoint p;
		float intensity;
		if(sscanf(line,"%f %f %f %f %f",&dummy,&p[0],&p[1],&p[2],&intensity)==5&&intensity!=0.0f)
			{
			/* Convert the read intensity to an RGB color according to the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=GLubyte(Math::floor(intensity*colorMask[i]+0.5f));
			p.value[3]=GLubyte(255);
			
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	}

void loadPointFileXyzrgb(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the ASCII input file: */
	Misc::File file(fileName,"rt");
	
	/* Read all lines from the input file: */
	while(!file.eof())
		{
		/* Read the next line from the file: */
		char line[256];
		file.gets(line,sizeof(line));
		
		/* Skip comment lines: */
		if(line[0]=='#')
			continue;
		
		/* Parse the point coordinates and RGB color from the line: */
		OctreePoint p;
		float rgb[3];
		if(sscanf(line,"%f %f %f %f %f %f",&p[0],&p[1],&p[2],&rgb[0],&rgb[1],&rgb[2])==6)
			{
			/* Process the RGB color according to the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=GLubyte(Math::floor(rgb[i]*colorMask[i]+0.5f));
			p.value[3]=GLubyte(255);
			
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	}

void loadPointFilePts(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the ASCII input file: */
	Misc::File file(fileName,"rt");
	
	/* Skip the header line (don't know what it means): */
	char line[256];
	file.gets(line,sizeof(line));
	
	/* Read all lines from the input file: */
	while(!file.eof())
		{
		/* Read the next line from the file: */
		file.gets(line,sizeof(line));
		
		/* Parse the point coordinates and RGB color from the line: */
		OctreePoint p;
		float dummy;
		float rgb[3];
		if(sscanf(line,"%f %f %f %f %f %f %f",&p[0],&p[1],&p[2],&dummy,&rgb[0],&rgb[1],&rgb[2])==7)
			{
			/* Process the RGB color according to the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=GLubyte(Math::floor(rgb[i]*colorMask[i]+0.5f));
			p.value[3]=GLubyte(255);
			
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	}

void loadPointFileTimeXyzi(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the ASCII input file: */
	Misc::File file(fileName,"rt");
	
	/* Read all lines from the input file: */
	while(!file.eof())
		{
		/* Read the next line from the file: */
		char line[256];
		file.gets(line,sizeof(line));
		
		/* Parse the point coordinates and intensity from the line: */
		float time;
		OctreePoint p;
		int returnNumber;
		float dummy;
		float intensity;
		if(sscanf(line,"%f %f %f %f %d %f %f",&time,&p[0],&p[1],&p[2],&returnNumber,&dummy,&intensity)==7&&intensity!=0.0f)
			{
			/* Convert the read intensity to an RGB color according to the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=GLubyte(Math::floor(intensity*colorMask[i]+0.5f));
			p.value[3]=GLubyte(255);
			
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	}

void loadPointFileAsc(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the ASCII input file: */
	Misc::File file(fileName,"rt");
	
	/* Read and otherwise ignore the header line: */
	char line[256];
	file.gets(line,sizeof(line));
	
	/* Read all lines from the input file: */
	while(!file.eof())
		{
		/* Read the next line from the file: */
		file.gets(line,sizeof(line));
		
		/* Parse the point coordinates and intensity from the line: */
		OctreePoint p;
		float date,time;
		int returnNumber,numReturns;
		float offNidar;
		float intensity;
		if(sscanf(line,"%f,%f,%f,%f,%f,%d,%d,%f,%f",&p[0],&p[1],&p[2],&time,&date,&returnNumber,&numReturns,&offNidar,&intensity)==9&&intensity!=0.0f)
			{
			/* Convert the read intensity to an RGB color according to the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=GLubyte(Math::floor(intensity*colorMask[i]+0.5f));
			p.value[3]=GLubyte(255);
			
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	}

void loadPointFileLas(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
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
	
	/* Read all points: */
	file.seekSet(pointDataOffset);
	for(unsigned int i=0;i<numPointRecords;++i)
		{
		/* Read this point: */
		OctreePoint p;
		
		/* Read the point position: */
		int pos[3];
		file.read(pos,3);
		for(int j=0;j<3;++j)
			p[j]=float(double(pos[j])*scale[j]+offset[j]);
		
		/* Read the point intensity: */
		float intensity=float(file.read<unsigned short>());
		for(int j=0;j<3;++j)
			p.value[j]=GLubyte(Math::floor(intensity*colorMask[j]+0.5f));
		p.value[3]=GLubyte(255);
		
		/* Ignore the rest of the point record: */
		char dummy[4];
		file.read(dummy,4);
		file.read<unsigned short>();
		if(pointDataFormat==1)
			file.read<double>();
		
		/* Store the point: */
		opl.addPoint(p);
		}
	}

/**********************************
Helper structure to load IDL files:
**********************************/

struct IDLFileRecord // Structure describing a record in an IDL file
	{
	/* Elements: */
	public:
	unsigned int galID[2]; // Galaxy ID (64-bit integer)
	unsigned int haloID[2]; // Halo ID (64-bit integer)
	unsigned int recordType; // 0=central, 1=sat, 2=orphan
	float position[3];
	float velocity[3];
	float spin[3];
	float ra;
	float dec;
	float z_obs;
	float z;
	float centralMVir;
	float mVir;
	float rVir;
	float vVir;
	float vMax;
	float velDisp;
	float stellarMass;
	float bulgeMass;
	float coldGas;
	float hotGas;
	float ejectedMass;
	float blackHoleMass;
	float sfr;
	float cooling;
	float heating;
	float appMagSdss[5];
	float appMagSdssBulge[5];
	float absMagSdss[5];
	float appMagDeep[4];
	float appMagDeepBulge[4];
	float absMagBvrik[5];
	float absMagBvrikBulge[5];
	float absMagBvrikNoDust[5];
	};

namespace Misc {

template <>
class EndiannessSwapper<IDLFileRecord>
	{
	/* Methods: */
	public:
	static void swap(IDLFileRecord& value)
		{
		swapEndianness(value.galID,2);
		swapEndianness(value.haloID,2);
		swapEndianness(value.recordType);
		swapEndianness(value.position,3);
		swapEndianness(value.velocity,3);
		swapEndianness(value.spin,3);
		swapEndianness(value.ra);
		swapEndianness(value.dec);
		swapEndianness(value.z_obs);
		swapEndianness(value.z);
		swapEndianness(value.centralMVir);
		swapEndianness(value.mVir);
		swapEndianness(value.rVir);
		swapEndianness(value.vVir);
		swapEndianness(value.vMax);
		swapEndianness(value.velDisp);
		swapEndianness(value.stellarMass);
		swapEndianness(value.bulgeMass);
		swapEndianness(value.coldGas);
		swapEndianness(value.hotGas);
		swapEndianness(value.ejectedMass);
		swapEndianness(value.blackHoleMass);
		swapEndianness(value.sfr);
		swapEndianness(value.cooling);
		swapEndianness(value.heating);
		swapEndianness(value.appMagSdss,5);
		swapEndianness(value.appMagSdssBulge,5);
		swapEndianness(value.absMagSdss,5);
		swapEndianness(value.appMagDeep,4);
		swapEndianness(value.appMagDeepBulge,4);
		swapEndianness(value.absMagBvrik,5);
		swapEndianness(value.absMagBvrikBulge,5);
		swapEndianness(value.absMagBvrikNoDust,5);
		};
	};

}

double angdiadistscaled(double z)
	{
	double h0=71.0;
	double oM=0.3;
	double oL=0.7;
	double oR=1.0-(oM+oL);
	
	double sum1=0.0;
	double dz=z/100.0;
	double id=1.0;
	for(int i=0;i<100;++i,id+=dz)
		{
		double ez=Math::sqrt((oM*id+oR)*id*id+oL);
		sum1+=dz/ez;
		}
	double dh=3.0e5/h0;
	double dc=dh*sum1;
	
	if(oR==0.0)
		return dh*sum1;
	else
		{
		double sqrtOR=Math::sqrt(Math::abs(oR));
		return dh*(1.0/sqrtOR)*sinh(sqrtOR*(sum1));
		}
	}

void loadPointFileIdl(OctreePointLoader& opl,const char* fileName,const float colorMask[3])
	{
	/* Open the IDL input file: */
	Misc::LargeFile file(fileName,"rb",Misc::LargeFile::LittleEndian);
	
	/* Read the IDL file header: */
	unsigned int numRecords=file.read<unsigned int>();
	
	/* Read all records: */
	float massMin=Math::Constants<float>::max;
	float massMax=Math::Constants<float>::min;
	float rgbMin[3],rgbMax[3];
	for(int j=0;j<3;++j)
		{
		rgbMin[j]=Math::Constants<float>::max;
		rgbMax[j]=Math::Constants<float>::min;
		}
	for(unsigned int i=0;i<numRecords;++i)
		{
		/* Read the next record: */
		IDLFileRecord record;
		file.read(record);
		
		/* Create a point: */
		OctreePoint p;
		#if 0
		for(int i=0;i<3;++i)
			p[i]=record.position[i];
		#else
		/* New formula using redshift to calculate galaxy position in Cartesian coordinates: */
		double z=3200.0*double(record.z);
		// double z=angdiadistscaled(record.z);
		p[0]=float(double(record.dec)*z);
		p[1]=float(double(record.ra)*z);
		p[2]=float(z);
		#endif
		
		if(massMin>record.mVir)
			massMin=record.mVir;
		if(massMax<record.mVir)
			massMax=record.mVir;
		float rgbFactor=(record.mVir/32565.4f)*0.5f+0.5f;
		
		/* Calculate false color from absolute SDSS magnitudes: */
		float rgb[3];
		rgb[0]=(record.absMagSdss[2]-record.absMagSdss[3]+1.13f)*colorMask[0];
		rgb[1]=((-record.absMagSdss[2]-14.62f)*0.3)*colorMask[1];
		rgb[2]=(record.absMagSdss[1]-record.absMagSdss[2]+0.73f)*colorMask[2];
		for(int j=0;j<3;++j)
			{
			if(rgbMin[j]>rgb[j])
				rgbMin[j]=rgb[j];
			if(rgbMax[j]<rgb[j])
				rgbMax[j]=rgb[j];
			if(rgb[j]<0.5f)
				p.value[j]=GLubyte(0);
			else if(rgb[j]>=254.5f)
				p.value[j]=GLubyte(255);
			else
				p.value[j]=GLubyte(Math::floor(rgb[j]+0.5f));
			}
		p.value[3]=GLubyte(255);
		
		// if(record.recordType==0)
			{
			/* Store the point: */
			opl.addPoint(p);
			}
		}
	
	std::cout<<"mVir range: "<<massMin<<" - "<<massMax<<std::endl;
	std::cout<<"RGB range: "<<rgbMin[0]<<" - "<<rgbMax[0]<<", "<<rgbMin[1]<<" - "<<rgbMax[1]<<", "<<rgbMin[2]<<" - "<<rgbMax[2]<<std::endl;
	}

/**************************************
Helper structures to load octree files:
**************************************/

struct OldLidarOctreeFileHeader // Header structure for old LiDAR octree files
	{
	/* Elements: */
	public:
	Point center; // Center of the cube containing all LiDAR points
	float radius; // Radius (half side length) of cube containing all LiDAR points
	unsigned int maxNumPointsPerNode; // Maximum number of points stored in each octree node
	
	/* Methods: */
	static size_t getSize(void) // Returns the file size of the header
		{
		return sizeof(Point)+sizeof(float)+sizeof(unsigned int);
		};
	void read(Misc::LargeFile& file) // Reads the header from file
		{
		file.read(center.getComponents(),3);
		file.read(radius);
		file.read(maxNumPointsPerNode);
		};
	};

struct OldLidarOctreeFileNode // Structure for a node in an old LiDAR octree file
	{
	/* Elements: */
	public:
	Misc::LargeFile::Offset childrenOffset; // Offset of the node's children in the octree file (or 0 if leaf node)
	unsigned int numPoints; // Number of points in the node
	Misc::LargeFile::Offset pointsOffset; // Offset of the node's points in the points file
	
	/* Methods: */
	static size_t getSize(void) // Returns the file size of an octree node
		{
		return sizeof(Misc::LargeFile::Offset)+sizeof(Misc::LargeFile::Offset)+sizeof(unsigned int);
		};
	void read(Misc::LargeFile& file) // Reads the octree node from file
		{
		file.read(childrenOffset);
		file.read(pointsOffset);
		file.read(numPoints);
		};
	};

void readSubtree(OctreePointLoader& opl,Misc::LargeFile& octFile,Misc::LargeFile& obinFile,const float colorMask[3])
	{
	/* Read the node's header from the octree file: */
	OldLidarOctreeFileNode ofn;
	ofn.read(octFile);
	
	if(ofn.childrenOffset!=0)
		{
		/* Recurse into the node's children: */
		Misc::LargeFile::Offset childOffset=ofn.childrenOffset;
		for(int childIndex=0;childIndex<8;++childIndex,childOffset+=OldLidarOctreeFileNode::getSize())
			{
			octFile.seekSet(childOffset);
			readSubtree(opl,octFile,obinFile,colorMask);
			}
		}
	else if(ofn.numPoints>0)
		{
		/* Read the node's points from the octree point file: */
		obinFile.seekSet(ofn.pointsOffset);
		for(unsigned int i=0;i<ofn.numPoints;++i)
			{
			LidarPoint p;
			obinFile.read<float>(p.getComponents(),3);
			obinFile.read<GLubyte>(p.value.getRgba(),4);
			opl.addPoint(p);
			}
		}
	}

void loadPointFileOctree(OctreePointLoader& opl,const char* fileNameStem,const float colorMask[3])
	{
	/* Open the input octree structure and octree point files: */
	char octFileName[1024],obinFileName[1024];
	snprintf(octFileName,sizeof(octFileName),"%s.oct",fileNameStem);
	snprintf(obinFileName,sizeof(obinFileName),"%s.obin",fileNameStem);
	Misc::LargeFile octFile(octFileName,"rb",Misc::LargeFile::LittleEndian);
	Misc::LargeFile obinFile(obinFileName,"rb",Misc::LargeFile::LittleEndian);
	
	/* Read the octree structure file header: */
	OldLidarOctreeFileHeader ofh;
	ofh.read(octFile);
	
	/* Read all original points from the octree: */
	readSubtree(opl,octFile,obinFile,colorMask);
	}

/**************
Helper classes:
**************/

enum PointFileType // Enumerated type for point file formats
	{
	AUTO, // Tries to autodetect input file format based on file name extension
	BIN, // Simple binary format: number of points followed by (x, y, z, i) tuples
	BINRGB, // Simple binary format: number of points followed by (x, y, z, r, g, b) tuples
	TXT, // Simple ASCII format: rows containing (x, y, z) tuples
	XYZI, // Simple ASCII format: rows containing (x, y, z, i) tuples
	CSV, // Generic CSV format with intensity values; needs column index mapping to work
	CSVRGB, // Generic CSV format with RGB values; needs column index mapping to work
	ALL, // Simple ASCII format: rows containing (GPS time, x, y, z, i) tuples (used in B4 data)
	XYZRGB, // Simple ASCII format: rows containing (x, y, z, r, g, b) tuples
	PTS, // Simple ASCII format: rows containing (x, y, z, something, r, g, b) tuples
	TIMEXYZI, // Simple ASCII format: rows containing (x, y, z, i) tuples and some ignored data columns
	ASC, // ASCII input file in comma-separated format (GEON format?)
	LAS, // Binary interchange format for airborne LiDAR data
	IDL, // Binary interchange format for astronomical applications
	OCTREE, // Point octree file in old file format
	OCTREENEW, // Point octree file in new file format (bypasses temporary octree creation)
	ILLEGAL
	};

int main(int argc,char* argv[])
	{
	/* Parse the command line and load all input files: */
	Misc::Timer loadTimer;
	const char* outputFileName=0;
	unsigned int maxPointsPerNode=1024;
	unsigned int subsamplingFactor=4;
	OctreePointLoader opl;
	const char* tempPointsFileNameTemplate="/tmp/LidarPreprocessorTempPointsXXXXXX";
	PointFileType pointFileType=AUTO;
	float colorMask[3]={1.0f,1.0f,1.0f};
	int csvColumnIndices[6];
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
			else if(strcasecmp(argv[i]+1,"sf")==0)
				{
				++i;
				if(i<argc)
					subsamplingFactor=(unsigned int)(atoi(argv[i]));
				else
					std::cerr<<"Dangling -sf flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"ooc")==0)
				{
				++i;
				if(i<argc)
					opl.setMemorySize(atoi(argv[i]),maxPointsPerNode);
				else
					std::cerr<<"Dangling -ooc flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"to")==0)
				{
				++i;
				if(i<argc)
					opl.setTempFileNameTemplate(argv[i]);
				else
					std::cerr<<"Dangling -to flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"tp")==0)
				{
				++i;
				if(i<argc)
					tempPointsFileNameTemplate=argv[i];
				else
					std::cerr<<"Dangling -tp flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"c")==0)
				{
				if(i+3<argc)
					{
					for(int j=0;j<3;++j)
						{
						++i;
						colorMask[j]=float(atof(argv[i]));
						}
					}
				else
					{
					i=argc;
					std::cerr<<"Dangling -c flag on command line"<<std::endl;
					}
				}
			else if(strcasecmp(argv[i]+1,"auto")==0)
				pointFileType=AUTO;
			else if(strcasecmp(argv[i]+1,"bin")==0)
				pointFileType=BIN;
			else if(strcasecmp(argv[i]+1,"binrgb")==0)
				pointFileType=BINRGB;
			else if(strcasecmp(argv[i]+1,"txt")==0)
				pointFileType=TXT;
			else if(strcasecmp(argv[i]+1,"xyzi")==0)
				pointFileType=XYZI;
			else if(strcasecmp(argv[i]+1,"csv")==0)
				{
				pointFileType=CSV;
				
				/* Read the column index mask: */
				for(int j=0;j<4;++j)
					csvColumnIndices[j]=atoi(argv[i+j+1]);
				i+=4;
				}
			else if(strcasecmp(argv[i]+1,"csvrgb")==0)
				{
				pointFileType=CSVRGB;
				
				/* Read the column index mask: */
				for(int j=0;j<6;++j)
					csvColumnIndices[j]=atoi(argv[i+j+1]);
				i+=6;
				}
			else if(strcasecmp(argv[i]+1,"all")==0)
				pointFileType=ALL;
			else if(strcasecmp(argv[i]+1,"xyzrgb")==0)
				pointFileType=XYZRGB;
			else if(strcasecmp(argv[i]+1,"pts")==0)
				pointFileType=PTS;
			else if(strcasecmp(argv[i]+1,"timexyzi")==0)
				pointFileType=TIMEXYZI;
			else if(strcasecmp(argv[i]+1,"asc")==0)
				pointFileType=ASC;
			else if(strcasecmp(argv[i]+1,"las")==0)
				pointFileType=LAS;
			else if(strcasecmp(argv[i]+1,"idl")==0)
				pointFileType=IDL;
			else if(strcasecmp(argv[i]+1,"oct")==0)
				pointFileType=OCTREE;
			else if(strcasecmp(argv[i]+1,"octnew")==0)
				pointFileType=OCTREENEW;
			else
				std::cerr<<"Unrecognized command line option "<<argv[i]<<std::endl;
			}
		else
			{
			PointFileType thisPointFileType=pointFileType;
			if(thisPointFileType==AUTO)
				{
				/* Find the extension of the input file: */
				const char* extPtr=0;
				for(const char* cPtr=argv[i];*cPtr!='\0';++cPtr)
					if(*cPtr=='.')
						extPtr=cPtr;
				
				/* Determine the file type: */
				if(strcasecmp(extPtr,".bin")==0)
					thisPointFileType=BIN;
				else if(strcasecmp(extPtr,".binrgb")==0)
					thisPointFileType=BINRGB;
				else if(strcasecmp(extPtr,".txt")==0)
					thisPointFileType=TXT;
				else if(strcasecmp(extPtr,".xyzi")==0)
					thisPointFileType=XYZI;
				else if(strcasecmp(extPtr,".all")==0)
					thisPointFileType=ALL;
				else if(strcasecmp(extPtr,".xyzrgb")==0)
					thisPointFileType=XYZRGB;
				else if(strcasecmp(extPtr,".pts")==0)
					thisPointFileType=PTS;
				else if(strcasecmp(extPtr,".xyzi")==0)
					thisPointFileType=XYZI;
				else if(strcasecmp(extPtr,".asc")==0)
					thisPointFileType=ASC;
				else if(strcasecmp(extPtr,".las")==0)
					thisPointFileType=LAS;
				else
					thisPointFileType=ILLEGAL;
				}
			
			switch(thisPointFileType)
				{
				case BIN:
					std::cout<<"Processing binary input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileBin(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case BINRGB:
					std::cout<<"Processing binary input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileBinRgb(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case TXT:
					std::cout<<"Processing ASCII input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileTxt(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case XYZI:
					std::cout<<"Processing ASCII input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileXyzi(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case CSV:
					std::cout<<"Processing CSV input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileGenericCsv(opl,argv[i],csvColumnIndices,colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case CSVRGB:
					std::cout<<"Processing CSV input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileGenericCsvRgb(opl,argv[i],csvColumnIndices,colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case ALL:
					std::cout<<"Processing ASCII input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileAll(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case XYZRGB:
					std::cout<<"Processing ASCII input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileXyzrgb(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case PTS:
					std::cout<<"Processing ASCII input file "<<argv[i]<<"..."<<std::flush;
					loadPointFilePts(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case TIMEXYZI:
					std::cout<<"Processing ASCII input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileTimeXyzi(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case ASC:
					std::cout<<"Processing ASCII input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileAsc(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case LAS:
					std::cout<<"Processing binary input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileLas(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case IDL:
					std::cout<<"Processing binary input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileIdl(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case OCTREE:
					std::cout<<"Processing point octree input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileOctree(opl,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case OCTREENEW:
					std::cout<<"Processing point octree input file "<<argv[i]<<"..."<<std::flush;
					opl.addOctree(argv[i]);
					std::cout<<" done."<<std::endl;
					break;
				
				default:
					std::cerr<<"Input file "<<argv[i]<<" has an unrecognized file name extension"<<std::endl;
				}
			}
		}
	
	/* Check if an output file name was given: */
	if(outputFileName==0)
		{
		std::cerr<<"Usage: "<<argv[0]<<"-o <output file name stem> -np <max points per node> [-c <color mask>] <input file 1> ... <input file n>"<<std::endl;
		return 1;
		}
	
	/* Finish reading points and check if out-of-core processing is required: */
	if(opl.finishReading())
		{
		loadTimer.elapse();
		
		/* Construct an octree with less than maxPointsPerNode points per leaf: */
		Misc::Timer createTimer;
		char tpfnt[1024];
		strcpy(tpfnt,tempPointsFileNameTemplate);
		PointOctreeCreator tree(tpfnt,opl.getTempOctrees(),opl.getMaxNumCachablePoints(),maxPointsPerNode,subsamplingFactor);
		
		/* Delete the temporary point octrees: */
		opl.deleteTempOctrees();
		createTimer.elapse();
		
		/* Write the octree structure and data to two files: */
		Misc::Timer writeTimer;
		char octFileName[1024],obinFileName[1024];
		snprintf(octFileName,sizeof(octFileName),"%s.oct",outputFileName);
		snprintf(obinFileName,sizeof(obinFileName),"%s.obin",outputFileName);
		tree.write(octFileName,obinFileName);
		writeTimer.elapse();
		
		std::cout<<"Time to load input data: "<<loadTimer.getTime()<<"s, time to create octree: "<<createTimer.getTime()<<"s, time to write final octree files: "<<writeTimer.getTime()<<"s"<<std::endl;
		}
	else
		{
		loadTimer.elapse();
		
		/* Construct an octree with less than maxPointsPerNode points per leaf: */
		Misc::Timer createTimer;
		LidarOctreeCreator tree(opl.getPoints(),opl.getNumPoints(),maxPointsPerNode,subsamplingFactor);
		createTimer.elapse();
		
		/* Write the octree structure and data to two files: */
		Misc::Timer writeTimer;
		char octFileName[1024],obinFileName[1024];
		snprintf(octFileName,sizeof(octFileName),"%s.oct",outputFileName);
		snprintf(obinFileName,sizeof(obinFileName),"%s.obin",outputFileName);
		tree.write(octFileName,obinFileName);
		writeTimer.elapse();
		
		std::cout<<"Time to load input data: "<<loadTimer.getTime()<<"s, time to create octree: "<<createTimer.getTime()<<"s, time to write final octree files: "<<writeTimer.getTime()<<"s"<<std::endl;
		}
	
	return 0;
	}
