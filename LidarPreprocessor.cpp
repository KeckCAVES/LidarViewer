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

#include "LidarTypes.h"
#include "PointAccumulator.h"
#include "LidarOctreeCreator.h"

void loadPointFileBin(PointAccumulator& pa,const char* fileName,const float colorMask[3])
	{
	/* Open the binary input file: */
	Misc::LargeFile file(fileName,"rb",Misc::LargeFile::LittleEndian);
	
	/* Read the number of points in the file: */
	size_t numPoints=file.read<unsigned int>();
	
	/* Read all points: */
	for(size_t i=0;i<numPoints;++i)
		{
		/* Read the point position and intensity from the input file: */
		float rp[4];
		file.read(rp,4);
		
		/* Store the point position: */
		LidarPoint p;
		for(int j=0;j<3;++j)
			p[j]=rp[j];
		
		/* Convert the point intensity to and RGB color according to the color mask: */
		for(int j=0;j<3;++j)
			p.value[j]=Color::Scalar(Math::floor(rp[3]*colorMask[j]+0.5f));
		p.value[3]=Color::Scalar(255);
		
		/* Store the point: */
		pa.addPoint(p);
		}
	}

void loadPointFileBinRgb(PointAccumulator& pa,const char* fileName,const float colorMask[3])
	{
	/* Open the binary input file: */
	Misc::LargeFile file(fileName,"rb",Misc::LargeFile::LittleEndian);
	
	/* Read the number of points in the file: */
	size_t numPoints=file.read<unsigned int>();
	
	/* Read all points: */
	for(size_t i=0;i<numPoints;++i)
		{
		/* Read the point position and color from the input file: */
		float rp[3];
		file.read(rp,3);
		Color::Scalar rcol[4];
		file.read(rcol,4);
		
		/* Store the point position: */
		LidarPoint p;
		for(int j=0;j<3;++j)
			p[j]=rp[j];
		
		/* Modify the RGB color according to the color mask: */
		for(int j=0;j<3;++j)
			p.value[j]=Color::Scalar(Math::floor(float(rcol[j])*colorMask[j]+0.5f));
		p.value[3]=Color::Scalar(255);
		
		/* Store the point: */
		pa.addPoint(p);
		}
	}

void loadPointFileLas(PointAccumulator& pa,const char* fileName,const float colorMask[3])
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
	size_t numPointRecords=file.read<unsigned int>();
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
	for(size_t i=0;i<numPointRecords;++i)
		{
		/* Read this point: */
		LidarPoint p;
		
		/* Read the point position: */
		int pos[3];
		file.read(pos,3);
		for(int j=0;j<3;++j)
			p[j]=float(double(pos[j])*scale[j]+offset[j]);
		
		/* Read the point intensity: */
		float intensity=float(file.read<unsigned short>());
		for(int j=0;j<3;++j)
			p.value[j]=Color::Scalar(Math::floor(intensity*colorMask[j]+0.5f));
		p.value[3]=Color::Scalar(255);
		
		/* Ignore the rest of the point record: */
		char dummy[4];
		file.read(dummy,4);
		file.read<unsigned short>();
		if(pointDataFormat==1)
			file.read<double>();
		
		/* Store the point: */
		pa.addPoint(p);
		}
	}

void loadPointFileXyzi(PointAccumulator& pa,const char* fileName,const float colorMask[3])
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
		
		/* Parse the point coordinates and intensity from the line: */
		LidarPoint p;
		float intensity;
		if(sscanf(line,"%f %f %f %f",&p[0],&p[1],&p[2],&intensity)==4)
			{
			/* Convert the read intensity to an RGB color according to the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=Color::Scalar(Math::floor(intensity*colorMask[i]+0.5f));
			p.value[3]=Color::Scalar(255);
			
			/* Store the point: */
			pa.addPoint(p);
			}
		}
	}

void loadPointFileXyzrgb(PointAccumulator& pa,const char* fileName,const float colorMask[3])
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
		LidarPoint p;
		float rgb[3];
		if(sscanf(line,"%f %f %f %f %f %f",&p[0],&p[1],&p[2],&rgb[0],&rgb[1],&rgb[2])==6)
			{
			/* Process the RGB color according to the color mask: */
			for(int i=0;i<3;++i)
				p.value[i]=Color::Scalar(Math::floor(rgb[i]*colorMask[i]+0.5f));
			p.value[3]=Color::Scalar(255);
			
			/* Store the point: */
			pa.addPoint(p);
			}
		}
	}

void loadPointFileGenericASCII(PointAccumulator& pa,const char* fileName,bool strictCsv,bool rgb,const int columnIndices[6],const float colorMask[3])
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
		componentValues[3]=255.0f;
		char* columnStart=line;
		for(int columnIndex=0;columnIndex<=maxColumnIndex;++columnIndex)
			{
			/* Skip whitespace: */
			while(*columnStart!='\0'&&isspace(*columnStart))
				++columnStart;
			
			/* Stop if end of line is reached: */
			if(*columnStart=='\0')
				break;
			
			/* Find the start and end of the current column's value: */
			char* valueStart;
			char* valueEnd;
			char* columnEnd;
			if(*columnStart=='\"')
				{
				/* Process a quoted value: */
				valueStart=columnStart+1;
				for(valueEnd=valueStart;*valueEnd!='\0'&&*valueEnd!='\"';++valueEnd)
					;
				columnEnd=valueEnd;
				if(*columnEnd=='\"')
					++columnEnd;
				}
			else if(strictCsv)
				{
				/* Process an unquoted value according to strict CSV file rules: */
				valueStart=columnStart;
				for(valueEnd=valueStart;*valueEnd!='\0'&&*valueEnd!=',';++valueEnd)
					;
				columnEnd=valueEnd;
				}
			else
				{
				/* Process an unquoted value according to relaxed CSV / space-separated rules: */
				valueStart=columnStart;
				for(valueEnd=valueStart;*valueEnd!='\0'&&*valueEnd!=','&&!isspace(*valueEnd);++valueEnd)
					;
				columnEnd=valueEnd;
				}
			
			/* Process the value if needed: */
			if(componentColumnIndices[columnIndex]>=0)
				{
				/* Process the value and check if anything was processed: */
				char* conversionEnd;
				componentValues[componentColumnIndices[columnIndex]]=strtof(valueStart,&conversionEnd);
				if(conversionEnd!=valueStart)
					++numParsedColumns;
				}
			
			/* Find the end of the column: */
			while(*columnEnd!='\0'&&isspace(*columnEnd))
				++columnEnd;
			
			/* Skip a potential separating comma: */
			if(*columnEnd==',')
				++columnEnd;
			
			/* Go to the next column: */
			columnStart=columnEnd;
			}
		
		/* Check if all required columns have been read: */
		if(numParsedColumns==numComponents)
			{
			/* Store the point position: */
			LidarPoint p;
			for(int i=0;i<3;++i)
				p[i]=componentValues[i];
			
			if(rgb)
				{
				/* Modify the read RGB color according to the color mask: */
				for(int i=0;i<3;++i)
					{
					p.value[i]=Color::Scalar(Math::floor(componentValues[3+i]*colorMask[i]+0.5f));
					if(intensityMin>componentValues[3+i])
						intensityMin=componentValues[3+i];
					if(intensityMax<componentValues[3+i])
						intensityMax=componentValues[3+i];
					}
				}
			else
				{
				/* Convert the read intensity to an RGB color according to the color mask: */
				for(int i=0;i<3;++i)
					p.value[i]=Color::Scalar(Math::floor(componentValues[3]*colorMask[i]+0.5f));
				
				if(intensityMin>componentValues[3])
					intensityMin=componentValues[3];
				if(intensityMax<componentValues[3])
					intensityMax=componentValues[3];
				}
			p.value[3]=Color::Scalar(255);
			
			/* Store the point: */
			pa.addPoint(p);
			}
		}
	if(rgb)
		std::cout<<"Point data RGB component range: "<<intensityMin<<", "<<intensityMax<<std::endl;
	else
		std::cout<<"Point data intensity range: "<<intensityMin<<", "<<intensityMax<<std::endl;
	
	/* Clean up: */
	delete[] componentColumnIndices;
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
		}
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

void loadPointFileIdl(PointAccumulator& pa,const char* fileName,const float colorMask[3])
	{
	/* Open the IDL input file: */
	Misc::LargeFile file(fileName,"rb",Misc::LargeFile::LittleEndian);
	
	/* Read the IDL file header: */
	size_t numRecords=file.read<unsigned int>();
	
	/* Read all records: */
	float massMin=Math::Constants<float>::max;
	float massMax=Math::Constants<float>::min;
	float rgbMin[3],rgbMax[3];
	for(int j=0;j<3;++j)
		{
		rgbMin[j]=Math::Constants<float>::max;
		rgbMax[j]=Math::Constants<float>::min;
		}
	for(size_t i=0;i<numRecords;++i)
		{
		/* Read the next record: */
		IDLFileRecord record;
		file.read(record);
		
		/* Create a point: */
		LidarPoint p;
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
				p.value[j]=Color::Scalar(0);
			else if(rgb[j]>=254.5f)
				p.value[j]=Color::Scalar(255);
			else
				p.value[j]=Color::Scalar(Math::floor(rgb[j]+0.5f));
			}
		p.value[3]=Color::Scalar(255);
		
		// if(record.recordType==0)
			{
			/* Store the point: */
			pa.addPoint(p);
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
	Scalar radius; // Radius (half side length) of cube containing all LiDAR points
	unsigned int maxNumPointsPerNode; // Maximum number of points stored in each octree node
	
	/* Constructors and destructors: */
	OldLidarOctreeFileHeader(Misc::LargeFile& file) // Reads the header from file
		{
		file.read(center.getComponents(),3);
		file.read(radius);
		file.read(maxNumPointsPerNode);
		}
	
	/* Methods: */
	static size_t getSize(void) // Returns the file size of the header
		{
		return sizeof(Point)+sizeof(Scalar)+sizeof(unsigned int);
		}
	};

struct OldLidarOctreeFileNode // Structure for a node in a LiDAR octree file
	{
	/* Elements: */
	public:
	Misc::LargeFile::Offset childrenOffset; // Offset of the node's children in the octree file (or 0 if leaf node)
	Scalar detailSize; // Detail size of this node, for proper LOD computation
	unsigned int numPoints; // Number of points in the node
	Misc::LargeFile::Offset pointsOffset; // Offset of the node's points in the points file
	
	/* Constructors and destructors: */
	OldLidarOctreeFileNode(Misc::LargeFile& file) // Reads the octree node from file
		{
		file.read(childrenOffset);
		file.read(detailSize);
		file.read(pointsOffset);
		file.read(numPoints);
		}
	
	/* Methods: */
	static size_t getSize(void) // Returns the file size of an octree node
		{
		return sizeof(Misc::LargeFile::Offset)+sizeof(Scalar)+sizeof(Misc::LargeFile::Offset)+sizeof(unsigned int);
		}
	};

void readOctreeFileSubtree(PointAccumulator& pa,Misc::LargeFile& octFile,Misc::LargeFile& obinFile,const float colorMask[3])
	{
	/* Read the node's header from the octree file: */
	OldLidarOctreeFileNode ofn(octFile);
	
	if(ofn.childrenOffset!=0)
		{
		/* Recurse into the node's children: */
		Misc::LargeFile::Offset childOffset=ofn.childrenOffset;
		for(int childIndex=0;childIndex<8;++childIndex,childOffset+=OldLidarOctreeFileNode::getSize())
			{
			octFile.seekSet(childOffset);
			readOctreeFileSubtree(pa,octFile,obinFile,colorMask);
			}
		}
	else if(ofn.numPoints>0)
		{
		/* Read the node's points from the octree point file: */
		obinFile.seekSet(ofn.pointsOffset);
		for(unsigned int i=0;i<ofn.numPoints;++i)
			{
			LidarPoint p;
			obinFile.read(p.getComponents(),3);
			obinFile.read(p.value.getRgba(),4);
			pa.addPoint(p);
			}
		}
	}

void loadPointFileOctree(PointAccumulator& pa,const char* fileNameStem,const float colorMask[3])
	{
	/* Open the input octree structure and octree point files: */
	char octFileName[1024];
	snprintf(octFileName,sizeof(octFileName),"%s.oct",fileNameStem);
	Misc::LargeFile octFile(octFileName,"rb",Misc::LargeFile::LittleEndian);
	char obinFileName[1024];
	snprintf(obinFileName,sizeof(obinFileName),"%s.obin",fileNameStem);
	Misc::LargeFile obinFile(obinFileName,"rb",Misc::LargeFile::LittleEndian);
	
	/* Read the octree structure file header: */
	OldLidarOctreeFileHeader ofh(octFile);
	
	/* Read all original points from the octree: */
	readOctreeFileSubtree(pa,octFile,obinFile,colorMask);
	}

/**************
Helper classes:
**************/

enum PointFileType // Enumerated type for point file formats
	{
	AUTO, // Tries to autodetect input file format based on file name extension
	BIN, // Simple binary format: number of points followed by (x, y, z, i) tuples
	BINRGB, // Simple binary format: number of points followed by (x, y, z, r, g, b) tuples
	LAS, // Standard binary LiDAR point cloud interchange format
	XYZI, // Simple ASCII format: rows containing space-separated (x, y, z, i) tuples
	XYZRGB, // Simple ASCII format: rows containing space-separated (x, y, z, r, g, b) tuples
	ASCII, // Generic ASCII format with intensity data
	ASCIIRGB, // Generic ASCII format with RGB data
	CSV, // Strict comma-separated values file with intensity data
	CSVRGB, // Strict comma-separated values file with RGB data
	IDL, // Redshift data file in IDL format
	OCTREE, // Point octree file in old file format
	ILLEGAL
	};

/****************
Helper functions:
****************/

bool readColumnIndexMask(int argc,char* argv[],int& argi,int columnIndices[6])
	{
	/* Initialize the column indices to invalid: */
	for(int i=0;i<6;++i)
		columnIndices[i]=-1;
	
	/* Read all parameters that are integers: */
	for(int i=0;i<6;++i)
		{
		++argi;
		if(argi<argc)
			{
			char* conversionEnd;
			int value=int(strtol(argv[argi],&conversionEnd,10));
			if(*conversionEnd=='\0')
				columnIndices[i]=value;
			else
				{
				--argi;
				break;
				}
			}
		}
	
	/* Check if the index array at least contains x, y, z components: */
	return columnIndices[0]>=0&&columnIndices[1]>=0&&columnIndices[2]>=0;
	}

int main(int argc,char* argv[])
	{
	/* Parse the command line and load all input files: */
	Misc::Timer loadTimer;
	const char* outputFileName=0;
	unsigned int maxPointsPerNode=1024;
	PointAccumulator pa;
	const char* tempPointsFileNameTemplate="/tmp/LidarPreprocessorTempPointsXXXXXX";
	PointFileType pointFileType=AUTO;
	float colorMask[3]={1.0f,1.0f,1.0f};
	int asciiColumnIndices[6];
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
			else if(strcasecmp(argv[i]+1,"ooc")==0)
				{
				++i;
				if(i<argc)
					pa.setMemorySize(atoi(argv[i]),pa.getMaxNumPointsPerNode());
				else
					std::cerr<<"Dangling -ooc flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"to")==0)
				{
				++i;
				if(i<argc)
					pa.setTempOctreeFileNameTemplate(argv[i]);
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
			else if(strcasecmp(argv[i]+1,"las")==0)
				pointFileType=LAS;
			else if(strcasecmp(argv[i]+1,"xyzi")==0)
				pointFileType=XYZI;
			else if(strcasecmp(argv[i]+1,"xyzrgb")==0)
				pointFileType=XYZRGB;
			else if(strcasecmp(argv[i]+1,"ascii")==0)
				{
				pointFileType=ASCII;
				
				/* Read the column index mask: */
				if(!readColumnIndexMask(argc,argv,i,asciiColumnIndices))
					{
					std::cerr<<"Invalid column indices for ASCII file"<<std::endl;
					pointFileType==ILLEGAL;
					}
				}
			else if(strcasecmp(argv[i]+1,"asciirgb")==0)
				{
				pointFileType=ASCIIRGB;
				
				/* Read the column index mask: */
				if(!readColumnIndexMask(argc,argv,i,asciiColumnIndices))
					{
					std::cerr<<"Invalid column indices for RGB ASCII file"<<std::endl;
					pointFileType==ILLEGAL;
					}
				}
			else if(strcasecmp(argv[i]+1,"csv")==0)
				{
				pointFileType=CSV;
				
				/* Read the column index mask: */
				if(!readColumnIndexMask(argc,argv,i,asciiColumnIndices))
					{
					std::cerr<<"Invalid column indices for CSV file"<<std::endl;
					pointFileType==ILLEGAL;
					}
				}
			else if(strcasecmp(argv[i]+1,"csvrgb")==0)
				{
				pointFileType=CSVRGB;
				
				/* Read the column index mask: */
				if(!readColumnIndexMask(argc,argv,i,asciiColumnIndices))
					{
					std::cerr<<"Invalid column indices for RGB CSV file"<<std::endl;
					pointFileType==ILLEGAL;
					}
				}
			else if(strcasecmp(argv[i]+1,"idl")==0)
				pointFileType=IDL;
			else if(strcasecmp(argv[i]+1,"oct")==0)
				pointFileType=OCTREE;
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
				else if(strcasecmp(extPtr,".las")==0)
					thisPointFileType=LAS;
				else if(strcasecmp(extPtr,".xyzi")==0)
					thisPointFileType=XYZI;
				else if(strcasecmp(extPtr,".xyzrgb")==0)
					thisPointFileType=XYZI;
				else if(strcasecmp(extPtr,".oct")==0)
					thisPointFileType=OCTREE;
				else
					thisPointFileType=ILLEGAL;
				}
			
			switch(thisPointFileType)
				{
				case BIN:
					std::cout<<"Processing binary input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileBin(pa,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case BINRGB:
					std::cout<<"Processing RGB binary input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileBinRgb(pa,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case LAS:
					std::cout<<"Processing binary input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileLas(pa,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case XYZI:
					std::cout<<"Processing XYZI input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileXyzi(pa,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case XYZRGB:
					std::cout<<"Processing XYZRGB input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileXyzrgb(pa,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case ASCII:
					std::cout<<"Processing generic ASCII input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileGenericASCII(pa,argv[i],false,false,asciiColumnIndices,colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case ASCIIRGB:
					std::cout<<"Processing generic RGB ASCII input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileGenericASCII(pa,argv[i],false,true,asciiColumnIndices,colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case CSV:
					std::cout<<"Processing generic CSV input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileGenericASCII(pa,argv[i],true,false,asciiColumnIndices,colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case CSVRGB:
					std::cout<<"Processing generic RGB CSV input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileGenericASCII(pa,argv[i],true,true,asciiColumnIndices,colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case IDL:
					std::cout<<"Processing redshift IDL input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileIdl(pa,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				case OCTREE:
					std::cout<<"Processing LiDAR octree input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileOctree(pa,argv[i],colorMask);
					std::cout<<" done."<<std::endl;
					break;
				
				default:
					std::cerr<<"Input file "<<argv[i]<<" has an unrecognized file format"<<std::endl;
				}
			}
		}
	
	/* Check if an output file name was given: */
	if(outputFileName==0)
		{
		std::cerr<<"Usage: "<<argv[0]<<"-o <output file name stem> -np <max points per node> <input file spec 1> ... <input file spec n>"<<std::endl;
		std::cerr<<"Input file spec: [-c <red> <green> <blue>] <format spec> <file name>"<<std::endl;
		std::cerr<<"Format spec: -AUTO"<<std::endl;
		std::cerr<<"             -BIN"<<std::endl;
		std::cerr<<"             -BINRGB"<<std::endl;
		std::cerr<<"             -LAS"<<std::endl;
		std::cerr<<"             -XYZI"<<std::endl;
		std::cerr<<"             -XYZRGB"<<std::endl;
		std::cerr<<"             -ASCII <x column> <y column> <z column> [<intensity column>]"<<std::endl;
		std::cerr<<"             -ASCIIRGB <x column> <y column> <z column> [<r column> <g column> <b column>]"<<std::endl;
		std::cerr<<"             -CSV <x column> <y column> <z column> [<intensity column>]"<<std::endl;
		std::cerr<<"             -CSVRGB <x column> <y column> <z column> [<r column> <g column> <b column>]"<<std::endl;
		std::cerr<<"             -IDL"<<std::endl;
		std::cerr<<"             -OCTREE"<<std::endl;
		return 1;
		}
	
	/* Finish reading points: */
	pa.finishReading();
	loadTimer.elapse();
	
	/* Construct an octree with less than maxPointsPerNode points per leaf: */
	Misc::Timer createTimer;
	LidarOctreeCreator tree(pa.getMaxNumCacheablePoints(),maxPointsPerNode,pa.getTempOctrees(),tempPointsFileNameTemplate);
	
	/* Delete the temporary point octrees: */
	pa.deleteTempOctrees();
	createTimer.elapse();
	
	/* Write the octree structure and data to the destination LiDAR file: */
	Misc::Timer writeTimer;
	tree.write(outputFileName);
	writeTimer.elapse();

	std::cout<<"Time to load input data: "<<loadTimer.getTime()<<"s, time to create octree: "<<createTimer.getTime()<<"s, time to write final octree files: "<<writeTimer.getTime()<<"s"<<std::endl;
	
	return 0;
	}
