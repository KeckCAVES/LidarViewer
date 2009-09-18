/***********************************************************************
New version of LiDAR data preprocessor.
Copyright (c) 2005-2009 Oliver Kreylos

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
#include <iomanip>
#include <string>
#include <vector>
#include <Misc/ThrowStdErr.h>
#include <Misc/Endianness.h>
#include <Misc/FileNameExtensions.h>
#include <Misc/File.h>
#include <Misc/LargeFile.h>
#include <Misc/FileCharacterSource.h>
#include <Threads/GzippedFileCharacterSource.h>
#include <Misc/ValueSource.h>
#include <Misc/Timer.h>
#include <Misc/StandardValueCoders.h>
#include <Misc/ConfigurationFile.h>
#include <Math/Math.h>
#include <Math/Constants.h>
#include <Geometry/OrthonormalTransformation.h>
#include <Geometry/GeometryValueCoders.h>

#include "LidarTypes.h"
#include "PointAccumulator.h"
#include "LidarProcessOctree.h"
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
		
		/* Convert the point intensity to an RGB color according to the color mask: */
		for(int j=0;j<3;++j)
			p.value[j]=Color::clampRound(rp[3]*colorMask[j]);
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
			p.value[j]=Color::clampRound(rcol[j]*colorMask[j]);
		p.value[3]=Color::Scalar(255);
		
		/* Store the point: */
		pa.addPoint(p);
		}
	}

void loadPointFileLas(PointAccumulator& pa,const char* fileName,const float colorMask[3],const double pointOffset[3])
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
	
	/* Check the LAS file header for sanity: */
	switch(pointDataFormat)
		{
		case 0:
			if(pointDataRecordLength!=20)
				return;
			break;
		
		case 1:
			if(pointDataRecordLength!=28)
				return;
			break;
		
		case 2:
			if(pointDataRecordLength!=26)
				return;
			break;
		
		case 3:
			if(pointDataRecordLength!=34)
				return;
			break;
		
		default:
			return;
		}
	
	/* Apply the LAS offset: */
	for(int i=0;i<3;++i)
		offset[i]+=pointOffset[i];
	
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
		
		/* Skip irrelevant information: */
		char dummy[4];
		file.read(dummy,4);
		file.read<unsigned short>();
		if(pointDataFormat&0x1)
			file.read<double>();
		if(pointDataFormat>=2)
			{
			/* Assign point color from stored RGB data: */
			unsigned short rgb[3];
			file.read(rgb,3);
			for(int j=0;j<3;++j)
				p.value[j]=Color::clampRound(rgb[j]*colorMask[j]);
			}
		else
			{
			/* Assign point color from intensity: */
			for(int j=0;j<3;++j)
				p.value[j]=Color::clampRound(intensity*colorMask[j]);
			}
		p.value[3]=Color::Scalar(255);
		
		/* Store the point: */
		pa.addPoint(p);
		}
	}

void loadPointFileXyzi(PointAccumulator& pa,const char* fileName,const float colorMask[3],const double pointOffset[3])
	{
	/* Open the ASCII input file: */
	Threads::GzippedFileCharacterSource file(fileName);
	Misc::ValueSource reader(file);
	reader.setPunctuation('#',true);
	reader.setPunctuation('\n',true);
	reader.skipWs();
	
	/* Read all lines from the input file: */
	LidarPoint::Scalar pMin[3],pMax[3];
	for(int i=0;i<3;++i)
		{
		pMin[i]=Math::Constants<LidarPoint::Scalar>::max;
		pMax[i]=Math::Constants<LidarPoint::Scalar>::min;
		}
	float intensityMin=Math::Constants<float>::max;
	float intensityMax=Math::Constants<float>::min;
	size_t lineNumber=1;
	while(!reader.eof())
		{
		/* Check for comment or empty lines: */
		if(reader.peekc()!='#'&&reader.peekc()!='\n')
			{
			try
				{
				/* Parse the point coordinates and intensity from the line: */
				LidarPoint p;
				for(int i=0;i<3;++i)
					{
					p[i]=LidarPoint::Scalar(reader.readNumber()+pointOffset[i]);
					if(pMin[i]>p[i])
						pMin[i]=p[i];
					if(pMax[i]<p[i])
						pMax[i]=p[i];
					}
				float intensity=float(reader.readNumber());
				if(intensityMin>intensity)
					intensityMin=intensity;
				if(intensityMax<intensity)
					intensityMax=intensity;
				
				/* Convert the read intensity to an RGB color according to the color mask: */
				for(int i=0;i<3;++i)
					p.value[i]=Color::clampRound(intensity*colorMask[i]);
				p.value[3]=Color::Scalar(255);
				
				/* Store the point: */
				pa.addPoint(p);
				}
			catch(Misc::ValueSource::NumberError err)
				{
				std::cerr<<"loadPointFileXyzi: Point parsing error in line "<<lineNumber<<std::endl;
				}
			}
		
		/* Skip the rest of the line: */
		reader.skipLine();
		++lineNumber;
		reader.skipWs();
		}
	
	std::cout.setf(std::ios::fixed);
	std::cout<<std::setprecision(6);
	std::cout<<"Bounding box: ["<<pMin[0]<<", "<<pMax[0]<<"] x ["<<pMin[1]<<", "<<pMax[1]<<"] x ["<<pMin[2]<<", "<<pMax[2]<<"]"<<std::endl;
	std::cout<<"Intensity range: "<<intensityMin<<", "<<intensityMax<<std::endl;
	}

void loadPointFileXyzrgb(PointAccumulator& pa,const char* fileName,const float colorMask[3],const double pointOffset[3])
	{
	/* Open the ASCII input file: */
	Threads::GzippedFileCharacterSource file(fileName);
	Misc::ValueSource reader(file);
	reader.setPunctuation('#',true);
	reader.setPunctuation('\n',true);
	reader.skipWs();
	
	/* Read all lines from the input file: */
	LidarPoint::Scalar pMin[3],pMax[3];
	float rgbMin[3],rgbMax[3];
	for(int i=0;i<3;++i)
		{
		pMin[i]=Math::Constants<LidarPoint::Scalar>::max;
		pMax[i]=Math::Constants<LidarPoint::Scalar>::min;
		rgbMin[i]=Math::Constants<float>::max;
		rgbMax[i]=Math::Constants<float>::min;
		}
	size_t lineNumber=1;
	while(!reader.eof())
		{
		/* Check for comment or empty lines: */
		if(reader.peekc()!='#'&&reader.peekc()!='\n')
			{
			try
				{
				/* Parse the point coordinates and color from the line: */
				LidarPoint p;
				for(int i=0;i<3;++i)
					{
					p[i]=LidarPoint::Scalar(reader.readNumber()+pointOffset[i]);
					if(pMin[i]>p[i])
						pMin[i]=p[i];
					if(pMax[i]<p[i])
						pMax[i]=p[i];
					}
				for(int i=0;i<3;++i)
					{
					float col=float(reader.readNumber())*colorMask[i];
					if(rgbMin[i]<col)
						rgbMin[i]=col;
					if(rgbMax[i]>col)
						rgbMax[i]=col;
					p.value[i]=Color::clampRound(col);
					}
				p.value[3]=Color::Scalar(255);
				
				/* Store the point: */
				pa.addPoint(p);
				}
			catch(Misc::ValueSource::NumberError err)
				{
				std::cerr<<"loadPointFileXyzrgb: Point parsing error in line "<<lineNumber<<std::endl;
				}
			}
		
		/* Skip the rest of the line: */
		reader.skipLine();
		++lineNumber;
		reader.skipWs();
		}
	
	std::cout.setf(std::ios::fixed);
	std::cout<<std::setprecision(6);
	std::cout<<"Bounding box: ["<<pMin[0]<<", "<<pMax[0]<<"] x ["<<pMin[1]<<", "<<pMax[1]<<"] x ["<<pMin[2]<<", "<<pMax[2]<<"]"<<std::endl;
	std::cout<<"RGB range: ["<<rgbMin[0]<<", "<<rgbMax[0]<<"] x ["<<rgbMin[1]<<", "<<rgbMax[1]<<"] x ["<<rgbMin[2]<<", "<<rgbMax[2]<<"]"<<std::endl;
	}

void loadPointFileGenericASCII(PointAccumulator& pa,const char* fileName,int numHeaderLines,bool strictCsv,bool rgb,const int columnIndices[6],const float colorMask[3])
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
	double* componentValues=new double[numComponents];
	
	/* Open the ASCII input file: */
	Threads::GzippedFileCharacterSource file(fileName);
	Misc::ValueSource reader(file);
	if(strictCsv)
		reader.setWhitespace("");
	reader.setPunctuation("#,\n");
	reader.skipWs();
	size_t lineNumber=1;
	
	/* Skip the header lines: */
	for(int i=0;i<numHeaderLines;++i)
		{
		reader.skipLine();
		reader.skipWs();
		++lineNumber;
		}
	
	/* Read all lines from the input file: */
	float intensityMin=Math::Constants<float>::max;
	float intensityMax=Math::Constants<float>::min;
	try
		{
		while(!reader.eof())
			{
			/* Check for comments or empty lines: */
			if(reader.peekc()!='#'&&reader.peekc()!='\n')
				{
				/* Read all columns: */
				for(int columnIndex=0;columnIndex<=maxColumnIndex;++columnIndex)
					{
					if(componentColumnIndices[columnIndex]>=0)
						{
						/* Read the next value: */
						componentValues[componentColumnIndices[columnIndex]]=reader.readNumber();
						}
					else
						reader.skipString();

					if(columnIndex<maxColumnIndex)
						{
						/* Check for separator: */
						if(reader.peekc()==',')
							reader.skipString();
						}
					}

				/* Store the point position: */
				LidarPoint p;
				for(int i=0;i<3;++i)
					p[i]=componentValues[i];

				if(rgb)
					{
					/* Modify the read RGB color according to the color mask: */
					for(int i=0;i<3;++i)
						{
						p.value[i]=Color::clampRound(componentValues[3+i]*colorMask[i]);
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
						p.value[i]=Color::clampRound(componentValues[3]*colorMask[i]);

					if(intensityMin>componentValues[3])
						intensityMin=componentValues[3];
					if(intensityMax<componentValues[3])
						intensityMax=componentValues[3];
					}
				p.value[3]=Color::Scalar(255);

				/* Store the point: */
				pa.addPoint(p);
				}

			/* Skip to the next line: */
			reader.skipLine();
			reader.skipWs();
			++lineNumber;
			}
		}
	catch(std::runtime_error err)
		{
		std::cerr<<"Caught exception "<<err.what()<<" in line "<<lineNumber<<" in input file "<<fileName<<std::endl;
		}
	
	if(rgb)
		std::cout<<"Point data RGB component range: "<<intensityMin<<", "<<intensityMax<<std::endl;
	else
		std::cout<<"Point data intensity range: "<<intensityMin<<", "<<intensityMax<<std::endl;
	
	/* Clean up: */
	delete[] componentColumnIndices;
	delete[] componentValues;
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
			p.value[j]=Color::clampRound(rgb[j]);
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

/*****************************************************
Helper class to load LiDAR files in new octree format:
*****************************************************/

class LidarFileLoader
	{
	/* Elements: */
	private:
	PointAccumulator& pa;
	float colorMask[3];
	
	/* Constructors and destructors: */
	public:
	LidarFileLoader(PointAccumulator& sPa,const float sColorMask[3])
		:pa(sPa)
		{
		for(int i=0;i<3;++i)
			colorMask[i]=sColorMask[i];
		}
	
	/* Methods: */
	void operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel)
		{
		/* Check if this node is a leaf: */
		if(node.isLeaf())
			{
			/* Add each node point to the point accumulator: */
			for(unsigned int i=0;i<node.getNumPoints();++i)
				pa.addPoint(node[i]);
			}
		}
	};


void loadLidarFile(PointAccumulator& pa,const char* lidarFileName,const float colorMask[3])
	{
	/* Open the LiDAR file: */
	LidarProcessOctree lpo(lidarFileName,size_t(64)*size_t(1024)*size_t(1024));
	LidarFileLoader lfl(pa,colorMask);
	lpo.processNodesPostfix(lfl);
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
	LIDAR, // LiDAR file in new octree file format
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
	/* Set default values for all parameters: */
	unsigned int memoryCacheSize=512;
	unsigned int tempOctreeMaxNumPointsPerNode=4096;
	std::string tempOctreeFileNameTemplate="/tmp/LidarPreprocessorTempOctree";
	unsigned int maxNumPointsPerNode=4096;
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
		tempPointFileNameTemplate=cfg.retrieveValue<std::string>("./tempPointFileNameTemplate",tempPointFileNameTemplate);
		}
	catch(std::runtime_error err)
		{
		/* Just ignore the error */
		}
	
	/* Initialize transient parameters: */
	const char* outputFileName=0;
	PointFileType pointFileType=AUTO;
	float colorMask[3]={1.0f,1.0f,1.0f};
	int asciiColumnIndices[6];
	double lasOffset[3]={0.0,0.0,0.0};
	int numHeaderLines=0;
	bool havePoints=false;
	
	/* Parse the command line and load all input files: */
	Misc::Timer loadTimer;
	PointAccumulator pa;
	pa.setMemorySize(memoryCacheSize,tempOctreeMaxNumPointsPerNode);
	pa.setTempOctreeFileNameTemplate(tempOctreeFileNameTemplate+"XXXXXX");
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
			else if(strcasecmp(argv[i]+1,"ooc")==0)
				{
				++i;
				if(i<argc)
					{
					memoryCacheSize=(unsigned int)(atoi(argv[i]));
					pa.setMemorySize(memoryCacheSize,tempOctreeMaxNumPointsPerNode);
					}
				else
					std::cerr<<"Dangling -ooc flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"to")==0)
				{
				++i;
				if(i<argc)
					{
					if(!havePoints)
						{
						tempOctreeFileNameTemplate=argv[i];
						pa.setTempOctreeFileNameTemplate(tempOctreeFileNameTemplate+"XXXXXX");
						}
					else
						std::cerr<<"Ignoring -to flag; must be specified before any input point sets are read"<<std::endl;
					}
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
			else if(strcasecmp(argv[i]+1,"transform")==0)
				{
				++i;
				if(i<argc)
					{
					/* Set the point accumulator's current transformation: */
					PointAccumulator::ONTransform transform=Misc::ValueCoder<PointAccumulator::ONTransform>::decode(argv[i],argv[i]+strlen(argv[i]),0);
					pa.setTransform(transform);
					}
				else
					std::cerr<<"Dangling -transform flag on command line"<<std::endl;
				}
			else if(strcasecmp(argv[i]+1,"auto")==0)
				pointFileType=AUTO;
			else if(strcasecmp(argv[i]+1,"bin")==0)
				pointFileType=BIN;
			else if(strcasecmp(argv[i]+1,"binrgb")==0)
				pointFileType=BINRGB;
			else if(strcasecmp(argv[i]+1,"las")==0)
				pointFileType=LAS;
			else if(strcasecmp(argv[i]+1,"lasoffset")==0)
				{
				for(int j=0;j<3;++j)
					{
					++i;
					lasOffset[j]=atof(argv[i]);
					}
				}
			else if(strcasecmp(argv[i]+1,"header")==0)
				{
				++i;
				numHeaderLines=atof(argv[i]);
				}
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
			else if(strcasecmp(argv[i]+1,"lidar")==0)
				pointFileType=LIDAR;
			else
				std::cerr<<"Unrecognized command line option "<<argv[i]<<std::endl;
			}
		else
			{
			PointFileType thisPointFileType=pointFileType;
			if(thisPointFileType==AUTO)
				{
				/* Find the extension of the input file: */
				const char* extPtr=Misc::getExtension(argv[i]);
				
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
				else if(strcasecmp(extPtr,".LiDAR")==0)
					thisPointFileType=LIDAR;
				else
					thisPointFileType=ILLEGAL;
				}
			
			switch(thisPointFileType)
				{
				case BIN:
					std::cout<<"Processing binary input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileBin(pa,argv[i],colorMask);
					havePoints=true;
					std::cout<<" done."<<std::endl;
					break;
				
				case BINRGB:
					std::cout<<"Processing RGB binary input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileBinRgb(pa,argv[i],colorMask);
					havePoints=true;
					std::cout<<" done."<<std::endl;
					break;
				
				case LAS:
					std::cout<<"Processing binary input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileLas(pa,argv[i],colorMask,lasOffset);
					havePoints=true;
					std::cout<<" done."<<std::endl;
					break;
				
				case XYZI:
					std::cout<<"Processing XYZI input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileXyzi(pa,argv[i],colorMask,lasOffset);
					havePoints=true;
					std::cout<<" done."<<std::endl;
					break;
				
				case XYZRGB:
					std::cout<<"Processing XYZRGB input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileXyzrgb(pa,argv[i],colorMask,lasOffset);
					havePoints=true;
					std::cout<<" done."<<std::endl;
					break;
				
				case ASCII:
					std::cout<<"Processing generic ASCII input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileGenericASCII(pa,argv[i],numHeaderLines,false,false,asciiColumnIndices,colorMask);
					havePoints=true;
					std::cout<<" done."<<std::endl;
					break;
				
				case ASCIIRGB:
					std::cout<<"Processing generic RGB ASCII input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileGenericASCII(pa,argv[i],numHeaderLines,false,true,asciiColumnIndices,colorMask);
					havePoints=true;
					std::cout<<" done."<<std::endl;
					break;
				
				case CSV:
					std::cout<<"Processing generic CSV input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileGenericASCII(pa,argv[i],numHeaderLines,true,false,asciiColumnIndices,colorMask);
					havePoints=true;
					std::cout<<" done."<<std::endl;
					break;
				
				case CSVRGB:
					std::cout<<"Processing generic RGB CSV input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileGenericASCII(pa,argv[i],numHeaderLines,true,true,asciiColumnIndices,colorMask);
					havePoints=true;
					std::cout<<" done."<<std::endl;
					break;
				
				case IDL:
					std::cout<<"Processing redshift IDL input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileIdl(pa,argv[i],colorMask);
					havePoints=true;
					std::cout<<" done."<<std::endl;
					break;
				
				case OCTREE:
					std::cout<<"Processing LiDAR octree input file "<<argv[i]<<"..."<<std::flush;
					loadPointFileOctree(pa,argv[i],colorMask);
					havePoints=true;
					std::cout<<" done."<<std::endl;
					break;
				
				case LIDAR:
					std::cout<<"Processing LiDAR input file "<<argv[i]<<"..."<<std::flush;
					loadLidarFile(pa,argv[i],colorMask);
					havePoints=true;
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
		std::cerr<<"Input file spec: [-c <red> <green> <blue>] [-header <number of header lines>] <format spec> <file name>"<<std::endl;
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
		std::cerr<<"             -OCT"<<std::endl;
		return 1;
		}
	
	/* Finish reading points: */
	pa.finishReading();
	loadTimer.elapse();
	
	/* Construct an octree with less than maxPointsPerNode points per leaf: */
	Misc::Timer createTimer;
	LidarOctreeCreator tree(pa.getMaxNumCacheablePoints(),maxNumPointsPerNode,pa.getTempOctrees(),tempPointFileNameTemplate+"XXXXXX");
	
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
