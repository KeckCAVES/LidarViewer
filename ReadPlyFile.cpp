/***********************************************************************
ReadPlyFile - Function to read 3D polygon files in PLY format.
Copyright (c) 2004-2010 Oliver Kreylos

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

#include "ReadPlyFile.h"

#include <string.h>
#include <string>
#include <vector>
#include <iostream>
#include <Misc/SelfDestructPointer.h>
#include <Misc/ThrowStdErr.h>
#include <Misc/File.h>
#include <Misc/FileCharacterSource.h>
#include <Misc/ValueSource.h>

#include "LidarTypes.h"
#include "PointAccumulator.h"

namespace {

/**********************************
Enumerated type for PLY file modes:
**********************************/

enum PlyFileMode
	{
	PLY_WRONGFORMAT,PLY_ASCII,PLY_BINARY
	};

/*********************************************
Enumerated type for basic PLY file data types:
*********************************************/

enum DataType
	{
	CHAR,UCHAR,SHORT,USHORT,INT,UINT,FLOAT,DOUBLE
	};

/******************************************************
Templatized class to define PLY file value type traits:
******************************************************/

template <int dataTypeParam>
class DataValueTypes
	{
	};

template <>
class DataValueTypes<CHAR>
	{
	/* Embedded classes: */
	public:
	static const int dataType=CHAR;
	typedef char FileType;
	typedef int MemoryType;
	};

template <>
class DataValueTypes<UCHAR>
	{
	/* Embedded classes: */
	public:
	static const int dataType=UCHAR;
	typedef unsigned char FileType;
	typedef unsigned int MemoryType;
	};

template <>
class DataValueTypes<SHORT>
	{
	/* Embedded classes: */
	public:
	static const int dataType=SHORT;
	typedef short FileType;
	typedef int MemoryType;
	};

template <>
class DataValueTypes<USHORT>
	{
	/* Embedded classes: */
	public:
	static const int dataType=USHORT;
	typedef unsigned short FileType;
	typedef unsigned int MemoryType;
	};

template <>
class DataValueTypes<INT>
	{
	/* Embedded classes: */
	public:
	static const int dataType=INT;
	typedef int FileType;
	typedef int MemoryType;
	};

template <>
class DataValueTypes<UINT>
	{
	/* Embedded classes: */
	public:
	static const int dataType=UINT;
	typedef unsigned int FileType;
	typedef unsigned int MemoryType;
	};

template <>
class DataValueTypes<FLOAT>
	{
	/* Embedded classes: */
	public:
	static const int dataType=FLOAT;
	typedef float FileType;
	typedef double MemoryType;
	};

template <>
class DataValueTypes<DOUBLE>
	{
	/* Embedded classes: */
	public:
	static const int dataType=DOUBLE;
	typedef double FileType;
	typedef double MemoryType;
	};

/************************************************************
Templatized helper class to read values from ASCII PLY files:
************************************************************/

template <typename MemoryTypeParam>
class AsciiFileReader
	{
	};

template <>
class AsciiFileReader<int>
	{
	/* Methods: */
	public:
	static int readValue(Misc::ValueSource& asciiFile)
		{
		return asciiFile.readInteger();
		}
	};

template <>
class AsciiFileReader<unsigned int>
	{
	/* Methods: */
	public:
	static unsigned int readValue(Misc::ValueSource& asciiFile)
		{
		return asciiFile.readUnsignedInteger();
		}
	};

template <>
class AsciiFileReader<double>
	{
	/* Methods: */
	public:
	static double readValue(Misc::ValueSource& asciiFile)
		{
		return asciiFile.readNumber();
		}
	};

/**********************************************
Base class for data values read from PLY files:
**********************************************/

class DataValue
	{
	/* Constructors and destructors: */
	public:
	virtual ~DataValue(void)
		{
		}
	
	/* Methods: */
	virtual size_t getFileSize(void) const =0;
	virtual size_t getMemorySize(void) const =0;
	virtual void read(Misc::File& binaryPlyFile) =0;
	virtual void read(Misc::ValueSource& asciiPlyFile) =0;
	virtual int getInt(void) const =0;
	virtual unsigned int getUnsignedInt(void) const =0;
	virtual double getDouble(void) const =0;
	};

/*****************************************************
Templatized class for data values read from PLY files:
*****************************************************/

template <int dataTypeParam>
class DataValueTemplate:public DataValue,public DataValueTypes<dataTypeParam>
	{
	/* Embedded classes: */
	public:
	typedef DataValue Base1;
	typedef DataValueTypes<dataTypeParam> Base2;
	
	/* Elements: */
	private:
	typename Base2::MemoryType value; // Data value
	
	/* Methods: */
	virtual size_t getFileSize(void) const
		{
		return sizeof(typename Base2::FileType);
		}
	virtual size_t getMemorySize(void) const
		{
		return sizeof(typename Base2::MemoryType);
		}
	virtual void read(Misc::File& binaryPlyFile)
		{
		value=typename Base2::MemoryType(binaryPlyFile.read<typename Base2::FileType>());
		}
	virtual void read(Misc::ValueSource& asciiPlyFile)
		{
		value=AsciiFileReader<typename Base2::MemoryType>::readValue(asciiPlyFile);
		}
	virtual int getInt(void) const
		{
		return int(value);
		}
	virtual unsigned int getUnsignedInt(void) const
		{
		return (unsigned int)(value);
		}
	virtual double getDouble(void) const
		{
		return double(value);
		}
	};

/******************************************
Factory class to create data value readers:
******************************************/

class DataValueFactory
	{
	/* Methods: */
	public:
	static DataValue* newDataValue(DataType dataType)
		{
		switch(dataType)
			{
			case CHAR:
				return new DataValueTemplate<CHAR>;
			
			case UCHAR:
				return new DataValueTemplate<UCHAR>;
			
			case SHORT:
				return new DataValueTemplate<SHORT>;
			
			case USHORT:
				return new DataValueTemplate<USHORT>;
			
			case INT:
				return new DataValueTemplate<INT>;
			
			case UINT:
				return new DataValueTemplate<UINT>;
			
			case FLOAT:
				return new DataValueTemplate<FLOAT>;
			
			case DOUBLE:
				return new DataValueTemplate<DOUBLE>;
			}
		}
	};

/*****************************
Class for PLY file properties:
*****************************/

class Property
	{
	/* Embedded classes: */
	public:
	enum PropertyType
		{
		SCALAR,LIST
		};
	class Value
		{
		/* Elements: */
		private:
		const Property* property; // Pointer to property definition for this value
		DataValue* scalar; // Pointer to scalar value for scalar properties
		DataValue* listSize; // List size value for list properties
		std::vector<DataValue*> listElements; // Vector of pointers to list elements for list properties
		
		/* Private methods: */
		void clearListElements(void)
			{
			for(std::vector<DataValue*>::iterator eIt=listElements.begin();eIt!=listElements.end();++eIt)
				delete *eIt;
			listElements.clear();
			}
		
		/* Constructors and destructors: */
		public:
		Value(const Property* sProperty) // Creates empty value structure for the given property
			:property(sProperty),
			 scalar(0),
			 listSize(0)
			{
			if(property->getPropertyType()==Property::LIST)
				{
				/* Allocate space for list size: */
				listSize=DataValueFactory::newDataValue(property->getListSizeType());
				}
			else
				{
				/* Allocate space for scalar: */
				scalar=DataValueFactory::newDataValue(property->getScalarType());
				}
			}
		Value(const Value& source)
			:property(source.property),
			 scalar(0),
			 listSize(0)
			{
			if(property->getPropertyType()==Property::LIST)
				{
				/* Allocate space for list size: */
				listSize=DataValueFactory::newDataValue(property->getListSizeType());
				}
			else
				{
				/* Allocate space for scalar: */
				scalar=DataValueFactory::newDataValue(property->getScalarType());
				}
			}
		~Value(void)
			{
			if(property->getPropertyType()==Property::LIST)
				{
				delete listSize;
				clearListElements();
				}
			else
				delete scalar;
			}
		
		/* Methods: */
		template <class PlyFileParam>
		void read(PlyFileParam& plyFile) // Reads value from binary or ASCII PLY file
			{
			if(property->getPropertyType()==Property::LIST)
				{
				/* Read list size: */
				listSize->read(plyFile);
				unsigned int listSizeValue=listSize->getUnsignedInt();
				
				/* Read all list elements: */
				clearListElements();
				listElements.reserve(listSizeValue);
				for(unsigned int i=0;i<listSizeValue;++i)
					{
					listElements[i]=DataValueFactory::newDataValue(property->getListElementType());
					listElements[i]->read(plyFile);
					}
				}
			else
				{
				/* Read scalar: */
				scalar->read(plyFile);
				}
			}
		const DataValue* getScalar(void) const
			{
			return scalar;
			}
		const DataValue* getListSize(void) const
			{
			return listSize;
			}
		const DataValue* getListElement(unsigned int listElementIndex) const
			{
			return listElements[listElementIndex];
			}
		};
	
	/* Elements: */
	private:
	PropertyType propertyType; // Type of this property (scalar/list)
	DataType scalarType; // Data type for scalar properties
	DataType listSizeType; // Data type for list sizes for list properties
	DataType listElementType; // Data type for list elements for list properties
	std::string name; // Property name
	
	/* Private methods: */
	static DataType parseDataType(const std::string& tag)
		{
		static const char* dataTypeTags[]={"char","uchar","short","ushort","int","uint","float","double"};
		static const DataType dataTypes[]={CHAR,UCHAR,SHORT,USHORT,INT,UINT,FLOAT,DOUBLE};
		int i;
		for(i=0;i<8;++i)
			if(tag==dataTypeTags[i])
				break;
		if(i>=8)
			Misc::throwStdErr("Unknown data type %s",tag.c_str());
		return dataTypes[i];
		}
	
	/* Constructors and destructors: */
	public:
	Property(Misc::ValueSource& plyFile)
		{
		/* Read the property type: */
		std::string tag=plyFile.readString();
		
		if(tag=="list")
			{
			/* Parse a list property: */
			propertyType=LIST;
			listSizeType=parseDataType(plyFile.readString());
			listElementType=parseDataType(plyFile.readString());
			}
		else
			{
			/* Parse a scalar property: */
			propertyType=SCALAR;
			scalarType=parseDataType(tag);
			}
		
		/* Read the property name: */
		name=plyFile.readLine();
		plyFile.skipWs();
		}
	
	/* Methods: */
	PropertyType getPropertyType(void) const // Returns property's type
		{
		return propertyType;
		}
	DataType getScalarType(void) const // Returns scalar property's scalar type
		{
		return scalarType;
		}
	DataType getListSizeType(void) const // Returns list property's size type
		{
		return listSizeType;
		}
	DataType getListElementType(void) const // Returns list property's element type
		{
		return listElementType;
		}
	std::string getName(void) const // Returns property's name
		{
		return name;
		}
	};

/***************************
Class for PLY file elements:
***************************/

class Element
	{
	/* Embedded classes: */
	public:
	typedef std::vector<Property> PropertyList;
	class Value
		{
		/* Embedded classes: */
		public:
		typedef std::vector<Property::Value> PropertyValueList;
		
		/* Elements: */
		private:
		const Element* element; // Pointer to element definition for this value
		PropertyValueList propertyValues; // Vector of values for the properties of this element
		
		/* Constructors and destructors: */
		public:
		Value(const Element* sElement)
			:element(sElement)
			{
			/* Initialize vector of property values: */
			for(PropertyList::const_iterator plIt=element->propertiesBegin();plIt!=element->propertiesEnd();++plIt)
				{
				propertyValues.push_back(Property::Value(&(*plIt)));
				}
			}
		Value(const Value& source)
			:element(source.element)
			{
			/* Initialize vector of property values: */
			for(PropertyList::const_iterator plIt=element->propertiesBegin();plIt!=element->propertiesEnd();++plIt)
				{
				propertyValues.push_back(Property::Value(&(*plIt)));
				}
			}
		
		/* Methods: */
		template <class PlyFileParam>
		void read(PlyFileParam& plyFile) // Reads element from binary or ASCII PLY file
			{
			for(PropertyValueList::iterator pvIt=propertyValues.begin();pvIt!=propertyValues.end();++pvIt)
				pvIt->read(plyFile);
			}
		const Property::Value& getValue(unsigned int propertyIndex) const
			{
			return propertyValues[propertyIndex];
			}
		};
	
	/* Elements: */
	private:
	std::string name; // Name of this element
	PropertyList properties; // Vector of properties of this element
	
	/* Constructors and destructors: */
	public:
	Element(std::string sName)
		:name(sName)
		{
		}
	
	/* Methods: */
	void addProperty(Misc::ValueSource& plyFile)
		{
		properties.push_back(Property(plyFile));
		}
	size_t getNumProperties(void) const
		{
		return properties.size();
		}
	PropertyList::const_iterator propertiesBegin(void) const
		{
		return properties.begin();
		}
	PropertyList::const_iterator propertiesEnd(void) const
		{
		return properties.end();
		}
	unsigned int getPropertyIndex(std::string propertyName) const
		{
		unsigned int result=0;
		for(PropertyList::const_iterator pIt=properties.begin();pIt!=properties.end();++pIt,++result)
			if(pIt->getName()==propertyName)
				break;
		return result;
		}
	};

/*******************************************************
Helper function to determine the file type of PLY files:
*******************************************************/

struct PlyFileHeader // Structure containing relevant information from a PLY file's header
	{
	/* Elements: */
	public:
	bool isPlyFile; // Flag if the file is a PLY file
	PlyFileMode plyFileMode; // ASCII or binary
	Misc::File::Endianness plyFileEndianness; // Endianness of binary PLY file
	Element vertex; // Vertex element descriptor
	unsigned int numVertices; // Number of vertex elements
	Element face; // Face element descriptor
	unsigned int numFaces; // Number of face elements
	
	/* Constructors and destructors: */
	PlyFileHeader(void)
		:isPlyFile(false),plyFileMode(PLY_WRONGFORMAT),plyFileEndianness(Misc::File::DontCare),
		 vertex("vertex"),numVertices(0),
		 face("face"),numFaces(0)
		{
		}
	};

PlyFileHeader* readPlyFileHeader(const char* plyFileName)
	{
	/* Open the PLY file in text mode: */
	Misc::FileCharacterSource plyFile(plyFileName);
	Misc::ValueSource ply(plyFile);
	ply.skipWs();
	
	/* Process the PLY file header: */
	Misc::SelfDestructPointer<PlyFileHeader> result(new PlyFileHeader());
	int currentElement=-1;
	bool haveEndHeader=false;
	while(!ply.eof())
		{
		/* Read the next tag: */
		std::string tag=ply.readString();
		if(tag=="ply")
			result->isPlyFile=true;
		else if(tag=="format")
			{
			/* Read the format type and version number: */
			std::string format=ply.readString();
			if(format=="ascii")
				result->plyFileMode=PLY_ASCII;
			else if(format=="binary_little_endian")
				{
				result->plyFileMode=PLY_BINARY;
				result->plyFileEndianness=Misc::File::LittleEndian;
				}
			else if(format=="binary_big_endian")
				{
				result->plyFileMode=PLY_BINARY;
				result->plyFileEndianness=Misc::File::BigEndian;
				}
			double version=ply.readNumber();
			if(version!=1.0)
				result->isPlyFile=false;
			}
		else if(tag=="comment")
			{
			/* Skip the rest of the line: */
			ply.skipLine();
			ply.skipWs();
			}
		else if(tag=="element")
			{
			/* Read the element type: */
			std::string elementType=ply.readString();
			if(elementType=="vertex")
				{
				/* Parse a vertex element: */
				currentElement=0;
				result->numVertices=ply.readUnsignedInteger();
				}
			else if(elementType=="face")
				{
				/* Parse a face element: */
				currentElement=1;
				result->numFaces=ply.readUnsignedInteger();
				}
			else
				{
				/* Parse an unknown element: */
				currentElement=-1;
				ply.skipLine();
				ply.skipWs();
				}
			}
		else if(tag=="property")
			{
			/* Parse a property: */
			switch(currentElement)
				{
				case 0: // Vertex element
					result->vertex.addProperty(ply);
					break;
				
				case 1: // Face element
					result->face.addProperty(ply);
					break;
				
				default:
					; // Can't happen
				}
			}
		else if(tag=="end_header")
			{
			haveEndHeader=true;
			break;
			}
		else
			{
			/* Skip the unknown tag: */
			ply.skipLine();
			ply.skipWs();
			}
		}
	if(!haveEndHeader)
		result->isPlyFile=false;
	
	return result.releaseTarget();
	}

template <class PlyFileParam>
void readPlyFileVertices(PlyFileParam& ply,PlyFileHeader& header,PointAccumulator& pa,const float colorMask[3])
	{
	/* Get the indices of all relevant vertex value components: */
	Element& vertex=header.vertex;
	unsigned int posIndex[3];
	posIndex[0]=vertex.getPropertyIndex("x");
	posIndex[1]=vertex.getPropertyIndex("y");
	posIndex[2]=vertex.getPropertyIndex("z");
	unsigned int colIndex[3];
	colIndex[0]=vertex.getPropertyIndex("red");
	colIndex[1]=vertex.getPropertyIndex("green");
	colIndex[2]=vertex.getPropertyIndex("blue");
	bool hasColor=colIndex[0]<vertex.getNumProperties()&&colIndex[1]<vertex.getNumProperties()&&colIndex[2]<vertex.getNumProperties();
	
	/* Read all vertices: */
	Element::Value vertexValue(&vertex);
	for(unsigned int i=0;i<header.numVertices;++i)
		{
		/* Read vertex element from file: */
		vertexValue.read(ply);

		/* Extract vertex coordinates from vertex element: */
		LidarPoint p;
		for(int j=0;j<3;++j)
			p[j]=LidarPoint::Scalar(vertexValue.getValue(posIndex[j]).getScalar()->getDouble());
		
		/* Extract vertex color from vertex element: */
		if(hasColor)
			{
			for(int j=0;j<3;++j)
				{
				float col=float(vertexValue.getValue(colIndex[j]).getScalar()->getDouble())*colorMask[j];
				p.value[j]=Color::clampRound(col);
				}
			p.value[3]=Color::Scalar(255);
			}
		else
			{
			for(int j=0;j<3;++j)
				p.value[j]=Color::clampRound(255.0f*colorMask[j]);
			p.value[3]=Color::Scalar(255);
			}
		
		/* Store the point: */
		pa.addPoint(p);
		}
	}

}

void readPlyFile(PointAccumulator& pa,const char* fileName,const float colorMask[3])
	{
	/* Read the PLY file's header: */
	PlyFileHeader* header=readPlyFileHeader(fileName);
	if(!header->isPlyFile||header->plyFileMode==PLY_WRONGFORMAT)
		{
		std::cerr<<"Error: File "<<fileName<<" is not a valid PLY file"<<std::endl;
		delete header;
		return;
		}
	
	if(header->plyFileMode==PLY_BINARY)
		{
		/* Open the PLY file in binary mode: */
		Misc::File plyFile(fileName,"rb",header->plyFileEndianness);
		
		/* Skip the header: */
		while(true)
			{
			char line[256];
			plyFile.gets(line,sizeof(line));
			if(strcmp(line,"end_header\n")==0)
				break;
			}
		
		/* Read all vertices: */
		readPlyFileVertices(plyFile,*header,pa,colorMask);
		}
	else
		{
		/* Open the PLY file in text mode: */
		Misc::FileCharacterSource plyFile(fileName);
		Misc::ValueSource ply(plyFile);
		ply.skipWs();
		
		/* Skip the header: */
		while(ply.readString()!="end_header")
			;
		
		/* Read all vertices: */
		readPlyFileVertices(ply,*header,pa,colorMask);
		}
	
	/* Clean up and return: */
	delete header;
	}
