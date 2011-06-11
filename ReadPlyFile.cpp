/***********************************************************************
ReadPlyFile - Function to read 3D polygon files in PLY format.
Copyright (c) 2004-2011 Oliver Kreylos

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
#include <IO/File.h>
#include <IO/OpenFile.h>
#include <IO/ValueSource.h>

#include "LidarTypes.h"
#include "PointAccumulator.h"

namespace {

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
	static int readValue(IO::ValueSource& asciiFile)
		{
		return asciiFile.readInteger();
		}
	};

template <>
class AsciiFileReader<unsigned int>
	{
	/* Methods: */
	public:
	static unsigned int readValue(IO::ValueSource& asciiFile)
		{
		return asciiFile.readUnsignedInteger();
		}
	};

template <>
class AsciiFileReader<double>
	{
	/* Methods: */
	public:
	static double readValue(IO::ValueSource& asciiFile)
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
	virtual DataValue* clone(void) const =0;
	virtual size_t getFileSize(void) const =0;
	virtual size_t getMemorySize(void) const =0;
	virtual void read(IO::File& binaryPlyFile) =0;
	virtual void read(IO::ValueSource& asciiPlyFile) =0;
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
	virtual DataValue* clone(void) const
		{
		return new DataValueTemplate(*this);
		}
	virtual size_t getFileSize(void) const
		{
		return sizeof(typename Base2::FileType);
		}
	virtual size_t getMemorySize(void) const
		{
		return sizeof(typename Base2::MemoryType);
		}
	virtual void read(IO::File& binaryPlyFile)
		{
		value=typename Base2::MemoryType(binaryPlyFile.read<typename Base2::FileType>());
		}
	virtual void read(IO::ValueSource& asciiPlyFile)
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
		PropertyType propertyType; // Type of the value's property
		DataValue* scalar; // Pointer to scalar value for scalar properties
		DataValue* listSize; // List size value for list properties
		std::vector<DataValue*> listElements; // Vector of pointers to list elements for list properties
		
		/* Constructors and destructors: */
		public:
		Value(const Property& property) // Creates empty value structure for the given property
			:propertyType(property.getPropertyType()),
			 scalar(0),
			 listSize(0)
			{
			if(propertyType==Property::SCALAR)
				{
				/* Allocate space for scalar: */
				scalar=DataValueFactory::newDataValue(property.getScalarType());
				}
			else
				{
				/* Allocate space for list size: */
				listSize=DataValueFactory::newDataValue(property.getListSizeType());
				
				/* Create one list element to get started: */
				listElements.push_back(DataValueFactory::newDataValue(property.getListElementType()));
				}
			}
		private:
		Value(const Value& source); // Prohibit copy constructor
		Value& operator=(const Value& source); // Prohibit assignment operator
		public:
		~Value(void)
			{
			delete scalar;
			delete listSize;
			for(std::vector<DataValue*>::iterator eIt=listElements.begin();eIt!=listElements.end();++eIt)
				delete *eIt;
			}
		
		/* Methods: */
		size_t getFileSize(void) const // Returns value's size in binary files (returns minimal file size for list values)
			{
			if(propertyType==SCALAR)
				return scalar->getFileSize();
			else
				return listSize->getFileSize();
			}
		void skip(IO::ValueSource& plyFile) // Skips value from ASCII PLY file (never used)
			{
			if(propertyType==Property::SCALAR)
				{
				/* Skip scalar: */
				scalar->read(plyFile);
				}
			else
				{
				/* Read list size: */
				listSize->read(plyFile);
				unsigned int listSizeValue=listSize->getUnsignedInt();
				
				/* Skip all list elements: */
				for(unsigned int i=0;i<listSizeValue;++i)
					listElements[0]->read(plyFile);
				}
			}
		void skip(IO::File& plyFile) // Skips value from binary PLY file
			{
			if(propertyType==Property::SCALAR)
				{
				/* Skip scalar: */
				plyFile.skip<char>(scalar->getFileSize());
				}
			else
				{
				/* Read list size: */
				listSize->read(plyFile);
				unsigned int listSizeValue=listSize->getUnsignedInt();
				
				/* Skip all list elements: */
				plyFile.skip<char>(listElements[0]->getFileSize()*listSizeValue);
				}
			}
		template <class PlyFileParam>
		void read(PlyFileParam& plyFile) // Reads value from binary or ASCII PLY file
			{
			if(propertyType==Property::SCALAR)
				{
				/* Read scalar: */
				scalar->read(plyFile);
				}
			else
				{
				/* Read list size: */
				listSize->read(plyFile);
				unsigned int listSizeValue=listSize->getUnsignedInt();
				
				/* Ensure the list storage is long enough: */
				unsigned int currentListSize=listElements.size();
				while(currentListSize<listSizeValue)
					{
					listElements.push_back(listElements[0]->clone());
					++currentListSize;
					}
				
				/* Read all list elements: */
				for(unsigned int i=0;i<listSizeValue;++i)
					listElements[i]->read(plyFile);
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
	Property(IO::ValueSource& plyFile)
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
		typedef std::vector<Property::Value*> PropertyValueList;
		
		/* Elements: */
		private:
		const Element& element; // Pointer to element definition for this value
		PropertyValueList propertyValues; // Vector of values for the properties of this element
		
		/* Constructors and destructors: */
		public:
		Value(const Element& sElement)
			:element(sElement)
			{
			/* Initialize vector of property values: */
			for(PropertyList::const_iterator plIt=element.propertiesBegin();plIt!=element.propertiesEnd();++plIt)
				propertyValues.push_back(new Property::Value(*plIt));
			}
		private:
		Value(const Value& source); // Prohibit copy constructor
		Value& operator=(const Value& source); // Prohibit assignment operator
		public:
		~Value(void)
			{
			/* Destroy vector of property values: */
			for(PropertyValueList::iterator pvIt=propertyValues.begin();pvIt!=propertyValues.end();++pvIt)
				delete *pvIt;
			}
		
		/* Methods: */
		public:
		size_t getFileSize(void) const // Returns the total file size of all property values (returns minimal file size for any included list values)
			{
			size_t result=0;
			for(PropertyValueList::const_iterator pvIt=propertyValues.begin();pvIt!=propertyValues.end();++pvIt)
				result+=(*pvIt)->getFileSize();
			return result;
			}
		template <class PlyFileParam>
		void skip(PlyFileParam& plyFile) // Skips element value from binary or ASCII PLY file
			{
			for(PropertyValueList::iterator pvIt=propertyValues.begin();pvIt!=propertyValues.end();++pvIt)
				(*pvIt)->skip(plyFile);
			}
		template <class PlyFileParam>
		void read(PlyFileParam& plyFile) // Reads element value from binary or ASCII PLY file
			{
			for(PropertyValueList::iterator pvIt=propertyValues.begin();pvIt!=propertyValues.end();++pvIt)
				(*pvIt)->read(plyFile);
			}
		const Property::Value& getValue(unsigned int propertyIndex) const // Returns value of one of the element's properties
			{
			return *propertyValues[propertyIndex];
			}
		};
	
	/* Elements: */
	private:
	std::string name; // Name of this element
	size_t numValues; // Number of values of this element in the file
	PropertyList properties; // Vector of properties of this element
	
	/* Constructors and destructors: */
	public:
	Element(std::string sName,size_t sNumValues)
		:name(sName),numValues(sNumValues)
		{
		}
	
	/* Methods: */
	bool isElement(const char* elementName) const // Returns true if the element's name matches the given string
		{
		return name==elementName;
		}
	size_t getNumValues(void) const
		{
		return numValues;
		}
	void addProperty(IO::ValueSource& plyFile)
		{
		properties.push_back(Property(plyFile));
		}
	bool hasListProperty(void) const // Returns true if the element has at least one list property
		{
		bool result=false;
		for(PropertyList::const_iterator pIt=properties.begin();!result&&pIt!=properties.end();++pIt,++result)
			result=pIt->getPropertyType()==Property::LIST;
		return result;
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
	unsigned int getPropertyIndex(const char* propertyName) const
		{
		unsigned int result=0;
		for(PropertyList::const_iterator pIt=properties.begin();pIt!=properties.end();++pIt,++result)
			if(pIt->getName()==propertyName)
				break;
		return result;
		}
	};

/**************************************
Helper class to parse PLY file headers:
**************************************/

struct PlyFileHeader // Structure containing relevant information from a PLY file's header
	{
	/* Embedded classes: */
	public:
	enum FileType
		{
		Unknown,Ascii,Binary
		};
	
	/* Elements: */
	private:
	bool valid; // Flag if the file is a valid PLY file (as much as determined by parsing the header)
	FileType fileType; // ASCII or binary
	IO::File::Endianness fileEndianness; // Endianness of binary PLY file
	std::vector<Element> elements; // List of elements in the file, in the order in which they appear in the file
	
	/* Constructors and destructors: */
	public:
	PlyFileHeader(void)
		:valid(false),fileType(Unknown),fileEndianness(IO::File::DontCare)
		{
		}
	PlyFileHeader(IO::File& plyFile)
		:valid(false),fileType(Unknown),fileEndianness(IO::File::DontCare)
		{
		parseHeader(plyFile);
		}
	
	/* Methods: */
	bool parseHeader(IO::File& plyFile) // Creates header structure by parsing a PLY file's header
		{
		/* Attach a new value source to the PLY file: */
		IO::ValueSource ply(plyFile);
		ply.skipWs();
		
		/* Process the PLY file header: */
		std::vector<Element>::iterator currentElement=elements.end();
		bool isPly=false;
		bool haveEndHeader=false;
		while(!ply.eof())
			{
			/* Read the next tag: */
			std::string tag=ply.readString();
			if(tag=="ply")
				isPly=true;
			else if(tag=="format")
				{
				/* Read the format type and version number: */
				std::string format=ply.readString();
				if(format=="ascii")
					fileType=Ascii;
				else if(format=="binary_little_endian")
					{
					fileType=Binary;
					fileEndianness=IO::File::LittleEndian;
					}
				else if(format=="binary_big_endian")
					{
					fileType=Binary;
					fileEndianness=IO::File::BigEndian;
					}
				else
					{
					/* Unknown format; bail out: */
					break;
					}
				double version=ply.readNumber();
				if(version!=1.0)
					break;
				}
			else if(tag=="comment")
				{
				/* Skip the rest of the line: */
				ply.skipLine();
				ply.skipWs();
				}
			else if(tag=="element")
				{
				/* Read the element type and number of elements: */
				std::string elementType=ply.readString();
				size_t numElements=ply.readUnsignedInteger();
				
				/* Append a new element: */
				elements.push_back(Element(elementType,numElements));
				currentElement=elements.end()-1;
				}
			else if(tag=="property")
				{
				if(currentElement!=elements.end())
					{
					/* Parse a property: */
					currentElement->addProperty(ply);
					}
				else
					{
					/* Skip the property: */
					ply.skipLine();
					ply.skipWs();
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
		
		/* Check if the header was read completely: */
		valid=isPly&&haveEndHeader&&fileType!=Unknown;
		return valid;
		}
	bool isValid(void) const // Returns true if the header described a valid PLY file
		{
		return valid;
		}
	FileType getFileType(void) const // Returns the file's type
		{
		return fileType;
		}
	IO::File::Endianness getFileEndianness(void) const // Returns the endianness for binary PLY files
		{
		return fileEndianness;
		}
	size_t getNumElements(void) const // Returns the number of elements in the PLY file
		{
		return elements.size();
		}
	const Element& getElement(size_t index) const // Returns the element of the given index
		{
		return elements[index];
		}
	};

template <class PlyFileParam>
void readPlyFileVertices(const Element& vertex,PlyFileParam& ply,PointAccumulator& pa,const float colorMask[3])
	{
	/* Get the indices of all relevant vertex value components: */
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
	Element::Value vertexValue(vertex);
	for(unsigned int i=0;i<vertex.getNumValues();++i)
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

void skipElement(const Element& element,IO::File& ply)
	{
	/* Check if the element has variable size: */
	Element::Value value(element);
	if(element.hasListProperty())
		{
		/* Need to skip each value separately: */
		for(size_t i=0;i<element.getNumValues();++i)
			value.skip(ply);
		}
	else
		{
		/* Calculate the file size of each value of the element: */
		size_t valueSize=value.getFileSize();
		ply.skip<char>(valueSize*element.getNumValues());
		}
	}

void skipElement(const Element& element,IO::ValueSource& ply)
	{
	/* Skip one line for each value of the element: */
	for(size_t i=0;i<element.getNumValues();++i)
		ply.skipLine();
	ply.skipWs();
	}

template <class PlyFileParam>
void readPlyFileElements(const PlyFileHeader& header,PlyFileParam& ply,PointAccumulator& pa,const float colorMask[3])
	{
	/* Process all elements in order: */
	for(size_t elementIndex=0;elementIndex<header.getNumElements();++elementIndex)
		{
		/* Get the next element: */
		const Element& element=header.getElement(elementIndex);
		
		/* Check if it's the vertex element: */
		if(element.isElement("vertex"))
			{
			/* Read the vertex element: */
			readPlyFileVertices(element,ply,pa,colorMask);
			}
		else
			{
			/* Skip the entire element: */
			skipElement(element,ply);
			}
		}
	}

}

void readPlyFile(PointAccumulator& pa,const char* fileName,const float colorMask[3])
	{
	/* Open the PLY file: */
	IO::AutoFile plyFile(IO::openFile(fileName));
	
	/* Read the PLY file's header: */
	PlyFileHeader header(*plyFile);
	if(!header.isValid())
		{
		std::cerr<<"Error: File "<<fileName<<" is not a valid PLY file"<<std::endl;
		return;
		}
	
	/* Read the PLY file in ASCII or binary mode: */
	if(header.getFileType()==PlyFileHeader::Ascii)
		{
		/* Attach a value source to the PLY file: */
		IO::ValueSource ply(*plyFile);
		
		/* Read the PLY file in ASCII mode: */
		readPlyFileElements(header,ply,pa,colorMask);
		}
	else if(header.getFileType()==PlyFileHeader::Binary)
		{
		/* Set the PLY file's endianness: */
		plyFile->setEndianness(header.getFileEndianness());
		
		/* Read the PLY file in binary mode: */
		readPlyFileElements(header,*plyFile,pa,colorMask);
		}
	}
