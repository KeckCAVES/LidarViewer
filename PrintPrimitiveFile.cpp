#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <Misc/ThrowStdErr.h>
#include <IO/File.h>
#include <IO/OpenFile.h>
#include <Geometry/Point.h>
#include <Geometry/Vector.h>
#include <Geometry/Endianness.h>
#include <Geometry/OutputOperators.h>

typedef Geometry::Point<double,3> Point;
typedef Geometry::Vector<double,3> Vector;

void processFile(const char* fileName,const std::vector<int> extractIndices,IO::File* outputFile)
	{
	/* Open the primitive file: */
	IO::AutoFile file(IO::openFile(fileName));
	file->setEndianness(IO::File::LittleEndian);
	
	/* Read the file header: */
	char header[40];
	file->read<char>(header,sizeof(header));
	if(strcmp(header,"LidarViewer primitive file v1.2       \n")!=0)
		{
		std::cerr<<"File "<<fileName<<" is not a valid version 1.2 primitive file"<<std::endl;
		return;
		}
	
	std::cout<<"Primitive file "<<fileName<<":"<<std::endl;
	
	/* Process all primitives in the file: */
	int primitiveIndex=0;
	while(!file->eof())
		{
		/* Check whether to save this primitive to the result file: */
		bool save=outputFile!=0;
		if(save)
			{
			save=false;
			for(std::vector<int>::const_iterator eaIt=extractIndices.begin();!save&&eaIt!=extractIndices.end();++eaIt)
				save=*eaIt==primitiveIndex;
			}
		
		/* Read the primitive type: */
		int primitiveType=file->read<int>();
		if(save)
			outputFile->write<int>(primitiveType);
		
		/* Read a primitive of the appropriate type: */
		switch(primitiveType)
			{
			case 0: // Plane primitive
				{
				size_t numPoints=file->read<unsigned int>();
				double rms=file->read<double>();
				std::cout<<"Plane primitive fitting "<<numPoints<<" points with RMS "<<rms<<":"<<std::endl;
				Vector normal=file->read<Vector>();
				double offset=file->read<double>();
				std::cout<<"  Plane equation: "<<normal<<" * X = "<<offset<<std::endl;
				std::cout<<"  Corner points:"<<std::endl;
				Point corners[4];
				for(int i=0;i<4;++i)
					{
					corners[i]=file->read<Point>();
					std::cout<<"    "<<corners[i]<<std::endl;
					}
				int numSegments[2];
				file->read(numSegments,2);
				std::cout<<"  Number of segments: "<<numSegments[0]<<" x "<<numSegments[1]<<std::endl;
				
				if(save)
					{
					outputFile->write<unsigned int>(numPoints);
					outputFile->write<double>(rms);
					outputFile->write<Vector>(normal);
					outputFile->write<double>(offset);
					outputFile->write<Point>(corners,4);
					outputFile->write<int>(numSegments,2);
					}
				
				break;
				}
			
			case 1: // Sphere primitive
				{
				size_t numPoints=file->read<unsigned int>();
				double rms=file->read<double>();
				std::cout<<"Sphere primitive fitting "<<numPoints<<" points with RMS "<<rms<<":"<<std::endl;
				Point center=file->read<Point>();
				double radius=file->read<double>();
				std::cout<<"  Center point: "<<center<<std::endl;
				std::cout<<"  Radius: "<<radius<<std::endl;
				
				if(save)
					{
					outputFile->write<unsigned int>(numPoints);
					outputFile->write<double>(rms);
					outputFile->write<Point>(center);
					outputFile->write<double>(radius);
					}
				
				break;
				}
			
			case 2: // Cylinder primitive
				{
				size_t numPoints=file->read<unsigned int>();
				double rms=file->read<double>();
				std::cout<<"Cylinder primitive fitting "<<numPoints<<" points with RMS "<<rms<<":"<<std::endl;
				Point center=file->read<Point>();
				Vector axis=file->read<Vector>();
				double radius=file->read<double>();
				double height=file->read<double>();
				int numSegments[2];
				file->read(numSegments,2);
				std::cout<<"  Center point: "<<center<<std::endl;
				std::cout<<"  Axis direction: "<<axis<<std::endl;
				std::cout<<"  Radius: "<<radius<<std::endl;
				std::cout<<"  Height: "<<height<<std::endl;
				std::cout<<"  Number of segments: "<<numSegments[0]<<" x "<<numSegments[1]<<std::endl;
				
				if(save)
					{
					outputFile->write<unsigned int>(numPoints);
					outputFile->write<double>(rms);
					outputFile->write<Point>(center);
					outputFile->write<Vector>(axis);
					outputFile->write<double>(radius);
					outputFile->write<double>(height);
					outputFile->write<int>(numSegments,2);
					}
				
				break;
				}
			
			case 3: // Line primitive
				{
				size_t numPoints=file->read<unsigned int>();
				double rms=file->read<double>();
				Point center=file->read<Point>();
				Vector axis=file->read<Vector>();
				double length=file->read<double>();
				std::cout<<"Line primitive fitting "<<numPoints<<" points with RMS "<<rms<<":"<<std::endl;
				std::cout<<"  Center point: "<<center<<std::endl;
				std::cout<<"  Axis direction: "<<axis<<std::endl;
				std::cout<<"  Length: "<<length<<std::endl;
				
				if(save)
					{
					outputFile->write<unsigned int>(numPoints);
					outputFile->write<double>(rms);
					outputFile->write<Point>(center);
					outputFile->write<Vector>(axis);
					outputFile->write<double>(length);
					}
				
				break;
				}
			
			case 4: // Point primitive
				{
				Point point=file->read<Point>();
				std::cout<<"Point primitive:"<<std::endl;
				std::cout<<"  Point: "<<point<<std::endl;
				
				if(save)
					{
					outputFile->write<Point>(point);
					}
				
				break;
				}
			
			default:
				std::cerr<<"Unknown primitive type "<<primitiveType<<" in primitive file "<<fileName<<std::endl;
				return;
			}
		std::cout<<std::endl;
		
		++primitiveIndex;
		}
	}

int main(int argc,char* argv[])
	{
	/* Process the command line: */
	IO::File* outputFile=0;
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"o")==0)
				{
				/* Close the previous output file: */
				delete outputFile;
				outputFile=0;
				
				/* Open a new output file: */
				++i;
				try
					{
					outputFile=IO::openFile(argv[i],IO::File::WriteOnly);
					outputFile->setEndianness(IO::File::LittleEndian);
					char header[40]="LidarViewer primitive file v1.2       \n";
					outputFile->write<char>(header,sizeof(header));
					}
				catch(std::runtime_error err)
					{
					std::cerr<<"Could not create output file "<<argv[i]<<" due to exception "<<err.what()<<std::endl;
					delete outputFile;
					outputFile=0;
					}
				}
			}
		else
			{
			/* Check if the user wants to extract primitives: */
			const char* fileName=argv[i];
			std::vector<int> extractIndices;
			while(i+1<argc&&isdigit(argv[i+1][0]))
				{
				++i;
				extractIndices.push_back(atoi(argv[i]));
				}
			
			/* Process the file: */
			processFile(fileName,extractIndices,outputFile);
			}
		}
	
	/* Close the last output file: */
	delete outputFile;
	
	return 0;
	}
