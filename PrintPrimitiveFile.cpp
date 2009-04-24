#include <string.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <Misc/ThrowStdErr.h>
#include <Misc/File.h>
#include <Geometry/Point.h>
#include <Geometry/Vector.h>
#include <Geometry/Endianness.h>

typedef Geometry::Point<double,3> Point;
typedef Geometry::Vector<double,3> Vector;

int main(int argc,char* argv[])
	{
	/* Open the primitive file: */
	Misc::File file(argv[1],"rb",Misc::File::LittleEndian);
	
	/* Read the file header: */
	char header[40];
	file.read<char>(header,sizeof(header));
	if(strcmp(header,"LidarViewer primitive file v1.2       \n")!=0)
		{
		std::cerr<<"File "<<argv[1]<<" is not a valid version 1.2 primitive file"<<std::endl;
		return 1;
		}
	
	/* Read all primitives stored in the file: */
	std::cout.precision(10);
	while(true)
		{
		/* Read the primitive type: */
		int primitiveType;
		try
			{
			primitiveType=file.read<int>();
			}
		catch(Misc::File::ReadError err)
			{
			/* Stop reading: */
			break;
			}
		
		/* Read a primitive of the appropriate type: */
		switch(primitiveType)
			{
			case 0: // Plane primitive
				{
				size_t numPoints=file.read<unsigned int>();
				double rms=file.read<double>();
				std::cout<<"Plane primitive fitting "<<numPoints<<" points with RMS "<<rms<<":"<<std::endl;
				Vector normal=file.read<Vector>();
				double offset=file.read<double>();
				std::cout<<"  Plane equation: ("<<normal[0]<<", "<<normal[1]<<", "<<normal[2]<<") * X = "<<offset<<std::endl;
				std::cout<<"  Corner points:"<<std::endl;
				for(int i=0;i<4;++i)
					{
					Point corner=file.read<Point>();
					std::cout<<"    ("<<corner[0]<<", "<<corner[1]<<", "<<corner[2]<<")"<<std::endl;
					}
				int numSegments[2];
				file.read(numSegments,2);
				std::cout<<"  Number of segments: "<<numSegments[0]<<" x "<<numSegments[1]<<std::endl;
				break;
				}
			
			case 1: // Sphere primitive
				{
				size_t numPoints=file.read<unsigned int>();
				double rms=file.read<double>();
				std::cout<<"Sphere primitive fitting "<<numPoints<<" points with RMS "<<rms<<":"<<std::endl;
				Point center=file.read<Point>();
				double radius=file.read<double>();
				std::cout<<"  Center point: ("<<center[0]<<", "<<center[1]<<", "<<center[2]<<")"<<std::endl;
				std::cout<<"  Radius: "<<radius<<std::endl;
				break;
				}
			
			case 2: // Cylinder primitive
				{
				size_t numPoints=file.read<unsigned int>();
				double rms=file.read<double>();
				std::cout<<"Cylinder primitive fitting "<<numPoints<<" points with RMS "<<rms<<":"<<std::endl;
				Point center=file.read<Point>();
				Vector axis=file.read<Vector>();
				double radius=file.read<double>();
				double height=file.read<double>();
				int numSegments[2];
				file.read(numSegments,2);
				std::cout<<"  Center point: ("<<center[0]<<", "<<center[1]<<", "<<center[2]<<")"<<std::endl;
				std::cout<<"  Axis direction: ("<<axis[0]<<", "<<axis[1]<<", "<<axis[2]<<")"<<std::endl;
				std::cout<<"  Radius: "<<radius<<std::endl;
				std::cout<<"  Height: "<<height<<std::endl;
				std::cout<<"  Number of segments: "<<numSegments[0]<<" x "<<numSegments[1]<<std::endl;
				break;
				}
			
			case 3: // Line primitive
				{
				Point center=file.read<Point>();
				Vector axis=file.read<Vector>();
				double length=file.read<double>();
				std::cout<<"Line primitive:"<<std::endl;
				std::cout<<"  Center point: ("<<center[0]<<", "<<center[1]<<", "<<center[2]<<")"<<std::endl;
				std::cout<<"  Axis direction: ("<<axis[0]<<", "<<axis[1]<<", "<<axis[2]<<")"<<std::endl;
				std::cout<<"  Length: "<<length<<std::endl;
				break;
				}
			
			case 4: // Point primitive
				{
				Point point=file.read<Point>();
				std::cout<<"Point primitive:"<<std::endl;
				std::cout<<"  Point: ("<<point[0]<<", "<<point[1]<<", "<<point[2]<<")"<<std::endl;
				break;
				}
			
			default:
				std::cerr<<"Unknown primitive type "<<primitiveType<<" in primitive file "<<argv[1]<<std::endl;
				return 1;
			}
		std::cout<<std::endl;
		}
	
	return 0;
	}

