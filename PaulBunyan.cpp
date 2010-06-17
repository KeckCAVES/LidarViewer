#include <string.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <Threads/Thread.h>
#include <Threads/Barrier.h>
#include <Math/Math.h>

#include "LidarTypes.h"
#include "LidarProcessOctree.h"

class Axe
	{
	/* Elements: */
	private:
	Scalar radius,radius2;
	Scalar tolerance;
	Scalar maxAngle;
	
	Point candidate;
	bool haveD0;
	Vector d0;
	Vector dl,dr;
	
	/* Constructors and destructors: */
	public:
	Axe(Scalar sRadius,Scalar sTolerance,Scalar sMaxAngle)
		:radius(sRadius),radius2(Math::sqr(radius)),tolerance(sTolerance),maxAngle(sMaxAngle)
		{
		}
	
	/* Methods: */
	void prepareTraversal(const Point& sCandidate)
		{
		candidate=sCandidate;
		haveD0=false;
		}
	Box getBox(void) const
		{
		Box result;
		for(int i=0;i<2;++i)
			{
			result.min[i]=candidate[i]-radius;
			result.max[i]=candidate[i]+radius;
			}
		result.min[2]=Box::full.min[2];
		result.max[2]=candidate[2]-tolerance;
		return result;
		}
	void operator()(const LidarPoint& lp)
		{
		/* Check if the point is inside the search cylinder: */
		Scalar r2=Math::sqr(lp[0]-candidate[0])+Math::sqr(lp[1]-candidate[1]);
		Scalar cyl=candidate[2]-tolerance-radius;
		if(r2<=radius2&&(lp[2]<=cyl||r2+Math::sqr(lp[2]-cyl)<=radius2))
			{
			/* Check if this is the first encountered point: */
			Vector d=lp-candidate;
			d[2]=Scalar(0);
			if(haveD0)
				{
				#if 1
				/* Check if the point is on the left or right side of d0: */
				if(d*d0>=Scalar(0))
					{
					/* Check the left vector: */
					if(d*dl>Scalar(0))
						dl=Vector(-d[1],d[0],0);
					}
				else
					{
					/* Check the right vector: */
					if(d*dr<Scalar(0))
						dr=Vector(-d[1],d[0],0);
					}
				#endif
				}
			else
				{
				/* Initialize the area finder: */
				d0=Vector(-d[1],d[0],0);
				d0.normalize();
				dl=dr=d0;
				haveD0=true;
				}
			}
		}
	bool isGround(void) const
		{
		if(haveD0)
			{
			#if 1
			/* Check that the angle formed by the left and right vectors is less than maxAngle: */
			double cosAlpha1=(dl*d0)/dl.mag();
			if(cosAlpha1>Scalar(1))
				cosAlpha1=Scalar(1);
			if(cosAlpha1<Scalar(-1))
				cosAlpha1=Scalar(-1);
			double cosAlpha2=(dr*d0)/dr.mag();
			if(cosAlpha2>Scalar(1))
				cosAlpha2=Scalar(1);
			if(cosAlpha2<Scalar(-1))
				cosAlpha2=Scalar(-1);
			return Math::acos(cosAlpha1)+Math::acos(cosAlpha2)<maxAngle;
			#else
			return false;
			#endif
			}
		else
			return true;
		}
	};

class PaulBunyan
	{
	/* Elements: */
	private:
	LidarProcessOctree& lpo; // The processed LiDAR octree
	Axe axe; // The axe
	Color* colorBuffer; // Array to hold colors for a node during processing
	Color* childColorBuffers[8]; // Arrays to hold colors for a node's children during processing
	LidarFile::Offset colorDataSize; // Size of each record in the color file
	LidarFile colorFile; // The file to which to write the color data
	FILE* pointFile; // Optional ASCII file to store ground points
	unsigned int numThreads; // Number of processing threads
	bool shutdownThreads; // Flag to shut down threads at the end
	Threads::Thread* calcThreads; // Array of threads to calculate normal vectors from point neighborhoods
	Threads::Barrier calcBarrier;
	Threads::Thread* subsampleThreads; // Array of threads to subsample normal vectors from node children
	Threads::Barrier subsampleBarrier;
	LidarProcessOctree::Node* currentNode; // Pointer to currently processed node
	LidarProcessOctree::Node* currentChildren[8]; // Pointer to children of currently processed node
	size_t numProcessedNodes; // Number of already processed nodes
	
	/* Private methods: */
	void* calcThreadMethod(unsigned int threadIndex);
	void* subsampleThreadMethod(unsigned int threadIndex);
	
	/* Constructors and destructors: */
	public:
	PaulBunyan(LidarProcessOctree& sLpo,const Axe& sAxe,const char* colorFileName,const char* pointFileName,unsigned int sNumThreads =1); // Creates an axe wielder with the given axe
	~PaulBunyan(void);
	
	/* Methods: */
	void operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel);
	};

/***************************
Methods of class PaulBunyan:
***************************/

void* PaulBunyan::calcThreadMethod(unsigned int threadIndex)
	{
	Threads::Thread::setCancelState(Threads::Thread::CANCEL_ENABLE);
	
	while(true)
		{
		/* Wait on the calculation barrier until there is a job: */
		calcBarrier.synchronize();
		if(shutdownThreads)
			break;
		
		/* Process the current node: */
		unsigned int firstI=(threadIndex*currentNode->getNumPoints())/numThreads;
		unsigned int lastI=((threadIndex+1)*currentNode->getNumPoints())/numThreads;
		for(unsigned int index=firstI;index<lastI;++index)
			{
			/* Apply the axe to the point: */
			axe.prepareTraversal((*currentNode)[index]);
			lpo.processPointsInBox(axe.getBox(),axe);
			
			/* Color the point based on its classification: */
			float intensity=float((*currentNode)[index].value[0])*0.299f+float((*currentNode)[index].value[1])*0.587f+float((*currentNode)[index].value[2])*0.114f;
			for(int i=0;i<4;++i)
				colorBuffer[index][i]=Color::Scalar(0);
			if(axe.isGround())
				{
				/* Color the point brown: */
				colorBuffer[index][0]=Color::Scalar(intensity*0.5f+0.5f);
				colorBuffer[index][1]=Color::Scalar(intensity*0.333f+0.5f);
				
				if(pointFile!=0)
					{
					/* Write the original point to a file: */
					fprintf(pointFile,"%.10g %.10g %.10g %3u %3u %3u\n",
					        (*currentNode)[index][0],(*currentNode)[index][1],(*currentNode)[index][2],
					        (*currentNode)[index].value[0],(*currentNode)[index].value[1],(*currentNode)[index].value[2]);
					}
				}
			else
				{
				/* Color the point green: */
				colorBuffer[index][1]=Color::Scalar(intensity+0.5f);
				}
			}
		
		/* Synchronize on the calculation barrier to signal job completion: */
		calcBarrier.synchronize();
		}
	
	return 0;
	}

namespace {

/**************
Helper classes:
**************/

class FindPoint // Class to find a point inside an octree node
	{
	/* Elements: */
	private:
	Point queryPoint; // The position of the point to find
	const LidarPoint* foundPoint; // The found LiDAR point
	
	/* Constructors and destructors: */
	public:
	FindPoint(const Point& sQueryPoint)
		:queryPoint(sQueryPoint),
		 foundPoint(0)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& lp)
		{
		if(Geometry::sqrDist(lp,queryPoint)==Scalar(0))
			foundPoint=&lp;
		}
	const Point& getQueryPoint(void) const
		{
		return queryPoint;
		}
	Scalar getQueryRadius2(void) const
		{
		return Scalar(0);
		}
	const LidarPoint* getFoundPoint(void) const
		{
		return foundPoint;
		}
	};

}

void* PaulBunyan::subsampleThreadMethod(unsigned int threadIndex)
	{
	Threads::Thread::setCancelState(Threads::Thread::CANCEL_ENABLE);
	
	while(true)
		{
		/* Wait on the subsampling barrier until there is a job: */
		subsampleBarrier.synchronize();
		if(shutdownThreads)
			break;
		
		/* Find the ancestor points of each node point and copy their classified colors: */
		unsigned int firstI=(threadIndex*currentNode->getNumPoints())/numThreads;
		unsigned int lastI=((threadIndex+1)*currentNode->getNumPoints())/numThreads;
		for(unsigned int i=firstI;i<lastI;++i)
			{
			/* Find the child node containing this point's ancestor: */
			int pointChildIndex=currentNode->getDomain().findChild((*currentNode)[i]);
			LidarProcessOctree::Node* pointChild=currentChildren[pointChildIndex];
			
			/* Find the point's ancestor: */
			FindPoint fp((*currentNode)[i]);
			lpo.processNodePointsDirected(pointChild,fp);
			if(fp.getFoundPoint()==0)
				Misc::throwStdErr("Things are fucked up!");
			
			/* Copy the ancestor's color: */
			colorBuffer[i]=childColorBuffers[pointChildIndex][fp.getFoundPoint()-pointChild->getPoints()];
			}
		
		/* Synchronize on the subsampling barrier to signal job completion: */
		subsampleBarrier.synchronize();
		}
	
	return 0;
	}

PaulBunyan::PaulBunyan(LidarProcessOctree& sLpo,const Axe& sAxe,const char* colorFileName,const char* pointFileName,unsigned int sNumThreads)
	:lpo(sLpo),axe(sAxe),
	 colorBuffer(new Color[lpo.getMaxNumPointsPerNode()]),
	 colorDataSize(sizeof(Color)),
	 colorFile(colorFileName,"w+b",LidarFile::LittleEndian),
	 pointFile(0),
	 numThreads(sNumThreads),shutdownThreads(false),
	 calcThreads(new Threads::Thread[numThreads]),calcBarrier(numThreads+1),
	 subsampleThreads(new Threads::Thread[numThreads]),subsampleBarrier(numThreads+1),
	 currentNode(0),
	 numProcessedNodes(0)
	{
	/* Allocate the child color buffers: */
	for(int i=0;i<8;++i)
		childColorBuffers[i]=new Color[lpo.getMaxNumPointsPerNode()];
	
	/* Write the color file's header: */
	LidarDataFileHeader dfh((unsigned int)(colorDataSize));
	dfh.write(colorFile);
	
	if(pointFileName!=0)
		{
		/* Open the point file: */
		pointFile=fopen(pointFileName,"wt");
		}
	
	/* Start the worker threads: */
	for(unsigned int i=0;i<numThreads;++i)
		{
		calcThreads[i].start(this,&PaulBunyan::calcThreadMethod,i);
		subsampleThreads[i].start(this,&PaulBunyan::subsampleThreadMethod,i);
		}
	}

PaulBunyan::~PaulBunyan(void)
	{
	/* Shut down all threads: */
	shutdownThreads=true;
	calcBarrier.synchronize();
	subsampleBarrier.synchronize();
	for(unsigned int i=0;i<numThreads;++i)
		{
		calcThreads[i].join();
		subsampleThreads[i].join();
		}
	delete[] calcThreads;
	delete[] subsampleThreads;
	
	delete[] colorBuffer;
	for(int i=0;i<8;++i)
		delete[] childColorBuffers[i];
	if(pointFile!=0)
		fclose(pointFile);
	}

void PaulBunyan::operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel)
	{
	currentNode=&node;
	
	if(currentNode->getNumPoints()>0)
		{
		if(currentNode->isLeaf())
			{
			/* Wake up the calculation threads: */
			calcBarrier.synchronize();
			
			/* Wait for their completion: */
			calcBarrier.synchronize();
			}
		else
			{
			/* Get pointers to the node's children and load their classified color arrays: */
			for(int childIndex=0;childIndex<8;++childIndex)
				{
				currentChildren[childIndex]=lpo.getChild(currentNode,childIndex);
				if(currentChildren[childIndex]->getNumPoints()>0)
					{
					colorFile.seekSet(LidarDataFileHeader::getFileSize()+colorDataSize*currentChildren[childIndex]->getDataOffset());
					colorFile.read(childColorBuffers[childIndex],currentChildren[childIndex]->getNumPoints());
					}
				}
			
			/* Wake up the subsampling threads: */
			subsampleBarrier.synchronize();
			
			/* Wait for their completion: */
			subsampleBarrier.synchronize();
			}
		}
	
	/* Write the node's colors to the color file: */
	colorFile.seekSet(LidarDataFileHeader::getFileSize()+colorDataSize*currentNode->getDataOffset());
	colorFile.write(colorBuffer,currentNode->getNumPoints());
	
	/* Update the progress counter: */
	++numProcessedNodes;
	int percent=int((numProcessedNodes*100+lpo.getNumNodes()/2)/lpo.getNumNodes());
	std::cout<<"\b\b\b\b"<<std::setw(3)<<percent<<"%"<<std::flush;
	}

int main(int argc,char* argv[])
	{
	const char* lidarFileName=0;
	int cacheSize=512;
	unsigned int numThreads=1;
	const char* colorFileName=0;
	const char* groundPointFileName=0;
	Scalar radius(1);
	Scalar tolerance(0);
	Scalar maxAngle(270);
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"cache")==0)
				{
				++i;
				cacheSize=atoi(argv[i]);
				}
			else if(strcasecmp(argv[i]+1,"threads")==0)
				{
				++i;
				numThreads=atoi(argv[i]);
				}
			else if(strcasecmp(argv[i]+1,"radius")==0)
				{
				++i;
				radius=Scalar(atof(argv[i]));
				}
			else if(strcasecmp(argv[i]+1,"tolerance")==0)
				{
				++i;
				tolerance=Scalar(atof(argv[i]));
				}
			else if(strcasecmp(argv[i]+1,"maxAngle")==0)
				{
				++i;
				maxAngle=Math::rad(Scalar(atof(argv[i])));
				}
			}
		else if(lidarFileName==0)
			lidarFileName=argv[i];
		else if(colorFileName==0)
			colorFileName=argv[i];
		else if(groundPointFileName==0)
			groundPointFileName=argv[i];
		}
	if(lidarFileName==0)
		{
		std::cerr<<"No LiDAR file name provided"<<std::endl;
		return 1;
		}
	if(colorFileName==0)
		colorFileName="Colors";
	
	/* Create a processing octree: */
	LidarProcessOctree lpo(lidarFileName,size_t(cacheSize)*size_t(1024*1024));
	
	/* Create an axe: */
	Axe axe(radius,tolerance,maxAngle);
	
	/* Hand the axe to Paul Bunyan: */
	std::string lidarColorFileName=lidarFileName;
	lidarColorFileName.push_back('/');
	lidarColorFileName.append(colorFileName);
	PaulBunyan paulBunyan(lpo,axe,lidarColorFileName.c_str(),groundPointFileName,numThreads);
	std::cout<<"Wielding axe...   0%"<<std::flush;
	lpo.processNodesPostfix(paulBunyan);
	std::cout<<std::endl;
	
	return 0;
	}
