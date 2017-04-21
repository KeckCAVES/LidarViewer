/***********************************************************************
LidarColorMapper - Post-processing filter to assign image colors to each
point in a LiDAR data set.
Copyright (c) 2009-2013 Oliver Kreylos

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

#include <utility>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <Misc/ThrowStdErr.h>
#include <Misc/HashTable.h>
#include <Misc/FileNameExtensions.h>
#include <IO/File.h>
#include <IO/OpenFile.h>
#include <IO/ValueSource.h>
#include <Geometry/Point.h>
#include <Geometry/Box.h>
#include <Geometry/AffineTransformation.h>
#include <Geometry/ProjectiveTransformation.h>
#include <Geometry/LambertConformalProjection.h>
#include <Geometry/UTMProjection.h>
#include <Geometry/OutputOperators.h>
#include <Images/RGBImage.h>
#include <Images/GetImageFileSize.h>
#include <Images/ReadImageFile.h>

#include "LidarTypes.h"
#include "LidarProcessOctree.h"

/*******************************************************************
Class to represent color images with 2D projection into world space:
*******************************************************************/

class Image2D
	{
	/* Embedded classes: */
	public:
	typedef Geometry::Point<double,2> Point; // Type for world and image points
	typedef Geometry::Box<double,2> Box; // Type for boxes in world and image space
	typedef Geometry::AffineTransformation<double,2> Transform; // Type for image transformations
	typedef Images::RGBImage::Color Color; // Type for image colors
	typedef std::pair<bool,Color> SampleResult; // Type for results of sampling an image
	
	/* Elements: */
	private:
	std::string fileName; // The image's file name
	unsigned int size[2]; // The image's width and height in pixels
	Transform transform; // The transformation from image coordinates to world coordinates
	Transform invTransform; // The transformation from world coordinates to image coordinates
	Box imageBox; // Bounding box of image in image coordinates
	Box worldBox; // Bounding box of image in world coordinates
	Images::RGBImage image; // Image data; invalid if image is not paged in
	
	/* Constructors and destructors: */
	public:
	Image2D(const char* sFileName); // Loads an image and world space information from a set of files
	
	/* Methods: */
	const unsigned int* getSize(void) const // Returns the image's size in pixels
		{
		return size;
		}
	const Box& getWorldBox(void) const // Returns the image's bounding box in world space
		{
		return worldBox;
		}
	void pageIn(void); // Ensures that the image's data resides in main memory
	void pageOut(void); // Releases the image's data from main memory
	SampleResult getColor(const Point& worldPos) const; // Returns the image's color for the given point in world coordinates
	};

/************************
Methods of class Image2D:
************************/

Image2D::Image2D(const char* sFileName)
	:fileName(sFileName)
	{
	/* Construct the file name of the world space metafile: */
	const char* extPtr=Misc::getExtension(fileName.c_str());
	std::string worldFileName=std::string(fileName.c_str(),extPtr)+".tfw";
	
	/* Read the world space metafile: */
	IO::ValueSource world(IO::openFile(worldFileName.c_str()));
	world.skipWs();
	for(int j=0;j<3;++j)
		for(int i=0;i<2;++i)
			transform.getMatrix()(i,j)=world.readNumber();
	
	/* Determine the image size: */
	Images::getImageFileSize(fileName.c_str(),size[0],size[1]);
	
	/* Flip the transformation: */
	for(int i=0;i<2;++i)
		{
		transform.getMatrix()(i,1)=-transform.getMatrix()(i,1);
		transform.getMatrix()(i,2)-=transform.getMatrix()(i,1)*double(size[1]-1);
		}
	
	/* Calculate the inverse transformation: */
	invTransform=Geometry::invert(transform);
	
	/* Calculate the bounding box in image space: */
	imageBox.min=Point::origin;
	imageBox.max=Point(double(size[0]-1),double(size[1]-1));
	
	/* Calculate the bounding box in world space: */
	worldBox=imageBox;
	worldBox.transform(transform);
	
	/* Print the world space bounding box: */
	std::cout<<"World space bounding box: "<<worldBox.min<<", "<<worldBox.max<<std::endl;
	}

void Image2D::pageIn(void)
	{
	if(!image.isValid())
		{
		/* Load the image file: */
		std::cout<<"\rPaging in image file "<<fileName<<"..."<<std::flush;
		image=Images::readImageFile(fileName.c_str());
		std::cout<<" done"<<std::endl;
		}
	}

void Image2D::pageOut(void)
	{
	/* Replace the image with an invalid image: */
	std::cout<<"\rPaging out image file "<<fileName<<std::endl;
	image=Images::RGBImage();
	}

Image2D::SampleResult Image2D::getColor(const Image2D::Point& worldPos) const
	{
	/* Transform the world position to image space: */
	Point imagePos=invTransform.transform(worldPos);
	
	/* Bail out if the image position is outside the image: */
	if(imagePos[0]<imageBox.min[0]||imagePos[0]>=imageBox.max[0]||imagePos[1]<imageBox.min[1]||imagePos[1]>=imageBox.max[1])
		return SampleResult(false,Color(0,0,0));
	
	/* Sample the image: */
	unsigned int cx=(unsigned int)imagePos[0];
	double dx=imagePos[0]-double(cx);
	unsigned int cy=(unsigned int)imagePos[1];
	double dy=imagePos[1]-double(cy);
	const Color* r0=image.getPixelRow(cy);
	const Color* r1=image.getPixelRow(cy+1);
	double p0[3],p1[3];
	for(int i=0;i<3;++i)
		{
		p0[i]=double(r0[cx][i])*(1.0-dx)+double(r0[cx+1][i])*dx;
		p1[i]=double(r1[cx][i])*(1.0-dx)+double(r1[cx+1][i])*dx;
		}
	Color result;
	for(int i=0;i<3;++i)
		{
		double v=p0[i]*(1.0-dy)+p1[i]*dy;
		if(v<0.5)
			result[i]=Color::Scalar(0);
		else if(v>=254.5)
			result[i]=Color::Scalar(255);
		else
			result[i]=Color::Scalar(v+0.5);
		}
	return SampleResult(true,result);
	}

/************************************************
Helper class to manage an LRU cache of 2D images:
************************************************/

class ImageCacher
	{
	/* Embedded classes: */
	private:
	struct LRUItem // Structure for least-recently-used list items
		{
		/* Elements: */
		Image2D* image; // Pointer to the image
		size_t imageSize; // Image memory footprint in bytes
		bool pagedIn; // Flag if the image's data are currently residing in main memory
		unsigned int pageInCounter; // Request counter value at which this image was last paged in
		LRUItem* pred; // Pointer to preceding entry in LRU list
		LRUItem* succ; // Pointer to succeeding entry in LRU list
		};
	
	/* Elements: */
	private:
	size_t maxMemory; // Allocated amount of memory in bytes
	Misc::HashTable<Image2D*,LRUItem*> lruMap; // Map from image pointers to least-recently-used list items
	LRUItem* lruHead; // Pointer to the least-recently-used image
	LRUItem* lruTail; // Pointer to the most-recently-used image
	size_t usedMemory; // Currently used amount of memory in bytes
	unsigned int requestCounter; // Counter to keep track of cache memory overflow
	unsigned int numPageInRequests; // Total number of times an image was loaded into memory
	
	/* Constructors and destructors: */
	public:
	ImageCacher(size_t sMaxMemory); // Creates an image cacher with the given memory size in bytes
	~ImageCacher(void);
	
	/* Methods: */
	void registerImage(Image2D* image); // Registers an image with the cache manager
	void requestImages(const std::vector<Image2D*>& images); // Requests in-memory access to the given set of images
	unsigned int getNumPageInRequests(void) const
		{
		return numPageInRequests;
		}
	};

/****************************
Methods of class ImageCacher:
****************************/

ImageCacher::ImageCacher(size_t sMaxMemory)
	:maxMemory(sMaxMemory),
	 lruMap(101),
	 lruHead(0),lruTail(0),
	 usedMemory(0),
	 requestCounter(0),
	 numPageInRequests(0)
	{
	}

ImageCacher::~ImageCacher(void)
	{
	/* Destroy all LRU list items: */
	for(Misc::HashTable<Image2D*,LRUItem*>::Iterator lmIt=lruMap.begin();!lmIt.isFinished();++lmIt)
		delete lmIt->getDest();
	}

void ImageCacher::registerImage(Image2D* image)
	{
	/* Create a new LRU list item: */
	LRUItem* newItem=new LRUItem;
	newItem->image=image;
	newItem->imageSize=size_t(image->getSize()[1])*size_t(image->getSize()[0])*sizeof(Image2D::Color);
	newItem->pagedIn=false;
	newItem->pred=0;
	newItem->succ=0;
	
	/* Store the LRU list item in the item map: */
	lruMap[image]=newItem;
	}

void ImageCacher::requestImages(const std::vector<Image2D*>& images)
	{
	/* Add all images in the list to the LRU cache in turn: */
	for(std::vector<Image2D*>::const_iterator iIt=images.begin();iIt!=images.end();++iIt)
		{
		/* Check if the image is already paged in: */
		LRUItem* item=lruMap[*iIt].getDest();
		if(item->pagedIn)
			{
			/* Move the image's LRU list item to the end of the LRU list: */
			if(lruTail!=item)
				{
				/* Unlink the item from its current place in the list: */
				if(item->pred!=0)
					item->pred->succ=item->succ;
				else
					lruHead=item->succ;
				if(item->succ!=0)
					item->succ->pred=item->pred;
				
				/* Link the item to the list's tail: */
				item->pred=lruTail;
				if(lruTail!=0)
					lruTail->succ=item;
				else
					lruHead=item;
				lruTail=item;
				item->succ=0;
				}
			}
		else
			{
			/* Page out images from the head of the LRU list until there is enough space in memory: */
			while(usedMemory+item->imageSize>maxMemory)
				{
				/* Page out the current least-recently-used item: */
				LRUItem* out=lruHead;
				if(out->pageInCounter==requestCounter)
					{
					/* We just paged this image in; there is not enough memory in the cache to satisfy the request list: */
					Misc::throwStdErr("ImageCacher: Not enough memory to satisfy image request");
					}
				out->image->pageOut();
				out->pagedIn=false;
				usedMemory-=out->imageSize;
				
				/* Unlink the item from the LRU list: */
				if(out->succ!=0)
					out->succ->pred=0;
				else
					lruTail=0;
				lruHead=out->succ;
				out->succ=0;
				}
			
			/* Page in the requested image: */
			item->image->pageIn();
			item->pagedIn=true;
			item->pageInCounter=requestCounter;
			usedMemory+=item->imageSize;
			++numPageInRequests;
			
			/* Link the item to the list's tail: */
			item->pred=lruTail;
			if(lruTail!=0)
				lruTail->succ=item;
			else
				lruHead=item;
			lruTail=item;
			item->succ=0;
			}
		}
	
	/* Go to the next request transaction: */
	++requestCounter;
	}

/*******************************************************
Octree traversal functor class to colorize LiDAR points:
*******************************************************/

class NodeColorSampler
	{
	/* Elements: */
	private:
	LidarProcessOctree& lpo; // The processed LiDAR octree
	const std::vector<Image2D*>& images; // The vector of images
	ImageCacher imageCacher; // The cache manager for images
	Color* colorBuffer; // Array to hold colors for a node during processing
	Color* childColorBuffers[8]; // Array of color arrays for a node's children during subsampling
	LidarFile::Offset colorDataSize; // Size of each record in the color file
	LidarFile colorFile; // The file to which to write the color data
	size_t numProcessedNodes; // Number of already processed nodes
	size_t nextProgressUpdate; // Number of processed nodes at which the progress indicator should be updated
	size_t numAssignedColors; // Number of LiDAR points to which colors could be assigned
	
	/* Constructors and destructors: */
	public:
	NodeColorSampler(LidarProcessOctree& sLpo,const std::vector<Image2D*>& sImages,size_t imageMemorysize,const char* colorFileName); // Creates a color sampler with the given parameters
	~NodeColorSampler(void);
	
	/* Methods: */
	void operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel);
	size_t getNumAssignedColors(void) const // Returns the number of LiDAR points that were assigned colors
		{
		return numAssignedColors;
		}
	const ImageCacher& getImageCacher(void) const
		{
		return imageCacher;
		}
	};

/*********************************
Methods of class NodeColorSampler:
*********************************/

NodeColorSampler::NodeColorSampler(LidarProcessOctree& sLpo,const std::vector<Image2D*>& sImages,size_t imageMemorySize,const char* colorFileName)
	:lpo(sLpo),
	 images(sImages),
	 imageCacher(imageMemorySize),
	 colorBuffer(new Color[lpo.getMaxNumPointsPerNode()]),
	 colorDataSize(sizeof(Color)),
	 colorFile(colorFileName,LidarFile::ReadWrite),
	 numProcessedNodes(0),nextProgressUpdate((lpo.getNumNodes()+199)/200),
	 numAssignedColors(0)
	{
	/* Allocate the color subsampling arrays: */
	for(int i=0;i<8;++i)
		childColorBuffers[i]=new Color[lpo.getMaxNumPointsPerNode()];
	
	/* Write the color file's header: */
	colorFile.setEndianness(Misc::LittleEndian);
	LidarDataFileHeader dfh((unsigned int)(colorDataSize));
	dfh.write(colorFile);
	
	/* Register all images in the image list: */
	for(std::vector<Image2D*>::const_iterator iIt=images.begin();iIt!=images.end();++iIt)
		imageCacher.registerImage(*iIt);
	}

NodeColorSampler::~NodeColorSampler(void)
	{
	delete[] colorBuffer;
	for(int i=0;i<8;++i)
		delete[] childColorBuffers[i];
	}

namespace {

/**************
Helper classes:
**************/

class NodePointFinder // Class to find a point inside an octree node
	{
	/* Elements: */
	private:
	Point queryPoint; // The position of the point to find
	const LidarPoint* foundPoint; // The found LiDAR point
	
	/* Constructors and destructors: */
	public:
	NodePointFinder(const Point& sQueryPoint)
		:queryPoint(sQueryPoint),
		 foundPoint(0)
		{
		}
	
	/* Methods: */
	void operator()(const LidarPoint& lp)
		{
		if(lp==queryPoint)
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

void NodeColorSampler::operator()(LidarProcessOctree::Node& node,unsigned int nodeLevel)
	{
	if(node.isLeaf())
		{
		if(node.getNumPoints()>0)
			{
			/* Find all images whose bounding boxes overlap this node: */
			std::vector<Image2D*> nodeImages;
			for(std::vector<Image2D*>::const_iterator iIt=images.begin();iIt!=images.end();++iIt)
				{
				const Image2D::Box& box=(*iIt)->getWorldBox();
				if(box.min[0]-lpo.getOffset()[0]<node.getDomain().getMax()[0]&&box.max[0]-lpo.getOffset()[0]>node.getDomain().getMin()[0]
		  		 &&box.min[1]-lpo.getOffset()[1]<node.getDomain().getMax()[1]&&box.max[1]-lpo.getOffset()[1]>node.getDomain().getMin()[1])
					nodeImages.push_back(*iIt);
				}
			
			/* Page in all found images: */
			imageCacher.requestImages(nodeImages);
			
			/* Assign a color to each LiDAR point in this node: */
			for(unsigned int i=0;i<node.getNumPoints();++i)
				{
				/* Copy the point's original color: */
				colorBuffer[i]=node[i].value;
				
				/* Lookup the point in all images overlapping this node: */
				Geometry::Point<double,3> pos;
				for(int j=0;j<3;++j)
					pos[j]=double(node[i][j]+lpo.getOffset()[j]);
				
				// Evil hack afoot!
				if(pos[2]<0.0)
					{
					if(pos[2]>=-525.0)
						{
						colorBuffer[i][0]=Color::Scalar((pos[2]+525.0)*240.0/525.0);
						colorBuffer[i][1]=Color::Scalar((pos[2]+525.0)*95.0/525.0+160.0);
						colorBuffer[i][2]=Color::Scalar(255);
						}
					else if(pos[2]>=-1050.0)
						{
						colorBuffer[i][0]=Color::Scalar(0);
						colorBuffer[i][1]=Color::Scalar((pos[2]+1050.0)*160.0/525.0);
						colorBuffer[i][2]=Color::Scalar(255);
						}
					else if(pos[2]>=-1575.0)
						{
						colorBuffer[i][0]=Color::Scalar(0);
						colorBuffer[i][1]=Color::Scalar(0);
						colorBuffer[i][2]=Color::Scalar((pos[2]+1575.0)*191.0/525.0+64.0);
						}
					else
						{
						colorBuffer[i][0]=Color::Scalar(0);
						colorBuffer[i][1]=Color::Scalar(0);
						colorBuffer[i][2]=Color::Scalar(64);
						}
					colorBuffer[i][3]=Color::Scalar(255);
					
					++numAssignedColors;
					}
				else
				for(std::vector<Image2D*>::iterator iIt=nodeImages.begin();iIt!=nodeImages.end();++iIt)
					{
					Image2D::SampleResult sr=(*iIt)->getColor(Image2D::Point(pos.getComponents()));
					if(sr.first)
						{
						/* Copy the sample result: */
						for(int j=0;j<3;++j)
							colorBuffer[i][j]=sr.second[j];
						colorBuffer[i][3]=Color::Scalar(255);
						
						++numAssignedColors;
						
						/* Bail out: */
						break;
						}
					}
				}
			}
		}
	else
		{
		/* Get pointers to the node's children and load their color arrays: */
		colorFile.flush();
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			LidarProcessOctree::Node* child=lpo.getChild(&node,childIndex);
			if(child->getNumPoints()>0)
				{
				colorFile.setReadPosAbs(LidarDataFileHeader::getFileSize()+colorDataSize*child->getDataOffset());
				colorFile.read(childColorBuffers[childIndex],child->getNumPoints());
				}
			}
		
		/* Find the direct ancestors of all LiDAR points in this node and copy their color values from the child arrays: */
		for(unsigned int i=0;i<node.getNumPoints();++i)
			{
			/* Find the child node containing this point's ancestor: */
			int pointChildIndex=node.getDomain().findChild(node[i]);
			LidarProcessOctree::Node* pointChild=lpo.getChild(&node,pointChildIndex);
			
			/* Find the point's ancestor: */
			NodePointFinder npf(node[i]);
			lpo.processNodePointsDirected(pointChild,npf);
			if(npf.getFoundPoint()==0)
				{
				/* This is an internal corruption in the octree file. Print a helpful and non-offensive error message: */
				Misc::throwStdErr("Fatal error: Octree file corrupted around position (%f, %f, %f)",node[i][0],node[i][1],node[i][2]);
				}
			
			/* Retrieve the ancestor's color: */
			colorBuffer[i]=childColorBuffers[pointChildIndex][npf.getFoundPoint()-pointChild->getPoints()];
			}
		}
	
	/* Write the node's colors to the color file: */
	colorFile.setWritePosAbs(LidarDataFileHeader::getFileSize()+colorDataSize*node.getDataOffset());
	colorFile.write(colorBuffer,node.getNumPoints());
	
	/* Update the progress counter: */
	++numProcessedNodes;
	if(numProcessedNodes>=nextProgressUpdate)
		{
		int percent=int((numProcessedNodes*100+lpo.getNumNodes()/2)/lpo.getNumNodes());
		std::cout<<"\rAssigning colors... "<<std::setw(3)<<percent<<"%"<<std::flush;
		nextProgressUpdate=((percent+1)*lpo.getNumNodes()-lpo.getNumNodes()/2+99)/100;
		}
	}

int main(int argc,char* argv[])
	{
	const char* lidarFileName=0;
	size_t octreeCacheSize=512;
	size_t imageCacheSize=512;
	const char* colorFileName=0;
	std::vector<Image2D*> images;
	for(int i=1;i<argc;++i)
		{
		if(argv[i][0]=='-')
			{
			if(strcasecmp(argv[i]+1,"octreeCache")==0)
				{
				++i;
				octreeCacheSize=size_t(atoi(argv[i]));
				}
			else if(strcasecmp(argv[i]+1,"imageCache")==0)
				{
				++i;
				imageCacheSize=size_t(atoi(argv[i]));
				}
			}
		else if(lidarFileName==0)
			lidarFileName=argv[i];
		else if(colorFileName==0)
			colorFileName=argv[i];
		else
			{
			/* Load an image: */
			Image2D* newImage=new Image2D(argv[i]);
			images.push_back(newImage);
			}
		}
	if(lidarFileName==0)
		{
		std::cerr<<"No LiDAR file name provided"<<std::endl;
		return 1;
		}
	if(colorFileName==0)
		colorFileName="Colors";
	if(images.empty())
		{
		std::cerr<<"No images provided"<<std::endl;
		return 1;
		}
	
	/* Create a processing octree: */
	LidarProcessOctree lpo(lidarFileName,octreeCacheSize*size_t(1024*1024));
	std::cout<<"LiDAR coordinate offset: "<<lpo.getOffset()[0]<<", "<<lpo.getOffset()[1]<<", "<<lpo.getOffset()[2]<<std::endl;
	
	/* Assign colors to all points in the octree: */
	std::string lidarColorFileName=lidarFileName;
	lidarColorFileName.push_back('/');
	lidarColorFileName.append(colorFileName);
	NodeColorSampler nodeColorSampler(lpo,images,imageCacheSize*size_t(1024*1024),lidarColorFileName.c_str());
	std::cout<<"Assigning colors...   0%"<<std::flush;
	lpo.processNodesPostfix(nodeColorSampler);
	std::cout<<std::endl;
	std::cout<<nodeColorSampler.getNumAssignedColors()<<" LiDAR points re-colored"<<std::endl;
	std::cout<<nodeColorSampler.getImageCacher().getNumPageInRequests()<<" images paged into main memory during processing"<<std::endl;
	
	/* Delete all images: */
	for(std::vector<Image2D*>::iterator iIt=images.begin();iIt!=images.end();++iIt)
		delete *iIt;
	
	return 0;
	}
