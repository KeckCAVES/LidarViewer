/***********************************************************************
LidarProcessOctree - Class to process multiresolution LiDAR point sets.
Copyright (c) 2008-2009 Oliver Kreylos

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

#include <string>
#include <iostream>
#include <Misc/Utility.h>
#include <Misc/ThrowStdErr.h>
#include <Math/Math.h>

#include "LidarProcessOctree.h"

/*****************************************
Methods of class LidarProcessOctree::Node:
*****************************************/

LidarProcessOctree::Node::Node(void)
	:parent(0),children(0),points(0),
	 processCounter(0),
	 lruPred(0),lruSucc(0)
	{
	}

LidarProcessOctree::Node::~Node(void)
	{
	/* Delete point array: */
	delete[] points;
	
	/* Delete children: */
	delete[] children;
	}

/***********************************
Methods of class LidarProcessOctree:
***********************************/

void LidarProcessOctree::subdivide(LidarProcessOctree::Node& node)
	{
	{
	#if ALLOW_THREADING
	Threads::Mutex::Lock lruLock(lruMutex);
	#endif
	++numSubdivideCalls;
	}
	
	#if ALLOW_THREADING
	Threads::Mutex::Lock nodeLock(node.mutex);
	#endif
	
	/* Bail out if the node's children have magically appeared since the last check: */
	if(node.children!=0)
		return;
	
	/* Check if the node cache is full: */
	{
	#if ALLOW_THREADING
	Threads::Mutex::Lock lruLock(lruMutex);
	#endif
	if(numCachedNodes+8>cacheSize)
		{
		/* Find an unused leaf node parent and remove its children: */
		Node* coarsenNode;
		#if ALLOW_THREADING
		for(coarsenNode=lruHead;coarsenNode!=0;coarsenNode=coarsenNode->lruSucc)
			{
			/* Check if the leaf parent is currently unused: */
			coarsenNode->mutex.lock();
			if(coarsenNode->processCounter==0)
				break;
			coarsenNode->mutex.unlock();
			}
		#else
		for(coarsenNode=lruHead;coarsenNode!=0&&coarsenNode->processCounter!=0;coarsenNode=coarsenNode->lruSucc)
			;
		#endif
		
		/* Remove the found node's children: */
		delete[] coarsenNode->children;
		coarsenNode->children=0;
		
		#if ALLOW_THREADING
		coarsenNode->mutex.unlock();
		#endif
		
		/* Remove the found node from the lru list: */
		if(coarsenNode->lruPred!=0)
			coarsenNode->lruPred->lruSucc=coarsenNode->lruSucc;
		else
			lruHead=coarsenNode->lruSucc;
		if(coarsenNode->lruSucc!=0)
			coarsenNode->lruSucc->lruPred=coarsenNode->lruPred;
		else
			lruTail=coarsenNode->lruPred;
		coarsenNode->lruPred=0;
		coarsenNode->lruSucc=0;
		
		/* Update the node cache: */
		numCachedNodes-=8;
		
		/* Check if the found node's parent is now a leaf parent node: */
		Node* parent=coarsenNode->parent;
		if(parent!=0)
			{
			/* Check if all children of the parent are leaves: */
			bool leafParent=true;
			for(int childIndex=0;childIndex<8&&leafParent;++childIndex)
				leafParent=parent->children[childIndex].children==0;
			
			if(leafParent)
				{
				/* Add the parent to the end of the leaf parent node list: */
				parent->lruPred=lruTail;
				parent->lruSucc=0;
				if(lruTail!=0)
					lruTail->lruSucc=parent;
				else
					lruHead=parent;
				lruTail=parent;
				}
			}
		}
	
	/* Check if the node's parent was a leaf parent: */
	Node* parent=node.parent;
	if(parent!=0&&(parent->lruPred!=0||parent->lruSucc!=0||lruHead==parent))
		{
		/* Remove the node's parent from the leaf parent node list: */
		if(parent->lruPred!=0)
			parent->lruPred->lruSucc=parent->lruSucc;
		else
			lruHead=parent->lruSucc;
		if(parent->lruSucc!=0)
			parent->lruSucc->lruPred=parent->lruPred;
		else
			lruTail=parent->lruPred;
		parent->lruPred=0;
		parent->lruSucc=0;
		}
	
	/* Add the node to the leaf parent node list: */
	node.lruPred=lruTail;
	node.lruSucc=0;
	if(lruTail!=0)
		lruTail->lruSucc=&node;
	else
		lruHead=&node;
	lruTail=&node;
	
	/* Update the node cache: */
	numCachedNodes+=8U;
	}
	
	/* Create and load the node's children: */
	{
	#if ALLOW_THREADING
	Threads::Mutex::Lock fileLock(fileMutex);
	#endif
	Node* children=new Node[8];
	indexFile.seekSet(node.childrenOffset);
	for(int childIndex=0;childIndex<8;++childIndex)
		{
		Node& child=children[childIndex];
		
		/* Read the child node's structure: */
		LidarOctreeFileNode ofn;
		ofn.read(indexFile);
		child.parent=&node;
		child.childrenOffset=ofn.childrenOffset;
		child.domain=Cube(node.domain,childIndex);
		child.numPoints=ofn.numPoints;
		child.dataOffset=ofn.dataOffset;
		child.detailSize=ofn.detailSize;
		
		if(child.numPoints>0)
			{
			/* Load the child node's points: */
			child.points=new LidarPoint[maxNumPointsPerNode]; // Always allocate maximum to prevent memory fragmentation
			pointsFile.seekSet(LidarDataFileHeader::getFileSize()+pointsRecordSize*child.dataOffset);
			pointsFile.read(child.points,child.numPoints);
			}
		
		++numLoadedNodes;
		}
	
	/* Install the node's children array: */
	node.children=children;
	}
	}

namespace {

/***********************************************************
Helper functions to load LiDAR files given a base file name:
***********************************************************/

std::string getLidarPartFileName(const char* lidarFileName,const char* partFileName)
	{
	std::string result=lidarFileName;
	result.push_back('/');
	result.append(partFileName);
	return result;
	}

}

LidarProcessOctree::LidarProcessOctree(const char* lidarFileName,size_t sCacheSize)
	:indexFile(getLidarPartFileName(lidarFileName,"Index").c_str(),"rb",LidarFile::LittleEndian),
	 pointsFile(getLidarPartFileName(lidarFileName,"Points").c_str(),"rb",LidarFile::LittleEndian),
	 numSubdivideCalls(0),numLoadedNodes(0),
	 lruHead(0),lruTail(0)
	{
	/* Read the octree file header: */
	LidarOctreeFileHeader ofh(indexFile);
	
	/* Initialize the root node's domain: */
	root.domain=ofh.domain;
	
	/* Initialize the tree structure: */
	maxNumPointsPerNode=ofh.maxNumPointsPerNode;
	
	/* Calculate the memory and GPU cache sizes: */
	size_t memNodeSize=sizeof(Node)+size_t(maxNumPointsPerNode)*sizeof(LidarPoint);
	cacheSize=(unsigned int)(sCacheSize/memNodeSize);
	if(cacheSize==0U)
		Misc::throwStdErr("LidarProcessOctree::LidarProcessOctree: Specified memory cache size too small");
	std::cout<<"Cache size: "<<cacheSize<<" memory nodes"<<std::endl;
	
	/* Read the root node's structure: */
	LidarOctreeFileNode rootfn;
	rootfn.read(indexFile);
	root.childrenOffset=rootfn.childrenOffset;
	root.numPoints=rootfn.numPoints;
	root.dataOffset=rootfn.dataOffset;
	root.detailSize=rootfn.detailSize;
	
	/* Get the total number of nodes by dividing the index file's size by the size of one octree node: */
	indexFile.seekEnd(0);
	LidarFile::Offset indexFileSize=indexFile.tell();
	numNodes=size_t((indexFileSize-LidarOctreeFileHeader::getFileSize())/LidarFile::Offset(LidarOctreeFileNode::getFileSize()));
	
	/* Read the point file's header: */
	LidarDataFileHeader dfh(pointsFile);
	pointsRecordSize=LidarFile::Offset(dfh.recordSize);
	
	if(root.numPoints>0)
		{
		/* Load the root node's points: */
		root.points=new LidarPoint[maxNumPointsPerNode]; // Always allocate maximum to prevent memory fragmentation
		pointsFile.seekSet(LidarDataFileHeader::getFileSize()+pointsRecordSize*root.dataOffset);
		pointsFile.read(root.points,root.numPoints);
		}
	
	++numLoadedNodes;
	
	/* Initialize the node cache: */
	numCachedNodes=1U;
	}

LidarProcessOctree::~LidarProcessOctree(void)
	{
	}

Point LidarProcessOctree::getRootCenter(void) const
	{
	return root.domain.getCenter();
	}

Scalar LidarProcessOctree::getRootSize(void) const
	{
	Scalar size=Scalar(0);
	for(int i=0;i<3;++i)
		size+=root.domain.getMax()[i]-root.domain.getMin()[i];
	return size/Scalar(6);
	}
