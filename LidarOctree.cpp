/***********************************************************************
LidarOctree - Class to render multiresolution LiDAR point sets.
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

#define PARANOIA 0

#include <string.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <Math/Math.h>
#include <Math/Constants.h>
#include <Geometry/Ray.h>
#include <Geometry/Box.h>
#include <GL/gl.h>
#include <GL/GLVertexTemplates.h>
#include <GL/GLContextData.h>
#include <GL/GLExtensionManager.h>
#include <GL/Extensions/GLARBPointParameters.h>
#include <GL/Extensions/GLARBVertexBufferObject.h>
#include <GL/GLFrustum.h>

#include "CoarseningHeap.h"

#include "LidarOctree.h"

/**********************************
Methods of class LidarOctree::Node:
**********************************/

Scalar LidarOctree::Node::intersectRay(const LidarOctree::Ray& ray,Scalar coneAngle2,Scalar lambda1,Scalar lambda2) const
	{
	if(children!=0)
		{
		#if 1
		Scalar lambdaMin=lambda2;
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			/* Make a box for the child's domain: */
			Vector domainSize(children[childIndex].size);
			Geometry::Box<Scalar,3> domain(children[childIndex].center-domainSize,children[childIndex].center+domainSize);
			std::pair<Scalar,Scalar> lambdas=domain.getRayParameters(ray);
			if(lambdas.first<lambda1)
				lambdas.first=lambda1;
			if(lambdas.second>lambda2)
				lambdas.second=lambda2;
			if(lambdas.first<lambdas.second)
				{
				Scalar lambda=children[childIndex].intersectRay(ray,coneAngle2,lambdas.first,lambdas.second);
				if(lambda<lambdas.second&&lambdaMin>lambda)
					lambdaMin=lambda;
				}
			}
		return lambdaMin;
		#else
		/* Determine the child node containing the ray's start point and the ray parameters at which the ray intersects the child node's separating planes: */
		int childIndex=0x0;
		Scalar planeLambdas[3];
		for(int i=0;i<3;++i)
			{
			if(ray.getOrigin()[i]>=center[i])
				{
				childIndex|=0x1<<i;
				if(ray.getDirection()[i]<Scalar(0))
					planeLambdas[i]=(center[i]-ray.getOrigin()[i])/ray.getDirection()[i];
				else
					planeLambdas[i]=Scalar(-1);
				}
			else
				{
				if(ray.getDirection()[i]>Scalar(0))
					planeLambdas[i]=(center[i]-ray.getOrigin()[i])/ray.getDirection()[i];
				else
					planeLambdas[i]=Scalar(-1);
				}
			}
		
		/* Traverse the children in the order they are intersected by the ray: */
		while(lambda1<lambda2)
			{
			/* Determine the next plane crossing: */
			int planeIndex=-1;
			Scalar nextLambda=lambda2;
			for(int i=0;i<3;++i)
				if(planeLambdas[i]>lambda1&&planeLambdas[i]<nextLambda)
					{
					planeIndex=i;
					nextLambda=planeLambdas[i];
					}
			
			/* Check the current child: */
			Scalar childLambda=children[childIndex].intersectRay(ray,coneAngle2,lambda1,nextLambda);
			
			/* If the current child reported an intersection, return it: */
			if(childLambda<nextLambda)
				return childLambda;
			
			/* Otherwise, go to the next child: */
			childIndex^=0x1<<planeIndex;
			lambda1=nextLambda;
			}
		
		/* No intersections found: */
		return lambda2;
		#endif
		}
	else
		{
		/* Intersect the ray with all points in this node: */
		Scalar lambdaMin=lambda2;
		for(unsigned int i=0;i<numPoints;++i)
			{
			/* Check if the point is closer than the previous one and inside the cone: */
			Vector sp;
			for(int j=0;j<3;++j)
				sp[j]=points[i].position[j]-ray.getOrigin()[j];
			Scalar x=sp*ray.getDirection();
			if(x>Scalar(0)&&x<lambdaMin)
				{
				Scalar y2=Geometry::sqr(Geometry::cross(sp,ray.getDirection()));
				if(y2/Math::sqr(x)<=coneAngle2)
					lambdaMin=x;
				}
			}
		return lambdaMin;
		}
	}

/**************************************
Methods of class LidarOctree::DataItem:
**************************************/

LidarOctree::DataItem::DataItem(unsigned int sCacheSize)
	:hasVertexBufferObjectExtension(GLARBVertexBufferObject::isSupported()),
	 cacheSize(sCacheSize-1),
	 cacheSlots(new CacheSlot[cacheSize]),
	 cacheNodeMap(cacheSize+cacheSize/4),
	 lruHead(0),lruTail(0),
	 bypassVertexBufferObjectId(0),
	 hasPointParametersExtension(GLARBPointParameters::isSupported()),
	 numRenderedNodes(0),numCacheMisses(0),numCacheBypasses(0),numRenderedPoints(0),numBypassedPoints(0)
	{
	/* Allocate the cache objects: */
	if(hasVertexBufferObjectExtension)
		{
		/* Initialize the vertex buffer object extension: */
		GLARBVertexBufferObject::initExtension();
		
		/* Create vertex buffer objects: */
		GLuint* vertexBufferObjectIds=new GLuint[cacheSize+1];
		glGenBuffersARB(cacheSize+1,vertexBufferObjectIds);
		for(unsigned int i=0;i<cacheSize;++i)
			cacheSlots[i].vertexBufferObjectId=vertexBufferObjectIds[i];
		bypassVertexBufferObjectId=vertexBufferObjectIds[cacheSize];
		delete[] vertexBufferObjectIds;
		}
	
	/* Initialize the LRU cache slot list: */
	lruHead=&cacheSlots[0];
	cacheSlots[0].pred=0;
	for(int i=1;i<cacheSize;++i)
		{
		cacheSlots[i-1].succ=&cacheSlots[i];
		cacheSlots[i].pred=&cacheSlots[i-1];
		}
	cacheSlots[cacheSize-1].succ=0;
	lruTail=&cacheSlots[cacheSize-1];
	
	if(hasPointParametersExtension)
		{
		/* Initialize the point parameters extension: */
		GLARBPointParameters::initExtension();
		}
	}

LidarOctree::DataItem::~DataItem(void)
	{
	if(hasVertexBufferObjectExtension)
		{
		/* Delete vertex buffer objects: */
		GLuint* vertexBufferObjectIds=new GLuint[cacheSize+1];
		for(unsigned int i=0;i<cacheSize;++i)
			vertexBufferObjectIds[i]=cacheSlots[i].vertexBufferObjectId;
		vertexBufferObjectIds[cacheSize]=bypassVertexBufferObjectId;
		glDeleteBuffersARB(cacheSize+1,vertexBufferObjectIds);
		delete[] vertexBufferObjectIds;
		}
	
	/* Delete the cache slots: */
	delete[] cacheSlots;
	}

/****************************
Methods of class LidarOctree:
****************************/

void LidarOctree::renderSubTree(const LidarOctree::Node* node,const LidarOctree::Frustum& frustum,LidarOctree::DataItem* dataItem) const
	{
	/* Bail out if the node is empty: */
	if(node->numPoints==0)
		return;
	
	/* Check if this node intersects the view frustum: */
	bool inside=true;
	for(int plane=0;plane<6&&inside;++plane)
		{
		const Frustum::Plane::Vector& normal=frustum.getFrustumPlane(plane).getNormal();
		
		/* Find the point on the node's bounding box which is closest to the frustum plane: */
		Point p=node->center;
		for(int i=0;i<3;++i)
			{
			if(normal[i]>Scalar(0))
				p[i]+=node->size;
			else
				p[i]-=node->size;
			}
		
		/* Check if the point is inside the view frustum: */
		inside=!frustum.getFrustumPlane(plane).contains(p);
		}
	
	/* Bail out if this node is not inside the view frustum: */
	if(!inside)
		return;
	
	/* Calculate the node's projected detail size: */
	Scalar projectedDetailSize=frustum.calcProjectedRadius(node->center,node->detailSize);
	if(projectedDetailSize>=Scalar(0))
		{
		/* Adjust the projected detail size based on the focus+context rule: */
		if(fncWeight>Scalar(0))
			{
			/* Calculate the distance from the node to the focus region: */
			Scalar fncDist=Geometry::dist(node->center,fncCenter)-node->radius;
			if(fncDist>fncRadius)
				{
				/* Adjust the node's projected detail size: */
				projectedDetailSize*=Math::pow(fncRadius/fncDist,fncWeight);
				}
			}
		}
	else
		projectedDetailSize=Math::Constants<Scalar>::max;
	
	/* Update node's rendering traversal state: */
	{
	Threads::Mutex::Lock nodeLock(node->nodeMutex);
	bool nodeStateChanged=false;
	if(node->renderPass!=renderPass)
		{
		/* Store the node's LOD value from this rendering pass: */
		node->maxLOD=projectedDetailSize;
		node->renderPass=renderPass;
		nodeStateChanged=true;
		}
	else
		{
		/* Update the node's maximum LOD value: */
		if(node->maxLOD<projectedDetailSize)
			{
			node->maxLOD=projectedDetailSize;
			nodeStateChanged=true;
			}
		}
	
	if(node->coarseningHeapIndex!=~0x0U&&nodeStateChanged)
		{
		/* Update the node's data in the coarsening heap: */
		{
		Threads::Mutex::Lock coarseningHeapLock(coarseningHeapMutex);
		coarseningHeap->move(const_cast<Node*>(node));
		}
		}
	}
	
	/* Check whether to render this node or its children: */
	if(projectedDetailSize>=maxRenderLOD)
		{
		if(node->children!=0)
			{
			/* Use the view point for a back-to-front traversal of the child nodes: */
			int childIndex=0x0;
			Point eye=frustum.getEye().toPoint();
			for(int i=0;i<3;++i)
				if(eye[i]>=node->center[i])
					childIndex|=0x1<<i;
			
			/* Render the node's children: */
			// for(int i=7;i>=0;--i) // Back-to-front rendering
			for(int i=0;i<8;++i) // Front-to-back rendering
				renderSubTree(&node->children[i^childIndex],frustum,dataItem);
			
			/* Done here... */
			return;
			}
		else if(node->childrenOffset!=Misc::LargeFile::Offset(0))
			{
			/* Try inserting this node into the node loader thread's request queue: */
			Threads::Mutex::Lock loadRequestLock(loadRequestMutex);
			if(node->subdivisionQueueIndex<subdivisionRequestQueueLength)
				{
				/* Update the node's position in the subdivision queue if necessary: */
				if(subdivisionRequestQueue[node->subdivisionQueueIndex].LOD<projectedDetailSize)
					{
					/* Adjust the node's position in the queue: */
					for(;node->subdivisionQueueIndex>0&&subdivisionRequestQueue[node->subdivisionQueueIndex-1].LOD<projectedDetailSize;--node->subdivisionQueueIndex)
						{
						subdivisionRequestQueue[node->subdivisionQueueIndex]=subdivisionRequestQueue[node->subdivisionQueueIndex-1];
						subdivisionRequestQueue[node->subdivisionQueueIndex].node->subdivisionQueueIndex=node->subdivisionQueueIndex;
						}
					
					/* Update the node's LOD value in the queue: */
					subdivisionRequestQueue[node->subdivisionQueueIndex].node=const_cast<Node*>(node);
					subdivisionRequestQueue[node->subdivisionQueueIndex].LOD=projectedDetailSize;
					}
				}
			else if(node->subdivisionQueueIndex==subdivisionRequestQueueLength&&subdivisionRequestQueue[subdivisionRequestQueueLength-1].LOD<projectedDetailSize)
				{
				/* Remove the node from the last queue slot: */
				if(subdivisionRequestQueue[subdivisionRequestQueueLength-1].node!=0)
					subdivisionRequestQueue[subdivisionRequestQueueLength-1].node->subdivisionQueueIndex=subdivisionRequestQueueLength;
				
				/* Insert the node into the subdivision request queue: */
				unsigned int insertIndex;
				for(insertIndex=subdivisionRequestQueueLength-1;insertIndex>0&&subdivisionRequestQueue[insertIndex-1].LOD<projectedDetailSize;--insertIndex)
					{
					if(subdivisionRequestQueue[insertIndex-1].node!=0)
						subdivisionRequestQueue[insertIndex-1].node->subdivisionQueueIndex=insertIndex;
					if(insertIndex<subdivisionRequestQueueLength)
						subdivisionRequestQueue[insertIndex]=subdivisionRequestQueue[insertIndex-1];
					}
				node->subdivisionQueueIndex=insertIndex;
				subdivisionRequestQueue[insertIndex].node=const_cast<Node*>(node);
				subdivisionRequestQueue[insertIndex].LOD=projectedDetailSize;
				
				/* Wake up the node loader thread: */
				loadRequestCond.signal();
				}
			}
		}
	
	/* Retrieve the cache slot containing this node's tile vertex array and image: */
	GLuint vertexBufferObjectId=0;
	DataItem::CacheSlot* slot=0;
	bool mustUploadData=false;
	DataItem::NodeHasher::Iterator cnIt=dataItem->cacheNodeMap.findEntry(node);
	if(cnIt.isFinished()) // Cache miss
		{
		/*******************************************************************
		We use an LRU cache replacement strategy, but do not replace cache
		slots that were already used in the same rendering pass to prevent
		cache thrashing. If the cache is full for a rendering pass,
		additional nodes will be rendered straight from the main memory data
		structures.
		*******************************************************************/
		
		/* Check if the head of the LRU list is from a previous rendering pass: */
		if(dataItem->lruHead->lastUsed!=renderPass)
			{
			/* Replace the cache slot: */
			slot=dataItem->lruHead;
			if(slot->node!=0)
				{
				/* Remove the currently cached node: */
				dataItem->cacheNodeMap.removeEntry(slot->node);
				}
			slot->node=node;
			dataItem->cacheNodeMap.setEntry(DataItem::NodeHasher::Entry(node,slot));
			
			/* Remember to upload the node's data later: */
			mustUploadData=true;
			}
		else
			++dataItem->numCacheBypasses;
		
		++dataItem->numCacheMisses;
		}
	else // Cache hit
		{
		/* Use the previously used cache slot: */
		slot=cnIt->getDest();
		
		/* Check if the cached point set is out-of-date: */
		if(slot->version!=node->pointsVersion)
			{
			/* Remember to upload the node's data later: */
			mustUploadData=true;
			}
		}
	
	if(slot!=0)
		{
		/* Install the cache slot's buffer object (and upload data if required): */
		if(dataItem->hasVertexBufferObjectExtension)
			{
			vertexBufferObjectId=slot->vertexBufferObjectId;
			glBindBufferARB(GL_ARRAY_BUFFER_ARB,vertexBufferObjectId);
			if(mustUploadData)
				glBufferDataARB(GL_ARRAY_BUFFER_ARB,node->numPoints*sizeof(Vertex),node->points,GL_DYNAMIC_DRAW_ARB);
			}
		
		/* Mark the cache slot as used and move it to the end of the LRU list: */
		slot->version=node->pointsVersion;
		slot->lastUsed=renderPass;
		if(slot->pred!=0)
			slot->pred->succ=slot->succ;
		else
			dataItem->lruHead=slot->succ;
		if(slot->succ!=0)
			slot->succ->pred=slot->pred;
		else
			dataItem->lruTail=slot->pred;
		slot->pred=dataItem->lruTail;
		dataItem->lruTail->succ=slot;
		slot->succ=0;
		dataItem->lruTail=slot;
		}
	else
		{
		/* Prepare for rendering directly from main memory: */
		if(dataItem->hasVertexBufferObjectExtension)
			{
			vertexBufferObjectId=dataItem->bypassVertexBufferObjectId;
			glBindBufferARB(GL_ARRAY_BUFFER_ARB,vertexBufferObjectId);
			glBufferDataARB(GL_ARRAY_BUFFER_ARB,node->numPoints*sizeof(Vertex),node->points,GL_STREAM_DRAW_ARB);
			}
		dataItem->numBypassedPoints+=node->numPoints;
		}
	
	/* Render this node's point set: */
	if(vertexBufferObjectId!=0)
		{
		/* Render from the vertex buffer: */
		glVertexPointer(static_cast<const Vertex*>(0));
		}
	else
		{
		/* Render straight from the node vertex array (only happens if GL_ARB_vertex_buffer_object not supported): */
		glVertexPointer(node->points);
		}
	glDrawArrays(GL_POINTS,0,node->numPoints);
	
	#if 0
	float s=float(node->size);
	glColor3f(1.0f,0.0f,0.0f);
	glBegin(GL_LINE_STRIP);
	glVertex3f(node->center[0]-s,node->center[1]-s,node->center[2]-s);
	glVertex3f(node->center[0]+s,node->center[1]-s,node->center[2]-s);
	glVertex3f(node->center[0]+s,node->center[1]+s,node->center[2]-s);
	glVertex3f(node->center[0]-s,node->center[1]+s,node->center[2]-s);
	glVertex3f(node->center[0]-s,node->center[1]-s,node->center[2]-s);
	glVertex3f(node->center[0]-s,node->center[1]-s,node->center[2]+s);
	glVertex3f(node->center[0]+s,node->center[1]-s,node->center[2]+s);
	glVertex3f(node->center[0]+s,node->center[1]+s,node->center[2]+s);
	glVertex3f(node->center[0]-s,node->center[1]+s,node->center[2]+s);
	glVertex3f(node->center[0]-s,node->center[1]-s,node->center[2]+s);
	glEnd();
	glBegin(GL_LINES);
	glVertex3f(node->center[0]+s,node->center[1]-s,node->center[2]-s);
	glVertex3f(node->center[0]+s,node->center[1]-s,node->center[2]+s);
	glVertex3f(node->center[0]+s,node->center[1]+s,node->center[2]-s);
	glVertex3f(node->center[0]+s,node->center[1]+s,node->center[2]+s);
	glVertex3f(node->center[0]-s,node->center[1]+s,node->center[2]-s);
	glVertex3f(node->center[0]-s,node->center[1]+s,node->center[2]+s);
	glEnd();
	#endif
	
	++dataItem->numRenderedNodes;
	dataItem->numRenderedPoints+=node->numPoints;
	}

void LidarOctree::interactWithSubTree(LidarOctree::Node* node,const LidarOctree::Interactor& interactor)
	{
	/* Bail out if the node is empty: */
	if(node->numPoints==0)
		return;
	
	/* Calculate the node's distance from the interactor: */
	Scalar interactorDist2=Geometry::sqrDist(interactor.center,node->center);
	if(interactorDist2>=Math::sqr(interactor.radius+node->radius))
		return;
	
	/* Calculate an appropriate LOD value for the node: */
	Scalar interactorLOD=(node->radius*maxRenderLOD)/interactor.radius; // This could use some tweaking
	
	/* Update node's rendering traversal state: */
	{
	Threads::Mutex::Lock nodeLock(node->nodeMutex);
	bool nodeStateChanged=false;
	if(node->renderPass!=renderPass)
		{
		/* Store the node's LOD value from this rendering pass: */
		node->maxLOD=interactorLOD;
		node->renderPass=renderPass;
		nodeStateChanged=true;
		}
	else
		{
		/* Update the node's maximum LOD value: */
		if(node->maxLOD<interactorLOD)
			{
			node->maxLOD=interactorLOD;
			nodeStateChanged=true;
			}
		}
	
	if(node->coarseningHeapIndex!=~0x0U&&nodeStateChanged)
		{
		/* Update the node's data in the coarsening heap: */
		{
		Threads::Mutex::Lock coarseningHeapLock(coarseningHeapMutex);
		coarseningHeap->move(node);
		}
		}
	}
	
	/* Check if the node should be subdivided: */
	if(node->children==0)
		{
		/* Subdivide if the node is not a leaf: */
		if(node->childrenOffset!=Misc::LargeFile::Offset(0))
			{
			/* Try inserting this node into the node loader thread's request queue: */
			Threads::Mutex::Lock loadRequestLock(loadRequestMutex);
			if(node->subdivisionQueueIndex<subdivisionRequestQueueLength)
				{
				/* Update the node's position in the subdivision queue if necessary: */
				if(subdivisionRequestQueue[node->subdivisionQueueIndex].LOD<interactorLOD)
					{
					/* Adjust the node's position in the queue: */
					for(;node->subdivisionQueueIndex>0&&subdivisionRequestQueue[node->subdivisionQueueIndex-1].LOD<interactorLOD;--node->subdivisionQueueIndex)
						{
						subdivisionRequestQueue[node->subdivisionQueueIndex]=subdivisionRequestQueue[node->subdivisionQueueIndex-1];
						subdivisionRequestQueue[node->subdivisionQueueIndex].node->subdivisionQueueIndex=node->subdivisionQueueIndex;
						}
					
					/* Update the node's LOD value in the queue: */
					subdivisionRequestQueue[node->subdivisionQueueIndex].node=node;
					subdivisionRequestQueue[node->subdivisionQueueIndex].LOD=interactorLOD;
					}
				}
			else if(node->subdivisionQueueIndex==subdivisionRequestQueueLength&&subdivisionRequestQueue[subdivisionRequestQueueLength-1].LOD<interactorLOD)
				{
				/* Remove the node from the last queue slot: */
				if(subdivisionRequestQueue[subdivisionRequestQueueLength-1].node!=0)
					subdivisionRequestQueue[subdivisionRequestQueueLength-1].node->subdivisionQueueIndex=subdivisionRequestQueueLength;
				
				/* Insert the node into the subdivision request queue: */
				unsigned int insertIndex;
				for(insertIndex=subdivisionRequestQueueLength-1;insertIndex>0&&subdivisionRequestQueue[insertIndex-1].LOD<interactorLOD;--insertIndex)
					{
					if(subdivisionRequestQueue[insertIndex-1].node!=0)
						subdivisionRequestQueue[insertIndex-1].node->subdivisionQueueIndex=insertIndex;
					if(insertIndex<subdivisionRequestQueueLength)
						subdivisionRequestQueue[insertIndex]=subdivisionRequestQueue[insertIndex-1];
					}
				node->subdivisionQueueIndex=insertIndex;
				subdivisionRequestQueue[insertIndex].node=node;
				subdivisionRequestQueue[insertIndex].LOD=interactorLOD;
				
				/* Wake up the node loader thread: */
				loadRequestCond.signal();
				}
			}
		}
	else
		{
		/* Recurse into the node's children: */
		for(int i=0;i<8;++i)
			interactWithSubTree(&node->children[i],interactor);
		}
	
	}

void LidarOctree::selectPoints(LidarOctree::Node* node,const LidarOctree::Interactor& interactor)
	{
	/* Check if the interactor's region of influence intersects the node's domain: */
	Scalar dist2=Scalar(0);
	for(int i=0;i<3;++i)
		{
		Scalar d;
		if((d=interactor.center[i]-(node->center[i]+node->size))>Scalar(0))
			dist2+=Math::sqr(d);
		else if((d=interactor.center[i]-(node->center[i]-node->size))<Scalar(0))
			dist2+=Math::sqr(d);
		}
	Scalar ir2=Math::sqr(interactor.radius);
	if(dist2<ir2)
		{
		/* Select all points inside the interactor's region of influence in this node: */
		bool pointsChanged=false;
		for(unsigned int i=0;i<node->numPoints;++i)
			{
			Scalar pdist2=Geometry::sqrDist(interactor.center,Point(node->points[i].position.getXyzw()));
			if(pdist2<ir2)
				{
				/* Create a selection mask if there is none already: */
				if(node->selectedPoints==0)
					{
					node->selectedPoints=new bool[maxNumPointsPerNode];
					for(unsigned int i=0;i<node->numPoints;++i)
						node->selectedPoints[i]=false;
					node->selectedPointColors=new Color[maxNumPointsPerNode];
					}
				
				/* Select this point: */
				if(!node->selectedPoints[i])
					{
					node->selectedPoints[i]=true;
					Color& col=node->points[i].color;
					node->selectedPointColors[i]=col;
					float intensity=float(col[0])*0.299f+float(col[1])*0.587f+float(col[2])*0.114f;
					if(intensity<127.5f)
						{
						col[0]=GLubyte(0);
						col[1]=GLubyte(intensity+127.5f);
						col[2]=GLubyte(0);
						}
					else
						{
						col[0]=GLubyte(intensity-127.5f);
						col[1]=GLubyte(255);
						col[2]=GLubyte(intensity-127.5f);
						}
					pointsChanged=true;
					}
				}
			}
		
		/* Check if the points array has to be invalidated: */
		if(pointsChanged)
			++node->pointsVersion;
		
		if(node->children!=0)
			{
			/* Recurse into the node's children: */
			for(int childIndex=0;childIndex<8;++childIndex)
				selectPoints(&node->children[childIndex],interactor);
			}
		}
	}

void LidarOctree::deselectPoints(LidarOctree::Node* node,const LidarOctree::Interactor& interactor)
	{
	/* Check if the interactor's region of influence intersects the node's domain: */
	Scalar dist2=Scalar(0);
	for(int i=0;i<3;++i)
		{
		Scalar d;
		if((d=interactor.center[i]-(node->center[i]+node->size))>Scalar(0))
			dist2+=Math::sqr(d);
		else if((d=interactor.center[i]-(node->center[i]-node->size))<Scalar(0))
			dist2+=Math::sqr(d);
		}
	Scalar ir2=Math::sqr(interactor.radius);
	if(dist2<ir2)
		{
		if(node->selectedPoints!=0)
			{
			/* Deselect all points inside the interactor's region of influence in this node: */
			bool pointsChanged=false;
			bool hasSelectedPoints=false;
			for(unsigned int i=0;i<node->numPoints;++i)
				{
				Scalar pdist2=Geometry::sqrDist(interactor.center,Point(node->points[i].position.getXyzw()));
				if(pdist2<ir2)
					{
					/* Deselect this point: */
					if(node->selectedPoints[i])
						{
						node->selectedPoints[i]=false;
						node->points[i].color=node->selectedPointColors[i];
						pointsChanged=true;
						}
					}
				hasSelectedPoints=hasSelectedPoints||node->selectedPoints[i];
				}
			
			/* Destroy the selection mask if there are no selected points: */
			if(!hasSelectedPoints)
				{
				delete[] node->selectedPoints;
				node->selectedPoints=0;
				delete[] node->selectedPointColors;
				node->selectedPointColors=0;
				}
			
			/* Check if the points array has to be invalidated: */
			if(pointsChanged)
				++node->pointsVersion;
			}
		
		if(node->children!=0)
			{
			/* Recurse into the node's children: */
			for(int childIndex=0;childIndex<8;++childIndex)
				deselectPoints(&node->children[childIndex],interactor);
			}
		}
	}

void LidarOctree::clearSelection(LidarOctree::Node* node)
	{
	if(node->selectedPoints!=0)
		{
		/* Deselect all selected points: */
		bool pointsChanged=false;
		for(unsigned int i=0;i<node->numPoints;++i)
			{
			if(node->selectedPoints[i])
				{
				GLubyte intensity=node->selectedPoints[i];
				node->points[i].color=node->selectedPointColors[i];
				pointsChanged=true;
				}
			}
		
		/* Destroy the selection mask: */
		delete[] node->selectedPoints;
		node->selectedPoints=0;
		delete[] node->selectedPointColors;
		node->selectedPointColors=0;
		
		/* Check if the points array has to be invalidated: */
		if(pointsChanged)
			++node->pointsVersion;
		}
	
	if(node->children!=0)
		{
		/* Recurse into the node's children: */
		for(int childIndex=0;childIndex<8;++childIndex)
			clearSelection(&node->children[childIndex]);
		}
	}

void LidarOctree::loadNodePoints(LidarOctree::Node* node)
	{
	/* Create the node's point array (always allocate the maximum size to prevent memory fragmentation): */
	node->points=new Vertex[maxNumPointsPerNode];
	
	/* Load the node's points into a temporary point buffer: */
	LidarPoint* pointsBuffer=new LidarPoint[maxNumPointsPerNode];
	pointsFile.seekSet(node->pointsOffset);
	pointsFile.read(pointsBuffer,node->numPoints);
	
	/* Convert the LiDAR points to render points: */
	for(unsigned int i=0;i<node->numPoints;++i)
		{
		/* Copy the point color: */
		node->points[i].color=pointsBuffer[i].value;
		
		#if RECENTER_OCTREE
		/* Offset the points so that the root node's center is the origin: */
		for(int j=0;j<3;++j)
			node->points[i].position[j]=Vertex::Position::Scalar(pointsBuffer[i][j]-dataCenter[j]);
		#else
		/* Copy the points: */
		for(int j=0;j<3;++j)
			node->points[i].position[j]=Vertex::Position::Scalar(pointsBuffer[i][j]);
		#endif
		}
	
	/* Clean up: */
	delete[] pointsBuffer;
	}

void* LidarOctree::nodeLoaderThreadMethod(void)
	{
	while(true)
		{
		/* Get the next subdivision request: */
		Node* node;
		{
		Threads::Mutex::Lock loadRequestLock(loadRequestMutex);
		while(subdivisionRequestQueue[0].node==0)
			loadRequestCond.wait(loadRequestMutex);
		node=subdivisionRequestQueue[0].node;
		
		/* Remove the node from the queue: */
		unsigned int i;
		for(i=1;i<subdivisionRequestQueueLength&&subdivisionRequestQueue[i].node!=0;++i)
			{
			subdivisionRequestQueue[i].node->subdivisionQueueIndex=i-1;
			subdivisionRequestQueue[i-1]=subdivisionRequestQueue[i];
			}
		subdivisionRequestQueue[i-1].node=0;
		subdivisionRequestQueue[i-1].LOD=Scalar(0);
		
		/* Lock the node until it has been loaded: */
		node->subdivisionQueueIndex=~0x0U;
		}
		
		/* Create the node's children: */
		Node* children=new Node[8];
		octreeFile.seekSet(node->childrenOffset);
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			LidarOctreeFileNode ofn;
			ofn.read(octreeFile);
			children[childIndex].parent=node;
			children[childIndex].childrenOffset=ofn.childrenOffset;
			children[childIndex].center=node->center;
			children[childIndex].size=Math::div2(node->size);
			for(int i=0;i<3;++i)
				if(childIndex&(1<<i))
					children[childIndex].center[i]+=children[childIndex].size;
				else
					children[childIndex].center[i]-=children[childIndex].size;
			children[childIndex].radius=Math::div2(node->radius);
			children[childIndex].numPoints=ofn.numPoints;
			children[childIndex].pointsOffset=ofn.pointsOffset;
			children[childIndex].detailSize=ofn.detailSize;
			children[childIndex].subdivisionQueueIndex=subdivisionRequestQueueLength;
			}
		
		/* Load the node's children's point data: */
		for(int childIndex=0;childIndex<8;++childIndex)
			loadNodePoints(&children[childIndex]);
		
		/* Add the new child nodes to the ready list: */
		{
		Threads::Mutex::Lock readyNodesLock(readyNodesMutex);
		readyNodes.push_back(children);
		}
		
		/* Notify a user program: */
		{
		Threads::Mutex::Lock treeUpdateFunctionLock(treeUpdateFunctionMutex);
		if(treeUpdateFunction!=0)
			treeUpdateFunction(treeUpdateFunctionArg);
		}
		}
	
	return 0;
	}

namespace {

/************************************************************
Helper functions to load octree files given a base file name:
************************************************************/

const char* getBaseNameExtension(const char* baseFileName)
	{
	/* Extract the base file name's extension: */
	const char* extPtr=0;
	const char* bfnPtr;
	for(bfnPtr=baseFileName;*bfnPtr!='\0';++bfnPtr)
		if(*bfnPtr=='.')
			extPtr=bfnPtr;
	if(extPtr==0)
		extPtr=bfnPtr;
	
	/* Check if this is the extension we're looking for: */
	if(strcasecmp(extPtr,".o")!=0&&strcasecmp(extPtr,".oct")!=0&&strcasecmp(extPtr,".obin")!=0)
		extPtr=bfnPtr;
	
	return extPtr;
	}

std::string createOctreeFileName(const char* baseFileName)
	{
	/* Copy the given base name without extension: */
	std::string result=std::string(baseFileName,getBaseNameExtension(baseFileName));
	
	/* Append the proper extension: */
	result.append(".oct");
	
	return result;
	}

std::string createPointsFileName(const char* baseFileName)
	{
	/* Copy the given base name without extension: */
	std::string result=std::string(baseFileName,getBaseNameExtension(baseFileName));
	
	/* Append the proper extension: */
	result.append(".obin");
	
	return result;
	}

}

LidarOctree::LidarOctree(const char* octreeFileName,unsigned int sCacheSize,unsigned int sGlCacheSize)
	:octreeFile(createOctreeFileName(octreeFileName).c_str(),"rb",Misc::LargeFile::LittleEndian),
	 pointsFile(createPointsFileName(octreeFileName).c_str(),"rb",Misc::LargeFile::LittleEndian),
	 maxRenderLOD(Math::sqrt(Scalar(2))),
	 fncWeight(0),
	 numCachedNodes(0),
	 renderPass(0U),
	 treeUpdateFunction(0),treeUpdateFunctionArg(0),
	 subdivisionRequestQueueLength(4),subdivisionRequestQueue(new SubdivisionRequest[subdivisionRequestQueueLength]),
	 coarseningHeap(0)
	{
	/* Read the octree file header: */
	LidarOctreeFileHeader ofh;
	ofh.read(octreeFile);
	
	/* Initialize the tree structure: */
	maxNumPointsPerNode=ofh.maxNumPointsPerNode;
	
	/* Calculate the memory and GPU cache sizes: */
	size_t memNodeSize=sizeof(Node)+maxNumPointsPerNode*sizeof(Vertex);
	size_t glNodeSize=maxNumPointsPerNode*sizeof(Vertex);
	cacheSize=sCacheSize/memNodeSize;
	if(cacheSize==0)
		Misc::throwStdErr("LidarOctree::LidarOctree: Specified memory cache size too small");
	glCacheSize=sGlCacheSize/glNodeSize;
	if(glCacheSize==0)
		Misc::throwStdErr("LidarOctree::LidarOctree: Specified GPU cache size too small");
	std::cout<<"Cache sizes: "<<cacheSize<<" memory nodes, "<<glCacheSize<<" GPU nodes"<<std::endl;
	
	/* Read the root node's structure: */
	LidarOctreeFileNode rootfn;
	rootfn.read(octreeFile);
	root.childrenOffset=rootfn.childrenOffset;
	
	/* Initialize the root node's domain: */
	#if RECENTER_OCTREE
	dataCenter=ofh.center;
	root.center=Point::origin;
	#else
	dataCenter=Point::origin;
	root.center=ofh.center;
	#endif
	root.size=ofh.radius;
	root.radius=root.size*Math::sqrt(Scalar(3));
	
	root.numPoints=rootfn.numPoints;
	root.pointsOffset=rootfn.pointsOffset;
	root.detailSize=rootfn.detailSize;
	root.subdivisionQueueIndex=subdivisionRequestQueueLength;
	
	/* Load the root's point data: */
	loadNodePoints(&root);
	numCachedNodes=1;
	
	/* Create the coarsening heap: */
	coarseningHeap=new CoarseningHeap<Node>((cacheSize+7U)/8U+1U); // Pessimistic estimate
	
	/* Start the node loader thread: */
	nodeLoaderThread.start(this,&LidarOctree::nodeLoaderThreadMethod);
	
	/* Initialize render state: */
	renderPass=1U;
	}

LidarOctree::~LidarOctree(void)
	{
	/* Shut down the node loader thread: */
	nodeLoaderThread.cancel();
	nodeLoaderThread.join();
	
	/* Delete the coarsening heap: */
	delete coarseningHeap;
	
	/* Delete the subdivision request queue: */
	delete[] subdivisionRequestQueue;
	}

void LidarOctree::initContext(GLContextData& contextData) const
	{
	/* Create a context data item: */
	DataItem* dataItem=new DataItem(glCacheSize);
	contextData.addDataItem(this,dataItem);
	}

void LidarOctree::setTreeUpdateFunction(LidarOctree::TreeUpdateFunction newTreeUpdateFunction,void* newTreeUpdateFunctionArg)
	{
	{
	Threads::Mutex::Lock treeUpdateFunctionLock(treeUpdateFunctionMutex);
	
	/* Install the new tree update function: */
	treeUpdateFunctionArg=newTreeUpdateFunctionArg;
	treeUpdateFunction=newTreeUpdateFunction;
	}
	}

void LidarOctree::setRenderQuality(Scalar qualityLevel)
	{
	/* Calculate the maximum render LOD based on the optimum and the quality level: */
	maxRenderLOD=Math::sqrt(Scalar(2));
	maxRenderLOD*=Math::pow(Scalar(0.5),qualityLevel);
	}

void LidarOctree::setFocusAndContext(const Point& newFncCenter,Scalar newFncRadius,Scalar newFncWeight)
	{
	fncCenter=newFncCenter;
	fncRadius=newFncRadius;
	fncWeight=newFncWeight;
	}

void LidarOctree::startRenderPass(void)
	{
	#if PARANOIA
	/* Check the coarsening heap: */
	if(!coarseningHeap->checkHeap())
		std::cerr<<"Inconsistent coarsening heap!"<<std::endl;
	#endif
	
	{
	/* Process all nodes recently loaded by the node loader thread: */
	Threads::Mutex::Lock readyNodesLock(readyNodesMutex);
	for(std::vector<Node*>::iterator rnIt=readyNodes.begin();rnIt!=readyNodes.end();++rnIt)
		{
		Node* children=*rnIt;
		Node* node=children->parent;
		
		/* Check if the node cache is full: */
		if(numCachedNodes+8>cacheSize)
			{
			/* Get the best coarsening candidate: */
			Node* coarsenNode=coarseningHeap->getTopNode();
			
			#if PARANOIA
			/* Check if this is a recipe for disaster: */
			bool disaster=false;
			for(int childIndex=0;childIndex<8&&!disaster;++childIndex)
				disaster=coarsenNode->children[childIndex].children!=0||coarsenNode->children[childIndex].subdivisionQueueIndex<subdivisionRequestQueueLength;
			if(disaster)
				std::cerr<<"Disaster!"<<std::endl;
			#endif
			
			if(coarsenNode!=0&&coarsenNode!=node->parent)
				{
				/* Check if subdividing the node will improve the overall tree: */
				if(coarsenNode->renderPass<node->renderPass||(coarsenNode->renderPass==node->renderPass&&coarsenNode->maxLOD<node->maxLOD))
					{
					/* Delete the coarsening candidate's children: */
					delete[] coarsenNode->children;
					coarsenNode->children=0;
					
					/* Update the number of cached nodes: */
					numCachedNodes-=8;
					
					/* Remove the coarsening candidate from the coarsening heap: */
					coarseningHeap->remove(coarsenNode);
					
					/* Check if the coarsening candidate's parent node is now a coarsening candidate: */
					if(coarsenNode->parent!=0)
						{
						bool canCoarsenParent=true;
						for(int childIndex=0;childIndex<8;++childIndex)
							{
							if(coarsenNode->parent->children[childIndex].children!=0)
								{
								canCoarsenParent=false;
								break;
								}
							}
						if(canCoarsenParent)
							{
							/* Insert the parent into the coarsening heap: */
							coarseningHeap->insert(coarsenNode->parent);
							}
						}
					}
				}
			}
		
		/* Check if the node can now be subdivided: */
		if(numCachedNodes+8<=cacheSize)
			{
			/* Mark the node as subdivided: */
			node->children=children;
			
			/* Remove the node's parent from the list of coarsening candidates: */
			if(node->parent!=0&&node->parent->coarseningHeapIndex!=~0x0U)
				coarseningHeap->remove(node->parent);
			
			/* Add the node to the list of coarsening candidates: */
			coarseningHeap->insert(node);
			
			/* Remove all new child nodes from the node cache: */
			for(int i=0;i<8;++i)
				children[i].pointsVersion=renderPass;
			
			/* Update the node cache: */
			numCachedNodes+=8;
			}
		else
			{
			/* Alas, all the effort was for nought -- delete the nodes just loaded: */
			delete[] children;
			}
		
		/* Unlock the node: */
		{
		Threads::Mutex::Lock loadRequestLock(loadRequestMutex);
		node->subdivisionQueueIndex=subdivisionRequestQueueLength;
		}
		}
	readyNodes.clear();
	}
	
	/* Bump up the render pass counter: */
	++renderPass;
	
	// std::cout<<numCachedNodes<<"   ";
	}

void LidarOctree::glRenderAction(const LidarOctree::Frustum& frustum,GLContextData& contextData) const
	{
	/* Retrieve the context data item: */
	DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
	
	/* Set up OpenGL state: */
	GLVertexArrayParts::enable(Vertex::getPartsMask());
	#if 0
	if(dataItem->hasPointParametersExtension)
		{
		/* Set up point parameters: */
		glPushAttrib(GL_POINT_BIT);
		GLfloat attenuation[3]={0.0f,0.0f,0.1f};
		glPointParameterfvARB(GL_POINT_DISTANCE_ATTENUATION_ARB,attenuation);
		glPointSize(10.0f);
		}
	#endif
	
	/* Render the tree: */
	dataItem->numRenderedNodes=0;
	dataItem->numCacheMisses=0;
	dataItem->numCacheBypasses=0;
	dataItem->numRenderedPoints=0;
	dataItem->numBypassedPoints=0;
	renderSubTree(&root,frustum,dataItem);
	
	/* Reset OpenGL state: */
	#if 0
	if(dataItem->hasPointParametersExtension)
		glPopAttrib();
	#endif
	if(dataItem->hasVertexBufferObjectExtension)
		glBindBufferARB(GL_ARRAY_BUFFER_ARB,0);
	GLVertexArrayParts::disable(Vertex::getPartsMask());
	
	/* Print performance counters: */
	// std::cout<<dataItem->numRenderedNodes<<", "<<dataItem->numCacheMisses<<", "<<dataItem->numCacheBypasses<<", "<<dataItem->numRenderedPoints<<", "<<dataItem->numBypassedPoints<<std::endl;
	}

Scalar LidarOctree::intersectRay(const LidarOctree::Ray& ray,Scalar coneAngle) const
	{
	/* Intersect the ray with the root domain to compute bounds on the ray parameter: */
	Vector rootSize(root.size);
	Geometry::Box<Scalar,3> rootDomain(root.center-rootSize,root.center+rootSize);
	std::pair<Scalar,Scalar> rayInterval=rootDomain.getRayParameters(ray);
	Scalar lambda1=rayInterval.first;
	if(lambda1<Scalar(0))
		lambda1=Scalar(0);
	Scalar lambda2=rayInterval.second;
	
	/* Bail out if the ray misses the domain entirely: */
	if(lambda1>=lambda2)
		return Scalar(-1);
	
	/* Recursively intersect the ray with all nodes: */
	Scalar lambda=root.intersectRay(ray,Math::sqr(coneAngle),lambda1,lambda2);
	if(lambda<lambda2)
		return lambda;
	else
		return Scalar(-1);
	}

void LidarOctree::interact(const LidarOctree::Interactor& interactor)
	{
	/* Interact with nodes recursively starting at the root: */
	interactWithSubTree(&root,interactor);
	}

void LidarOctree::selectPoints(const LidarOctree::Interactor& interactor)
	{
	/* Select points recursively starting at the root: */
	selectPoints(&root,interactor);
	}

void LidarOctree::deselectPoints(const LidarOctree::Interactor& interactor)
	{
	/* Deselect points recursively starting at the root: */
	deselectPoints(&root,interactor);
	}

void LidarOctree::clearSelection(void)
	{
	/* Clear the selection starting at the root: */
	clearSelection(&root);
	}
