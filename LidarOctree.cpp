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

LidarOctree::Node::~Node(void)
	{
	/* Delete point array: */
	if(haveNormals)
		delete[] static_cast<NVertex*>(points);
	else
		delete[] static_cast<Vertex*>(points);
	
	/* Delete selected point flag mask: */
	delete[] selectedPoints;
	delete[] selectedPointColors;
	
	/* Delete children: */
	delete[] children;
	}

namespace {

/***************
Helper function:
***************/

template <class VertexParam>
inline
Scalar
intersectRayWithPoints(
	const LidarOctree::Ray& ray,
	Scalar lambdaMin,
	Scalar coneAngle2,
	const VertexParam* points,
	unsigned int numPoints)
	{
	/* Intersect the ray with all points in this node: */
	for(unsigned int i=0;i<numPoints;++i)
		{
		/* Check if the point is closer than the previous one and inside the cone: */
		Vector sp=points[i].position-ray.getOrigin();
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

Scalar LidarOctree::Node::intersectRay(const LidarOctree::Ray& ray,Scalar coneAngle2,Scalar lambda1,Scalar lambda2) const
	{
	if(children!=0)
		{
		#if 1
		Scalar lambdaMin=lambda2;
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			/* Convert the child's domain into a box: */
			Geometry::Box<Scalar,3> childDomainBox(children[childIndex].domain.getMin(),children[childIndex].domain.getMax());
			
			/* Intersect the ray with the child's domain: */
			std::pair<Scalar,Scalar> lambdas=childDomainBox.getRayParameters(ray);
			if(lambdas.first<lambda1)
				lambdas.first=lambda1;
			if(lambdas.second>lambda2)
				lambdas.second=lambda2;
			
			/* Recurse into the child if the ray intersects its domain: */
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
			Scalar center=domain.getCenter(i);
			if(ray.getOrigin()[i]>=center)
				{
				childIndex|=0x1<<i;
				if(ray.getDirection()[i]<Scalar(0))
					planeLambdas[i]=(center-ray.getOrigin()[i])/ray.getDirection()[i];
				else
					planeLambdas[i]=Scalar(-1);
				}
			else
				{
				if(ray.getDirection()[i]>Scalar(0))
					planeLambdas[i]=(center-ray.getOrigin()[i])/ray.getDirection()[i];
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
		if(haveNormals)
			return intersectRayWithPoints(ray,lambda2,coneAngle2,static_cast<const NVertex*>(points),numPoints);
		else
			return intersectRayWithPoints(ray,lambda2,coneAngle2,static_cast<const Vertex*>(points),numPoints);
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
		Point p;
		for(int i=0;i<3;++i)
			p[i]=normal[i]>Scalar(0)?node->domain.getMax()[i]:node->domain.getMin()[i];
		
		/* Check if the point is inside the view frustum: */
		inside=!frustum.getFrustumPlane(plane).contains(p);
		}
	
	/* Bail out if this node is not inside the view frustum: */
	if(!inside)
		return;
	
	/* Calculate the node's projected detail size: */
	Point nodeCenter=node->domain.getCenter();
	Scalar projectedDetailSize=frustum.calcProjectedRadius(nodeCenter,node->detailSize);
	if(projectedDetailSize>=Scalar(0))
		{
		/* Adjust the projected detail size based on the focus+context rule: */
		if(fncWeight>Scalar(0))
			{
			/* Calculate the distance from the node to the focus region: */
			Scalar fncDist=Geometry::dist(nodeCenter,fncCenter)-node->radius;
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
			/* Use the view point for a view-ordered traversal of the child nodes: */
			int childIndex=0x0;
			Point eye=frustum.getEye().toPoint();
			for(int i=0;i<3;++i)
				if(eye[i]>=nodeCenter[i])
					childIndex|=0x1<<i;
			
			/* Render the node's children: */
			// for(int i=7;i>=0;--i) // Back-to-front rendering
			for(int i=0;i<8;++i) // Front-to-back rendering
				renderSubTree(&node->children[i^childIndex],frustum,dataItem);
			
			/* Done here... */
			return;
			}
		else if(node->childrenOffset!=LidarFile::Offset(0))
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
	
	/* Retrieve the cache slot containing this node's vertex array: */
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
				{
				size_t arraySize=node->numPoints*(node->haveNormals?sizeof(NVertex):sizeof(Vertex));
				glBufferDataARB(GL_ARRAY_BUFFER_ARB,arraySize,node->points,GL_DYNAMIC_DRAW_ARB);
				}
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
			size_t arraySize=node->numPoints*(node->haveNormals?sizeof(NVertex):sizeof(Vertex));
			glBufferDataARB(GL_ARRAY_BUFFER_ARB,arraySize,node->points,GL_STREAM_DRAW_ARB);
			}
		dataItem->numBypassedPoints+=node->numPoints;
		}
	
	/* Render this node's point set: */
	if(vertexBufferObjectId!=0)
		{
		/* Render from the vertex buffer: */
		if(node->haveNormals)
			glVertexPointer(static_cast<const NVertex*>(0));
		else
			glVertexPointer(static_cast<const Vertex*>(0));
		}
	else
		{
		/* Render straight from the node vertex array (only happens if GL_ARB_vertex_buffer_object not supported): */
		if(node->haveNormals)
			glVertexPointer(static_cast<const NVertex*>(node->points));
		else
			glVertexPointer(static_cast<const Vertex*>(node->points));
		}
	glDrawArrays(GL_POINTS,0,node->numPoints);
	
	#if 0
	glColor3f(1.0f,0.0f,0.0f);
	glBegin(GL_LINE_STRIP);
	glVertex3f(node->domain.getMin()[0],node->domain.getMin()[1],node->domain.getMin()[2]);
	glVertex3f(node->domain.getMax()[0],node->domain.getMin()[1],node->domain.getMin()[2]);
	glVertex3f(node->domain.getMax()[0],node->domain.getMax()[1],node->domain.getMin()[2]);
	glVertex3f(node->domain.getMin()[0],node->domain.getMax()[1],node->domain.getMin()[2]);
	glVertex3f(node->domain.getMin()[0],node->domain.getMin()[1],node->domain.getMin()[2]);
	glVertex3f(node->domain.getMin()[0],node->domain.getMin()[1],node->domain.getMax()[2]);
	glVertex3f(node->domain.getMax()[0],node->domain.getMin()[1],node->domain.getMax()[2]);
	glVertex3f(node->domain.getMax()[0],node->domain.getMax()[1],node->domain.getMax()[2]);
	glVertex3f(node->domain.getMin()[0],node->domain.getMax()[1],node->domain.getMax()[2]);
	glVertex3f(node->domain.getMin()[0],node->domain.getMin()[1],node->domain.getMax()[2]);
	glEnd();
	glBegin(GL_LINES);
	glVertex3f(node->domain.getMax()[0],node->domain.getMin()[1],node->domain.getMin()[2]);
	glVertex3f(node->domain.getMax()[0],node->domain.getMin()[1],node->domain.getMax()[2]);
	glVertex3f(node->domain.getMax()[0],node->domain.getMax()[1],node->domain.getMin()[2]);
	glVertex3f(node->domain.getMax()[0],node->domain.getMax()[1],node->domain.getMax()[2]);
	glVertex3f(node->domain.getMin()[0],node->domain.getMax()[1],node->domain.getMin()[2]);
	glVertex3f(node->domain.getMin()[0],node->domain.getMax()[1],node->domain.getMax()[2]);
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
	Scalar interactorDist2=node->domain.sqrDist(interactor.center);
	if(interactorDist2>=Math::sqr(interactor.radius))
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
		if(node->childrenOffset!=LidarFile::Offset(0))
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
	Scalar dist2=node->domain.sqrDist(interactor.center);
	Scalar ir2=Math::sqr(interactor.radius);
	if(dist2<ir2)
		{
		{
		Threads::Mutex::Lock selectionLock(node->selectionMutex);
		
		bool selectionEmpty=node->selectedPoints==0;
		
		/* Select points in this node: */
		if(node->haveNormals)
			selectPointsInNode<NVertex>(node,interactor);
		else
			selectPointsInNode<Vertex>(node,interactor);
		
		/* Check if this node just had its first points selected: */
		if(selectionEmpty&&node->selectedPoints!=0)
			{
			/* Check if the node's parent is in the coarsening heap: */
			if(node->parent!=0&&node->parent->coarseningHeapIndex!=~0x0)
				{
				Threads::Mutex::Lock coarseningHeapLock(coarseningHeapMutex);
				coarseningHeap->remove(node->parent);
				}
			}
		}
		
		if(node->children!=0)
			{
			/* Recurse into the node's children: */
			for(int childIndex=0;childIndex<8;++childIndex)
				selectPoints(&node->children[childIndex],interactor);
			}
		}
	}

bool LidarOctree::deselectPoints(LidarOctree::Node* node,const LidarOctree::Interactor& interactor)
	{
	/* Check if the interactor's region of influence intersects the node's domain: */
	Scalar dist2=node->domain.sqrDist(interactor.center);
	Scalar ir2=Math::sqr(interactor.radius);
	if(dist2<ir2)
		{
		if(node->selectedPoints!=0)
			{
			Threads::Mutex::Lock selectionLock(node->selectionMutex);
			
			/* Deselect points in this node: */
			if(node->haveNormals)
				deselectPointsInNode<NVertex>(node,interactor);
			else
				deselectPointsInNode<Vertex>(node,interactor);
			}
		
		if(node->children!=0)
			{
			/* Recurse into the node's children: */
			bool canCoarsen=true;
			for(int childIndex=0;childIndex<8;++childIndex)
				canCoarsen=deselectPoints(&node->children[childIndex],interactor)&&canCoarsen;
			
			if(canCoarsen&&node->coarseningHeapIndex==~0x0)
				{
				/* Insert the node into the coarsening heap: */
				Threads::Mutex::Lock coarseningHeapLock(coarseningHeapMutex);
				coarseningHeap->insert(node);
				}
			}
		}
	
	return node->children==0&&node->selectedPoints==0;
	}

bool LidarOctree::clearSelection(LidarOctree::Node* node)
	{
	if(node->selectedPoints!=0)
		{
		Threads::Mutex::Lock selectionLock(node->selectionMutex);
		
		/* Deselect all selected points: */
		bool pointsChanged=false;
		if(node->haveNormals)
			{
			NVertex* points=static_cast<NVertex*>(node->points);
			for(unsigned int i=0;i<node->numPoints;++i)
				{
				if(node->selectedPoints[i])
					{
					GLubyte intensity=node->selectedPoints[i];
					points[i].color=node->selectedPointColors[i];
					pointsChanged=true;
					}
				}
			}
		else
			{
			Vertex* points=static_cast<Vertex*>(node->points);
			for(unsigned int i=0;i<node->numPoints;++i)
				{
				if(node->selectedPoints[i])
					{
					GLubyte intensity=node->selectedPoints[i];
					points[i].color=node->selectedPointColors[i];
					pointsChanged=true;
					}
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
		bool canCoarsen=true;
		for(int childIndex=0;childIndex<8;++childIndex)
			canCoarsen=clearSelection(&node->children[childIndex])&&canCoarsen;
		
		if(canCoarsen&&node->coarseningHeapIndex==~0x0)
			{
			/* Insert the node into the coarsening heap: */
			Threads::Mutex::Lock coarseningHeapLock(coarseningHeapMutex);
			coarseningHeap->insert(node);
			}
		}
	
	return node->children==0;
	}

void LidarOctree::loadNodePoints(LidarOctree::Node* node)
	{
	if(node->haveNormals)
		{
		/* Create the node's point array (always allocate the maximum size to prevent memory fragmentation): */
		NVertex* points=new NVertex[maxNumPointsPerNode];
		node->points=points;
		
		/* Load the node's points into a temporary point buffer: */
		LidarPoint* pointsBuffer=new LidarPoint[maxNumPointsPerNode];
		pointsFile.seekSet(LidarDataFileHeader::getFileSize()+pointsRecordSize*node->dataOffset);
		pointsFile.read(pointsBuffer,node->numPoints);
		Vector* normalsBuffer=new Vector[maxNumPointsPerNode];
		normalsFile->seekSet(LidarDataFileHeader::getFileSize()+normalsRecordSize*node->dataOffset);
		normalsFile->read(normalsBuffer,node->numPoints);
		
		/* Convert the LiDAR points to render points: */
		for(unsigned int i=0;i<node->numPoints;++i)
			{
			/* Copy the point color: */
			for(int j=0;j<4;++j)
				points[i].color[j]=pointsBuffer[i].value[j];
			
			/* Copy the normal vector: */
			points[i].normal=normalsBuffer[i];
			
			/* Copy the point position: */
			points[i].position=pointsBuffer[i];
			#if RECENTER_OCTREE
			/* Offset the points so that the root node's center is the origin: */
			points[i].position-=pointOffset;
			#endif
			}
		
		/* Clean up: */
		delete[] pointsBuffer;
		delete[] normalsBuffer;
		}
	else
		{
		/* Create the node's point array (always allocate the maximum size to prevent memory fragmentation): */
		Vertex* points=new Vertex[maxNumPointsPerNode];
		node->points=points;
		
		/* Load the node's points into a temporary point buffer: */
		LidarPoint* pointsBuffer=new LidarPoint[maxNumPointsPerNode];
		pointsFile.seekSet(LidarDataFileHeader::getFileSize()+pointsRecordSize*node->dataOffset);
		pointsFile.read(pointsBuffer,node->numPoints);
		
		/* Convert the LiDAR points to render points: */
		for(unsigned int i=0;i<node->numPoints;++i)
			{
			/* Copy the point color: */
			for(int j=0;j<4;++j)
				points[i].color[j]=pointsBuffer[i].value[j];
			
			/* Copy the point position: */
			points[i].position=pointsBuffer[i];
			#if RECENTER_OCTREE
			/* Offset the points so that the root node's center is the origin: */
			points[i].position-=pointOffset;
			#endif
			}
		
		/* Clean up: */
		delete[] pointsBuffer;
		}
	}

template <class VertexParam>
inline
void
LidarOctree::selectCloseNeighbors(
	LidarOctree::Node* node,
	unsigned int left,
	unsigned int right,
	int splitDimension,
	const VertexParam& point,
	Scalar maxDist)
	{
	/* Calculate the index of the current point: */
	unsigned int mid=(left+right)>>1;
	VertexParam& np=static_cast<VertexParam*>(node->points)[mid];
	
	int childSplitDimension=splitDimension+1;
	if(childSplitDimension==3)
		childSplitDimension=0;
	
	/* Traverse into child closer to query point: */
	if(point.position[splitDimension]<np.position[splitDimension])
		{
		/* Traverse left child: */
		if(left<mid)
			selectCloseNeighbors(node,left,mid-1,childSplitDimension,point,maxDist);
		
		/* Process the current point: */
		if(Geometry::sqrDist(point.position,np.position)<=Math::sqr(maxDist))
			{
			/* Select the point: */
			selectPoint<VertexParam>(node,mid);
			}
		
		/* Traverse the right child: */
		if(point.position[splitDimension]+maxDist>=np.position[splitDimension]&&right>mid)
			selectCloseNeighbors(node,mid+1,right,childSplitDimension,point,maxDist);
		}
	else
		{
		/* Traverse right child: */
		if(right>mid)
			selectCloseNeighbors(node,mid+1,right,childSplitDimension,point,maxDist);
		
		/* Process the current point: */
		if(Geometry::sqrDist(point.position,np.position)<=Math::sqr(maxDist))
			{
			/* Select the point: */
			selectPoint<VertexParam>(node,mid);
			}
		
		/* Traverse the left child: */
		if(point.position[splitDimension]-maxDist<=np.position[splitDimension]&&left<mid)
			selectCloseNeighbors(node,left,mid-1,childSplitDimension,point,maxDist);
		}
	}

template <class VertexParam>
inline
void
LidarOctree::propagateSelectedPoints(
	LidarOctree::Node* node,
	LidarOctree::Node* children)
	{
	VertexParam* nodePoints=static_cast<VertexParam*>(node->points);
	for(int childIndex=0;childIndex<8;++childIndex)
		if(children[childIndex].numPoints!=0)
			{
			VertexParam* childPoints=static_cast<VertexParam*>(children[childIndex].points);
			
			/* Process each selected point in the node: */
			for(unsigned int i=0;i<node->numPoints;++i)
				if(node->selectedPoints[i])
					{
					/* Select the point's close neighbors in the child node: */
					selectCloseNeighbors<VertexParam>(&children[childIndex],0,children[childIndex].numPoints-1,0,nodePoints[i],node->detailSize);
					}
			}
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
		indexFile.seekSet(node->childrenOffset);
		for(int childIndex=0;childIndex<8;++childIndex)
			{
			LidarOctreeFileNode ofn;
			ofn.read(indexFile);
			children[childIndex].parent=node;
			children[childIndex].childrenOffset=ofn.childrenOffset;
			children[childIndex].domain=Cube(node->domain,childIndex);
			children[childIndex].radius=Math::div2(node->radius);
			children[childIndex].numPoints=ofn.numPoints;
			children[childIndex].haveNormals=node->haveNormals;
			children[childIndex].dataOffset=ofn.dataOffset;
			children[childIndex].detailSize=ofn.detailSize;
			children[childIndex].subdivisionQueueIndex=subdivisionRequestQueueLength;
			}
		
		/* Load the node's children's point data: */
		for(int childIndex=0;childIndex<8;++childIndex)
			loadNodePoints(&children[childIndex]);
		
		{
		Threads::Mutex::Lock nodeSelectionLock(node->selectionMutex);
		if(node->selectedPoints!=0)
			{
			/* Propagate selected points from the node to its children: */
			if(node->haveNormals)
				propagateSelectedPoints<NVertex>(node,children);
			else
				propagateSelectedPoints<Vertex>(node,children);
			}
		}
		
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

LidarOctree::LidarOctree(const char* lidarFileName,unsigned int sCacheSize,unsigned int sGlCacheSize)
	:indexFile(getLidarPartFileName(lidarFileName,"Index").c_str(),"rb",LidarFile::LittleEndian),
	 pointsFile(getLidarPartFileName(lidarFileName,"Points").c_str(),"rb",LidarFile::LittleEndian),
	 normalsFile(0),
	 maxRenderLOD(Math::sqrt(Scalar(2))),
	 fncWeight(0),
	 numCachedNodes(0),
	 renderPass(0U),
	 treeUpdateFunction(0),treeUpdateFunctionArg(0),
	 subdivisionRequestQueueLength(4),subdivisionRequestQueue(new SubdivisionRequest[subdivisionRequestQueueLength]),
	 coarseningHeap(0)
	{
	/* Check if there is a normal vector file: */
	try
		{
		normalsFile=new LidarFile(getLidarPartFileName(lidarFileName,"Normals").c_str(),"rb",LidarFile::LittleEndian);
		}
	catch(std::runtime_error err)
		{
		}
	
	/* Read the octree file header: */
	LidarOctreeFileHeader ofh(indexFile);
	
	/* Initialize the root node's domain: */
	#if RECENTER_OCTREE
	pointOffset=ofh.domain.getCenter()-Point::origin;
	root.domain=Cube(ofh.domain.getMin()-pointOffset,ofh.domain.getMax()-pointOffset);
	#else
	pointOffset=Vector::zero;
	root.domain=ofh.domain;
	#endif
	Scalar rootRadius2=Scalar(0);
	for(int i=0;i<3;++i)
		rootRadius2+=Math::sqr(root.domain.getMax()[i]-root.domain.getMin()[i]);
	root.radius=Math::div2(Math::sqrt(rootRadius2));
	
	/* Initialize the tree structure: */
	maxNumPointsPerNode=ofh.maxNumPointsPerNode;
	root.haveNormals=normalsFile!=0;
	
	/* Calculate the memory and GPU cache sizes: */
	size_t vertexSize=root.haveNormals?sizeof(NVertex):sizeof(Vertex);
	size_t memNodeSize=sizeof(Node)+maxNumPointsPerNode*vertexSize;
	size_t glNodeSize=maxNumPointsPerNode*vertexSize;
	cacheSize=sCacheSize/memNodeSize;
	if(cacheSize==0)
		Misc::throwStdErr("LidarOctree::LidarOctree: Specified memory cache size too small");
	glCacheSize=sGlCacheSize/glNodeSize;
	if(glCacheSize==0)
		Misc::throwStdErr("LidarOctree::LidarOctree: Specified GPU cache size too small");
	std::cout<<"Cache sizes: "<<cacheSize<<" memory nodes, "<<glCacheSize<<" GPU nodes"<<std::endl;
	
	/* Read the root node's structure: */
	LidarOctreeFileNode rootfn;
	rootfn.read(indexFile);
	root.childrenOffset=rootfn.childrenOffset;
	root.numPoints=rootfn.numPoints;
	root.dataOffset=rootfn.dataOffset;
	root.detailSize=rootfn.detailSize;
	root.subdivisionQueueIndex=subdivisionRequestQueueLength;
	
	/* Read the point file's header: */
	LidarDataFileHeader dfh(pointsFile);
	pointsRecordSize=LidarFile::Offset(dfh.recordSize);
	
	if(root.haveNormals)
		{
		/* Read the normal file's header: */
		LidarDataFileHeader dfh(*normalsFile);
		normalsRecordSize=LidarFile::Offset(dfh.recordSize);
		}
	
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
	
	/* Delete the optional normals file: */
	delete normalsFile;
	}

void LidarOctree::initContext(GLContextData& contextData) const
	{
	/* Create a context data item: */
	DataItem* dataItem=new DataItem(glCacheSize);
	contextData.addDataItem(this,dataItem);
	}

Point LidarOctree::getDomainCenter(void) const
	{
	return root.domain.getCenter();
	}

Scalar LidarOctree::getDomainRadius(void) const
	{
	#if 1
	return root.radius;
	#else
	Scalar radius2=Scalar(0);
	for(int i=0;i<3;++i)
		radius2+=Math::sqr(root.domain.getMax()[i]-root.domain.getMin()[i]);
	return Math::div2(Math::sqrt(radius2));
	#endif
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
			Threads::Mutex::Lock coarseningHeapLock(coarseningHeapMutex);
			// DEBUGGING
			//std::cout<<coarseningHeap->getNumItems()<<std::endl;
			Node* coarsenNode=coarseningHeap->getTopNode();
			
			if(coarsenNode!=0&&coarsenNode!=node->parent)
				{
				/* Check if subdividing the node will improve the overall tree: */
				if(coarsenNode->renderPass<node->renderPass||(coarsenNode->renderPass==node->renderPass&&coarsenNode->maxLOD<node->maxLOD))
					{
					// DEBUGGING
					/* Check if something bad happened: */
					bool isRemovable=true;
					for(int i=0;i<8&&isRemovable;++i)
						isRemovable=coarsenNode->children[i].children==0&&coarsenNode->children[i].selectedPoints==0;
					if(!isRemovable)
						std::cout<<"Removing non-removable nodes!"<<std::endl;
					
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
							Node& sibling=coarsenNode->parent->children[childIndex];
							if(sibling.children!=0||sibling.selectedPoints!=0)
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
			Threads::Mutex::Lock coarseningHeapLock(coarseningHeapMutex);
			
			/* Mark the node as subdivided: */
			node->children=children;
			
			/* Remove the node's parent from the list of coarsening candidates: */
			if(node->parent!=0&&node->parent->coarseningHeapIndex!=~0x0U)
				coarseningHeap->remove(node->parent);
			
			/* Check if the just refined node can be coarsened again: */
			bool canCoarsen=true;
			for(int i=0;i<8&&canCoarsen;++i)
				canCoarsen=node->children[i].selectedPoints==0;
			if(canCoarsen)
				{
				/* Add the node to the list of coarsening candidates: */
				coarseningHeap->insert(node);
				}
			
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
	if(root.haveNormals)
		GLVertexArrayParts::enable(NVertex::getPartsMask());
	else
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
	if(root.haveNormals)
		GLVertexArrayParts::disable(NVertex::getPartsMask());
	else
		GLVertexArrayParts::disable(Vertex::getPartsMask());
	
	/* Print performance counters: */
	// std::cout<<dataItem->numRenderedNodes<<", "<<dataItem->numCacheMisses<<", "<<dataItem->numCacheBypasses<<", "<<dataItem->numRenderedPoints<<", "<<dataItem->numBypassedPoints<<std::endl;
	}

Scalar LidarOctree::intersectRay(const LidarOctree::Ray& ray,Scalar coneAngle) const
	{
	/* Convert the root's domain into a box: */
	Geometry::Box<Scalar,3> rootDomainBox(root.domain.getMin(),root.domain.getMax());
	
	/* Intersect the ray with the root's domain: */
	std::pair<Scalar,Scalar> lambdas=rootDomainBox.getRayParameters(ray);
	if(lambdas.first<Scalar(0))
		lambdas.first=Scalar(0);
	
	/* Bail out if the ray misses the domain entirely: */
	if(lambdas.first>=lambdas.second)
		return Scalar(-1);
	
	/* Recursively intersect the ray with all nodes: */
	Scalar lambda=root.intersectRay(ray,Math::sqr(coneAngle),lambdas.first,lambdas.second);
	if(lambda<lambdas.second)
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
