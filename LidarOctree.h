/***********************************************************************
LidarOctree - Class to render multiresolution LiDAR point sets.
Copyright (c) 2005-2013 Oliver Kreylos

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

#ifndef LIDAROCTREE_INCLUDED
#define LIDAROCTREE_INCLUDED

#include <vector>
#include <Misc/HashTable.h>
#include <Threads/Mutex.h>
#include <Threads/Cond.h>
#include <Threads/Thread.h>
#include <Math/Math.h>
#include <Math/Constants.h>
#include <GL/gl.h>
#include <GL/GLColor.h>
#include <GL/GLObject.h>
#define GLGEOMETRY_NONSTANDARD_TEMPLATES
#include <GL/GLGeometryVertex.h>

#include "LidarTypes.h"
#include "Cube.h"
#include "LidarFile.h"

/* Flag whether to shift the center of the octree's domain to the coordinate system's origin: */
#define RECENTER_OCTREE 1

/* Forward declarations: */
namespace Geometry {
template <class ScalarParam,int dimensionParam>
class Ray;
}
template <class ScalarParam>
class GLFrustum;
template <class NodeParam>
class CoarseningHeap;
class PointBasedLightingShader;

class LidarOctree:public GLObject
	{
	/* Embedded classes: */
	public:
	typedef void (*TreeUpdateFunction)(void*); // Type for callback functions when image tree changes asynchronously
	
	struct Interactor // Structure describing the influence of an interaction tool on the octree
		{
		/* Elements: */
		public:
		Point center; // The interactor's center of influence
		Scalar radius; // The interactor's radius of influence
		
		/* Constructors and destructors: */
		Interactor(const Point& sCenter,const Scalar sRadius)
			:center(sCenter),radius(sRadius)
			{
			}
		};
	
	typedef std::vector<Interactor> InteractorList; // Type for lists of interactors
	typedef Geometry::Ray<Scalar,3> Ray; // Type for rays for intersection tests
	typedef GLFrustum<Scalar> Frustum; // Data type to represent view frusta
	typedef GLColor<GLubyte,4> Color; // Type for point colors
	
	struct ConeIntersection // Helper structure containing necessary state for cone intersection queries
		{
		/* Elements: */
		public:
		Ray ray; // The apex and axis direction of the intersection cone
		Scalar d2; // The squared length of the ray's direction vector
		Scalar coneAngleCos2; // The squared cosine of the intersection cone's opening angle
		Scalar testLambda1,testLambda2; // Range of valid adjusted ray parameters for the cone query
		Scalar testLambdaMin; // The smallest adjusted ray parameter already encountered
		Scalar lambdaMin; // The true ray parameter associated with the smallest adjusted ray parameter
		
		/* Constructors and destructors: */
		ConeIntersection(const Ray& sRay,Scalar sConeAngleCos) // Creates a cone query for the given axis ray and opening angle cosine
			:ray(sRay),d2(Geometry::sqr(ray.getDirection())),
			 coneAngleCos2(Math::sqr(sConeAngleCos)),
			 testLambda1(0),testLambda2(Math::Constants<Scalar>::max),
			 testLambdaMin(testLambda2),lambdaMin(testLambda2)
			{
			}
		ConeIntersection(const Ray& sRay,Scalar sConeAngleCos,Scalar sTestLambda1,Scalar sTestLambda2) // Ditto, with adjusted ray parameter limits
			:ray(sRay),d2(Geometry::sqr(ray.getDirection())),
			 coneAngleCos2(Math::sqr(sConeAngleCos)),
			 testLambda1(sTestLambda1),testLambda2(sTestLambda2),
			 testLambdaMin(testLambda2),lambdaMin(testLambda2)
			{
			}
		
		/* Methods: */
		bool isValid(void) const // Returns true if a valid intersection was found
			{
			return testLambdaMin<testLambda2;
			}
		Scalar getParameter(void) const // Returns the closest intersection ray parameter
			{
			return lambdaMin;
			}
		Point getIntersection(void) const // Returns the intersection point; assumes valid intersection was found
			{
			return ray(lambdaMin);
			}
		};
	
	private:
	typedef GLGeometry::Vertex<void,0,GLubyte,4,void,float,3> Vertex; // Type for rendered points
	typedef GLGeometry::Vertex<void,0,GLubyte,4,float,float,3> NVertex; // Type for rendered points with normal vectors
	
	struct Node // Structure for Octree nodes
		{
		/* Elements: */
		public:
		Node* parent; // Pointer to parent of this node (0 for root)
		LidarFile::Offset childrenOffset; // Offset of the node's children in the octree file (0 if node is leaf node)
		Node* children; // Pointer to an array of eight child nodes (0 if node is not subdivided)
		Cube domain; // This node's domain
		Scalar radius; // This node's radius
		unsigned int numPoints; // Number of LiDAR points belonging to this node
		bool haveNormals; // Flag whether the node's points have normal vectors associated with them
		LidarFile::Offset dataOffset; // Offset of node's data in the LiDAR data file(s)
		Scalar detailSize; // Detail size of this node, for proper LOD computation
		void* points; // Pointer to the LiDAR points belonging to this node
		unsigned int pointsVersion; // Version counter for points array to invalidate render cache on changes
		Threads::Mutex selectionMutex; // Mutex protecting the node's selection state
		bool* selectedPoints; // Array holding flags for selected points
		Vertex::Color* selectedPointColors; // Array holding the original colors of selected points
		
		/* State touched during rendering traversals: */
		mutable Threads::Mutex nodeMutex; // Mutex protecting node state touched during rendering traversal
		mutable unsigned int renderPass; // Counter value of last rendering pass that entered this node
		mutable Scalar maxLOD; // Maximum LOD value of node during current rendering pass
		mutable unsigned int subdivisionQueueIndex; // Node's index in the node loader's request queue, == queue size if no request pending, == ~0x0U if locked by node loader
		mutable unsigned int coarseningHeapIndex; // Node's index in the heap of coarsening candidates; == ~0x0U if not a coarsening candidate
		
		/* Constructors and destructors: */
		Node(void) // Creates a leaf node without points
			:parent(0),children(0),
			 points(0),pointsVersion(0),
			 selectedPoints(0),selectedPointColors(0),
			 renderPass(0U),
			 subdivisionQueueIndex(~0x0U-1),
			 coarseningHeapIndex(~0x0U)
			{
			}
		~Node(void); // Destroys a node and its subtree
		
		/* Methods: */
		void intersectCone(ConeIntersection& cone) const; // Recursively intersects a cone with this node's subtree
		};
	
	struct SubdivisionRequest // Structure to store subdivision request for nodes
		{
		/* Elements: */
		public:
		Node* node; // Node to be subdivided
		Scalar LOD; // LOD value of the node in the most recent rendering pass
		
		/* Constructors and destructors: */
		SubdivisionRequest(void)
			:node(0),LOD(0.0f)
			{
			}
		};
	
	struct DataItem:public GLObject::DataItem
		{
		/* Embedded classes: */
		public:
		struct CacheSlot // Structure for a slot in the node cache
			{
			/* Elements: */
			public:
			GLuint vertexBufferObjectId; // ID of vertex buffer object for this cache slot (only used if vertex array objects supported)
			const Node* node; // Pointer to node currently cached in this slot (0 if empty)
			unsigned int version; // Version of point data currently stored in cache
			unsigned int lastUsed; // Counter value of last rendering pass in which this cache slot was used (to prohibit cache thrashing)
			CacheSlot* pred; // Pointer to predecessor in LRU cache slot list
			CacheSlot* succ; // Pointer to successor in LRU cache slot list
			
			/* Constructors and destructors: */
			CacheSlot(void) // Creates empty cache slot
				:vertexBufferObjectId(0),node(0),version(0),
				lastUsed(0),pred(0),succ(0)
				{
				}
			};
		
		typedef Misc::HashTable<const Node*,CacheSlot*> NodeHasher; // Hash table to map from nodes to the cache slots containing their points
		
		/* Elements: */
		bool hasVertexBufferObjectExtension; // Flag if vertex buffer objects are supported
		unsigned int cacheSize; // Number of node vertex arrays to keep in graphics card memory
		CacheSlot* cacheSlots; // Array of cache slots
		NodeHasher cacheNodeMap; // Hash table of nodes whose points are currently held in graphics card memory
		CacheSlot* lruHead; // Pointer to first (least recently used) cache slot in LRU list
		CacheSlot* lruTail; // Pointer to last cache slot in LRU list
		GLuint bypassVertexBufferObjectId; // ID of the vertex buffer object used to render nodes bypassing the cache
		bool hasPointParametersExtension; // Flag if point parameters are supported
		unsigned int numRenderedNodes; // Total number of nodes rendered in the last render pass
		unsigned int numCacheMisses; // Total number of cache misses in the last render pass
		unsigned int numCacheBypasses; // Total number of cache bypasses in the last render pass
		unsigned int numRenderedPoints; // Total number of points rendered in the last render pass
		unsigned int numBypassedPoints; // Total number of points rendered without using the cache
		
		/* Constructors and destructors: */
		DataItem(unsigned int sCacheSize);
		virtual ~DataItem(void);
		};
	
	/* Elements: */
	LidarFile indexFile; // The file containing the octree's structural data
	LidarFile pointsFile; // The file containing the octree's point data
	LidarFile::Offset pointsRecordSize; // Record size of points file
	LidarFile* normalsFile; // Pointer to an optional file containing normal vectors for each point
	LidarFile::Offset normalsRecordSize; // Record size of normals file
	LidarFile* colorsFile; // Pointer to an optional file containing colors for each point
	LidarFile::Offset colorsRecordSize; // Record size of colors file
	unsigned int maxNumPointsPerNode; // Maximum number of points per node
	Vector pointOffset; // Offset from original LiDAR point positions to re-centered point positions
	Node root; // The octree's root node (offset to center around the origin)
	Scalar maxRenderLOD; // Maximum LOD value at which a node will be rendered
	Point fncCenter; // Center point of focus region
	Scalar fncRadius; // Radius of focus region
	Scalar fncWeight; // Weight for focus+context LOD adjustment
	float baseSurfelSize; // Base splat size for leaf nodes
	float surfelScale; // Scale factor applied to adjusted splat size
	unsigned int cacheSize; // Maximum number of nodes to hold in memory at any time
	unsigned int glCacheSize; // Maximum number of nodes to hold in graphics card memory at any time
	unsigned int numCachedNodes; // Number of nodes currently in the memory cache
	unsigned int renderPass; // Render pass counter
	
	/* Elements for communication with the node loading thread: */
	Threads::Mutex treeUpdateFunctionMutex; // Mutex protecting the tree update function's state
	TreeUpdateFunction treeUpdateFunction; // Callback function called whenever the tree is updated asynchronously
	void* treeUpdateFunctionArg; // Arbitrary argument for tree update function
	Threads::Thread nodeLoaderThread; // Handle for the node loading thread
	mutable Threads::Mutex loadRequestMutex; // Mutex protecting the tile loader thread state
	mutable Threads::Cond loadRequestCond; // Condition variable to signal a new tile loading request
	mutable unsigned int subdivisionRequestQueueLength; // Maximum number of pending subdivision requests
	mutable SubdivisionRequest* subdivisionRequestQueue; // Queue of subdivision requests to the node loader thread
	Threads::Mutex readyNodesMutex; // Mutex protecting the list of ready nodes
	std::vector<Node*> readyNodes; // List of ready child nodes
	mutable Threads::Mutex coarseningHeapMutex; // Mutex protecting the coarsening heap during rendering
	mutable CoarseningHeap<Node>* coarseningHeap; // Heap of nodes that are candidates for coarsening
	
	/* Private methods: */
	void renderSubTree(const Node* node,const Frustum& frustum,PointBasedLightingShader& pbls,DataItem* dataItem) const;
	void interactWithSubTree(Node* node,const Interactor& interactor); // Prepares a subtree for interaction with an interactor
	template <class VertexParam>
	bool selectPoint(Node* node,unsigned int pointIndex); // Selects the given point in the given node; returns true if selection changed
	template <class VertexParam>
	void selectPointsInNode(Node* node,const Interactor& interactor); // Selects points in the given node
	void selectPoints(Node* node,const Interactor& interactor); // Selects points in the given subtree
	template <class VertexParam>
	void deselectPointsInNode(Node* node,const Interactor& interactor); // Deselects points in the given node
	bool deselectPoints(Node* node,const Interactor& interactor); // Deselects points in the given subtree; returns true if node is removable (i.e., has no children or selected points)
	bool clearSelection(Node* node); // Clears selected points in the given subtree; returns true if node is removable (i.e., has no children or selected points)
	template <class VertexParam,class PointProcessorParam>
	void processPointsInBox(const Box& box,const Node* node,PointProcessorParam& dpp) const; // Processes points inside the given box in the given subtree
	template <class VertexParam,class DirectedPointProcessorParam>
	void processPointsDirected(const Node* node,DirectedPointProcessorParam& dpp) const; // Processes points in the given subtree in order of increasing distance from the processor's query point
	template <class VertexParam,class PointProcessorParam>
	void processSelectedPoints(const Node* node,PointProcessorParam& pp) const; // Processes selected points in the given subtree
	template <class PointNormalProcessorParam>
	void processSelectedPointsWithNormals(const Node* node,PointNormalProcessorParam& pp) const; // Processes selected points in the given subtree
	template <class PointNormalProcessorParam>
	void processSelectedPointsWithNullNormals(const Node* node,PointNormalProcessorParam& pp) const; // Processes selected points in the given subtree
	template <class VertexParam,class ColoringPointProcessorParam>
	void colorSelectedPoints(Node* node,ColoringPointProcessorParam& cpp); // Colors selected points in the given subtree
	void loadNodePoints(Node* node); // Loads the point array of the given node
	template <class VertexParam>
	void selectCloseNeighbors(Node* node,unsigned int left,unsigned int right,int splitDimension,const VertexParam& point,Scalar maxDist);
	template <class VertexParam>
	void propagateSelectedPoints(Node* node,Node* children); // Propagates selected points from the given node into its just loaded child nodes
	void* nodeLoaderThreadMethod(void); // Thread method to subdivide octree nodes without interrupting rendering
	
	/* Constructors and destructors: */
	public:
	LidarOctree(const char* lidarFileName,size_t sCacheSize,size_t sGlCacheSize); // Creates octree from the given LiDAR file; cache sizes are in bytes
	virtual ~LidarOctree(void);
	
	/* Methods: */
	virtual void initContext(GLContextData& contextData) const;
	bool hasNormalVectors(void) const // Returns true if the octree has normal vectors for each point
		{
		return root.haveNormals;
		}
	Point getDomainCenter(void) const; // Returns the center of the octree's domain
	Scalar getDomainRadius(void) const; // Returns the radius of the octree's domain
	const Vector& getPointOffset(void) const // Returns an offset vector to transform LiDAR octree coordinates to source point coordinates
		{
		return pointOffset;
		};
	void setTreeUpdateFunction(TreeUpdateFunction newTreeUpdateFunction,void* newTreeUpdateFunctionArg);
	void setRenderQuality(Scalar qualityLevel); // Sets the quality level for rendering, 0 is normal quality
	void setFocusAndContext(const Point& newFncCenter,Scalar newFncRadius,Scalar newFncWeight); // Adjusts focus+context LOD adjustment parameters
	void setBaseSurfelSize(float newBaseSurfelSize,float newSurfelScale); // Sets the splat size for leaf nodes
	void startRenderPass(void); // Starts the next rendering pass of the point octree
	void glRenderAction(const Frustum& frustum,PointBasedLightingShader& pbls,GLContextData& contextData) const;
	void intersectCone(ConeIntersection& cone) const; // Intersects a cone with all points in the current octree
	void interact(const Interactor& interactor); // Prepares an octree for interaction with an interactor
	void selectPoints(const Interactor& interactor); // Selects all points inside the interactor's region of influence
	void deselectPoints(const Interactor& interactor); // Deselects all points inside the interactor's region of influence
	void clearSelection(void); // Clears the set of selected points
	template <class PointProcessorParam>
	void processPointsInBox(const Box& box,PointProcessorParam& dpp) const // Processes points in leaf nodes that are inside the given box
		{
		/* Process nodes from the root: */
		if(root.haveNormals)
			processPointsInBox<NVertex,PointProcessorParam>(box,&root,dpp);
		else
			processPointsInBox<Vertex,PointProcessorParam>(box,&root,dpp);
		}
	template <class DirectedPointProcessorParam>
	void processPointsDirected(DirectedPointProcessorParam& dpp) const // Processes points in leaf nodes in order of approximately increasing distance from the processor's query point
		{
		/* Process nodes from the root: */
		if(root.haveNormals)
			processPointsDirected<NVertex,DirectedPointProcessorParam>(&root,dpp);
		else
			processPointsDirected<Vertex,DirectedPointProcessorParam>(&root,dpp);
		}
	template <class PointProcessorParam>
	void processSelectedPoints(PointProcessorParam& pp) const // Processes the set of selected points in leaf nodes with the given point processor
		{
		/* Process nodes from the root: */
		if(root.haveNormals)
			processSelectedPoints<NVertex,PointProcessorParam>(&root,pp);
		else
			processSelectedPoints<Vertex,PointProcessorParam>(&root,pp);
		};
	template <class PointNormalProcessorParam>
	void processSelectedPointsWithNormals(PointNormalProcessorParam& pp) const // Processes the set of selected points in leaf nodes with the given point processor
		{
		/* Process nodes from the root: */
		if(root.haveNormals)
			processSelectedPointsWithNormals<PointNormalProcessorParam>(&root,pp);
		else
			processSelectedPointsWithNullNormals<PointNormalProcessorParam>(&root,pp);
		};
	template <class ColoringPointProcessorParam>
	void colorSelectedPoints(ColoringPointProcessorParam& cpp) // Ditto, but allows point processor to change the selection color of selected points
		{
		/* Process nodes from the root: */
		if(root.haveNormals)
			colorSelectedPoints<NVertex,ColoringPointProcessorParam>(&root,cpp);
		else
			colorSelectedPoints<Vertex,ColoringPointProcessorParam>(&root,cpp);
		};
	};

#ifndef LIDAROCTREE_IMPLEMENTATION
#include "LidarOctree.icpp"
#endif

#endif
