/***********************************************************************
CoarseningHeap - Helper class to store terrain tree nodes that can be
removed from the node cache, in order of least recent/most finely
resolved first.
Copyright (c) 2007-2008 Oliver Kreylos

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

#ifndef COARSENINGHEAP_INCLUDED
#define COARSENINGHEAP_INCLUDED

template <class NodeParam>
class CoarseningHeap
	{
	/* Embedded classes: */
	private:
	typedef NodeParam Node; // Type of tree nodes
	
	struct HeapItem // Structure to store heap items
		{
		/* Elements: */
		public:
		Node* node; // Pointer to node stored in item
		unsigned int renderPass; // Counter value of last rendering pass that touched the node
		float LOD; // LOD value of node in last rendering pass that touched it
		
		/* Constructors and destructors: */
		HeapItem(void)
			:node(0),renderPass(0U),LOD(0.0f)
			{
			};
		
		/* Methods: */
		friend bool operator<=(const HeapItem& item1,const HeapItem& item2) // Comparison operator for heap ordering
			{
			if(item1.renderPass==item2.renderPass)
				return item1.LOD<=item2.LOD;
			else
				return item1.renderPass<item2.renderPass;
			};
		};
	
	/* Elements: */
	HeapItem* items; // Array of heap items
	unsigned int numItems; // Current number of items in heap
	
	/* Constructors and destructors: */
	public:
	CoarseningHeap(unsigned int maxNumItems) // Creates heap for given number of items; caller's responsibility to never insert more items
		:items(new HeapItem[maxNumItems]),
		 numItems(0U)
		{
		};
	~CoarseningHeap(void) // Destroys heap
		{
		delete[] items;
		};
	
	/* Methods: */
	unsigned int getNumItems(void) const // Returns current number of items in heap
		{
		return numItems;
		};
	Node* getTopNode(void) const // Returns the node most appropriate for coarsening
		{
		return numItems>0U?items[0].node:0;
		};
	void insert(Node* newNode) // Inserts the given node into the heap; takes renderPass and LOD from node, updates node's coarseningHeapIndex
		{
		/* Create a new heap item: */
		HeapItem newItem;
		newItem.node=newNode;
		newItem.renderPass=newNode->renderPass;
		newItem.LOD=newNode->maxLOD;
		
		/* Insert new items at bottom of heap and let them percolate up: */
		unsigned int insertionPos=numItems;
		while(insertionPos>0U)
			{
			unsigned int parent=(insertionPos-1U)>>1;
			if(items[parent]<=newItem)
				break;
			items[insertionPos]=items[parent];
			items[insertionPos].node->coarseningHeapIndex=insertionPos;
			insertionPos=parent;
			}
		items[insertionPos]=newItem;
		items[insertionPos].node->coarseningHeapIndex=insertionPos;
		++numItems;
		};
	void remove(Node* node) // Removes given node from heap
		{
		/* Remove node by inserting item from bottom of heap at node's position: */
		--numItems;
		unsigned int insertionPos=node->coarseningHeapIndex;
		node->coarseningHeapIndex=~0x0U;
		
		/* Start by letting the bottom item percolate up the heap: */
		while(insertionPos>0U)
			{
			unsigned int parent=(insertionPos-1U)>>1;
			if(items[parent]<=items[numItems])
				break;
			items[insertionPos]=items[parent];
			items[insertionPos].node->coarseningHeapIndex=insertionPos;
			insertionPos=parent;
			}
		
		/* Then let the bottom item trickle down to its final position: */
		while(insertionPos<numItems)
			{
			unsigned int minIndex=numItems;
			unsigned int child=(insertionPos<<1)+1;
			if(child<=numItems&&!(items[minIndex]<=items[child]))
				minIndex=child;
			++child;
			if(child<=numItems&&!(items[minIndex]<=items[child]))
				minIndex=child;
			items[insertionPos]=items[minIndex];
			items[insertionPos].node->coarseningHeapIndex=insertionPos;
			insertionPos=minIndex;
			}
		};
	void move(Node* node) // Updates a node's heap data and moves the node to its proper position in the heap
		{
		/* Update the node's data: */
		unsigned int insertionPos=node->coarseningHeapIndex;
		items[insertionPos].renderPass=node->renderPass;
		items[insertionPos].LOD=node->maxLOD;
		
		/* Start by letting the node's item percolate up the heap: */
		while(insertionPos>0U)
			{
			unsigned int parent=(insertionPos-1U)>>1;
			if(items[parent]<=items[insertionPos])
				break;
			HeapItem temp=items[insertionPos];
			items[insertionPos]=items[parent];
			items[insertionPos].node->coarseningHeapIndex=insertionPos;
			items[parent]=temp;
			items[parent].node->coarseningHeapIndex=parent;
			insertionPos=parent;
			}
		
		/* Then let the node's item trickle down to its final position: */
		while(true)
			{
			unsigned int minIndex=insertionPos;
			unsigned int child=(insertionPos<<1)+1;
			if(child<numItems&&!(items[minIndex]<=items[child]))
				minIndex=child;
			++child;
			if(child<numItems&&!(items[minIndex]<=items[child]))
				minIndex=child;
			if(minIndex==insertionPos)
				break;
			HeapItem temp=items[insertionPos];
			items[insertionPos]=items[minIndex];
			items[insertionPos].node->coarseningHeapIndex=insertionPos;
			items[minIndex]=temp;
			items[minIndex].node->coarseningHeapIndex=minIndex;
			insertionPos=minIndex;
			}
		};
	bool checkHeap(void) const // Checks the heap for internal consistency
		{
		for(unsigned int i=0;i<numItems;++i)
			{
			const Node* node=items[i].node;
			
			/* Check if the node's heap index matches its position: */
			if(node->coarseningHeapIndex!=i)
				return false;
			
			/* Check if the node's data in the heap matches the node's data: */
			if(node->renderPass!=items[i].renderPass||node->maxLOD!=items[i].LOD)
				return false;
			
			/* Check the heap invariant: */
			unsigned int child=(i<<1)+1U;
			if(child<numItems&&!(items[i]<=items[child]))
				return false;
			++child;
			if(child<numItems&&!(items[i]<=items[child]))
				return false;
			}
		
		return true;
		};
	};

#endif
