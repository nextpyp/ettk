#ifndef DISCRITE_PRIORITY_QUEUE_H
#define DISCRITE_PRIORITY_QUEUE_H

#define PREALLOCATE_SIZE 128

#ifndef NULL
#define NULL (0)
#endif

template<class T> struct DPQNode;
template<class T> struct DPQSeedNode;

template<class T, class Dtype=float>
class DiscritePriorityQueue {
public:

	DiscritePriorityQueue(int numBins, Dtype factor, Dtype initPriority=0);

	// insert a new element to the priority queue. The priority
	// inserted may not be small than the smallest priority last poped from
	// queue, and may not be larger than the smallest priority last poped
	// from queue + the numBins
	//inline void insert(long priority, const T& data);
	inline void insert(Dtype priority, const T& data);

	// can insert any priority larger the last poped
	inline void insertSeed(Dtype priority, const T& data);

	// gives the element with the smallest priority. Returns true if such
	// element what found. False if queue is empty.
	inline bool pop(T& data);

	inline bool isEmpty() { return numElements==0; }

	///////////////////////////////////////////////////////////////////////
private:
	int numBins;

	// create the cyclic table
	DPQNode<T>** tail;
	int minPos;
	long minPriority;
	Dtype factor;

	int numElements;

	// fast memory allocation and access
	DPQNode<T>* pool;
	DPQNode<DPQNode<T>* >* preAllocateList;
	int preAllocateSize;
	void preAllocate();

	// the case where seeds are far distance
	DPQSeedNode<T>* seedNodes;

	inline DPQNode<T>* getNodeFromPool();
	inline void returnNodeToPool(DPQNode<T>* node);
	void addSeedsFromNodes();



public:
	//DEBUG
	//int debugStatPops;
	//int debugStatSkips;
	//int debugNewCount;


public:
	~DiscritePriorityQueue();
};

//////////////////////////// IMPLEMENTATION ////////////////////////////////

template<class T>
struct DPQNode {
	DPQNode<T>* next;
	T data;
};

template<class T>
struct DPQSeedNode {
	DPQSeedNode* next;
	T data;
	long discretePriority;
};


#include <assert.h>

template<class T, class Dtype>
DiscritePriorityQueue<T,Dtype>::DiscritePriorityQueue(int numBins_, Dtype factor_, Dtype initPriority) 
: numBins(numBins_), minPos(0),  numElements(0), pool(NULL), seedNodes(NULL),
preAllocateSize(PREALLOCATE_SIZE), preAllocateList(NULL) {

	factor=(numBins-2)/factor_;  // the 2 is for Beto, better 1
	minPriority=(long)(initPriority*factor);

	tail = new DPQNode<T>*[numBins];

	for(int i=0; i<numBins; ++i) {
		tail[i]=NULL;
	}

	preAllocate();
}


template<class T, class Dtype>
DiscritePriorityQueue<T,Dtype>::~DiscritePriorityQueue() {
	delete[] tail;
	
	while (preAllocateList!=NULL) {
		delete[] preAllocateList->data;
		preAllocateList=preAllocateList->next;
	}

	while(seedNodes!=NULL) {
		DPQSeedNode<T>* oldNode=seedNodes;
		seedNodes=seedNodes->next;
		delete oldNode;
	}
}

template<class T, class Dtype>
void DiscritePriorityQueue<T,Dtype>::insertSeed(Dtype priority, const T& data) {
	long discretePriority=(long)(priority*factor);
	int reletivePriority= discretePriority - (long)minPriority;
	assert(reletivePriority>=0);

	// if out of range for discrete queue use extern sorted linked list
	if (reletivePriority>=numBins) {
		// create a new node for the Seed link list
		DPQSeedNode<T>* newNode=new DPQSeedNode<T>;
		newNode->discretePriority=discretePriority;
		newNode->data=data;
		
		// add to link list
		if (seedNodes==NULL) {
			newNode->next=NULL;
			seedNodes=newNode;
		} else {
			DPQSeedNode<T>* ptr=seedNodes;
			while (ptr->next!=NULL && ptr->next->discretePriority <discretePriority) {
				ptr=ptr->next;
			}

			if (ptr==seedNodes && ptr->discretePriority>discretePriority) {
				seedNodes=newNode;
				newNode->next=ptr;
			} else {
				newNode->next=ptr->next;
				ptr->next=newNode;
			}
		}
	} else insert(priority, data);
}



template<class T, class Dtype>
void DiscritePriorityQueue<T,Dtype>::insert(Dtype priority, const T& data) {

	int reletivePriority= (long)(priority*factor) - (long)minPriority;

	assert(reletivePriority>=0 && reletivePriority<numBins);

	// find the cyclic position of the reletive priority
	int pos= reletivePriority+minPos;
	if (pos>=numBins) pos-=numBins;

	DPQNode<T>* newNode= getNodeFromPool(); 

	DPQNode<T>* tailPos=*(tail+pos);
	
	if (tailPos!=NULL) {
		newNode->next=tailPos->next;
		tailPos->next=newNode;
	} else {
		newNode->next=newNode;
	}

	tail[pos]=newNode;

	newNode->data=data;

	++numElements;
}


template<class T, class Dtype>
bool DiscritePriorityQueue<T,Dtype>::pop(T& data) {
	if (isEmpty()) {
		if (seedNodes!=NULL) {
			minPriority=seedNodes->discretePriority;
			addSeedsFromNodes();
		} else return false;
	}

	bool skippedLevel=false;
	while ( tail[minPos]==NULL) {
		++minPos;
		if (minPos==numBins) minPos=0;
		++minPriority;
		skippedLevel=true;
	}

	if (skippedLevel) {
		addSeedsFromNodes();
	}

	data=tail[minPos]->next->data;
	DPQNode<T>* head=tail[minPos]->next;
	if (tail[minPos]==head) {
		tail[minPos]=NULL;
	} else {
		tail[minPos]->next=head->next;
	}

	returnNodeToPool(head);
	--numElements;
	return true;
}

template<class T, class Dtype>
void DiscritePriorityQueue<T,Dtype>::addSeedsFromNodes() {
	while(seedNodes!=NULL && minPriority+numBins-1>=seedNodes->discretePriority) {
		insert(seedNodes->discretePriority/factor,seedNodes->data);
		DPQSeedNode<T>* oldNode=seedNodes;
		seedNodes=seedNodes->next;
		delete oldNode;
	}
}

template<class T, class Dtype>
DPQNode<T>* DiscritePriorityQueue<T,Dtype>::getNodeFromPool() {
	DPQNode<T>* node;
	if (pool==NULL) {
		preAllocate();
	}
	node= pool;
	pool=pool->next;
	return node;
}

template<class T, class Dtype>
void DiscritePriorityQueue<T,Dtype>::returnNodeToPool(DPQNode<T>* node) {
	node->next=pool;
	pool=node;
}

template<class T, class Dtype>
void DiscritePriorityQueue<T,Dtype>::preAllocate() {
	DPQNode<DPQNode<T>* >* newAlloc = new DPQNode<DPQNode<T>* >;
	newAlloc->next=preAllocateList;
	preAllocateList=newAlloc;
	newAlloc->data= new DPQNode<T>[preAllocateSize];

	DPQNode<T>* temp= newAlloc->data;
	for (int i=1; i<preAllocateSize; ++i) {
		temp->next=pool;
		pool=temp;
		++temp;
	}

	preAllocateSize<<=2;
}


#endif
