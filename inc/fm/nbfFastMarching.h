#ifndef FILE_nbfFastMarching
#define FILE_nbfFastMarching

// Class nbfFastMarching.
//
// Implements 2D and 3D fast marching method with arbitrary weights.
// 

#include <nbfArray.h>

#include <vector>

#include <queue>

// #include <fm/DiscritePriorityQueue.h>

template< class Pixel, const int Dim >
struct greaterPD : public binary_function< TinyVector< Pixel, Dim + 1 >, TinyVector< Pixel, Dim + 1 >, bool >
{
  bool operator()( const TinyVector< Pixel, Dim + 1 >& x, const TinyVector< Pixel, Dim + 1 >& y) const 
    {
      return x[ Dim ] > y[ Dim ];
    }
};

template< class T, class Cmp, class C = vector<T> >
class PriorityQueue : public std::priority_queue< T,C,Cmp>
{
public :
  void clear(){
    this->c.clear();
   };
  explicit PriorityQueue() : std::priority_queue<T,C,Cmp>(){};
  
};

template< class Pixel, int const Dim >
class nbfFastMarching : public nbfArray< Pixel, Dim >
{
public:

	static const int ACCEPTED_SET = 0, FAR_SET = 1, TRIAL_SET = 2;

	// constructor takes weight array as input
	nbfFastMarching( Array< Pixel, Dim > & );
	nbfFastMarching( Array< Pixel, Dim > &, TinyVector< int, Dim > & );

	//~nbfFastMarching(){this->pqueue->Delete();};
	~nbfFastMarching(){};

	void construct( Array< Pixel, Dim > & );

	// fast marching initalization
	void inicialize();

	// set individual point and corresponding distance as initial set
	void setAliveSet( TinyVector< int, Dim >, Pixel );

	// set multiple points and corresponding distances as initial set
	void setAliveSet( vector< TinyVector< int, Dim > > &,  vector< Pixel > & = 0 );

	// set new set of weights
	void setWeights( Array< Pixel, Dim > & );

	// run fast marching, write result on argument
	virtual TinyVector< int, Dim > execute( Array< Pixel, Dim > & );

	// set freezing threshold value
	void setFreezing( Pixel th ){ this->dth = th; };

	void updateQueue( TinyVector< int, Dim > & );

	// Stopping conditions

	// set stoping point coordinates
	void setStopPoint( TinyVector< int, Dim > & );

	// unset stoping point flag
	void unSetStopPoint();

	// set stoping distance
	void setStopDistance( Pixel );

	// unset stoping distance flag
	void unSetStopDistance();

	// set stoping border at given dimension
	void setStopBorder( int );

	// unset stoping border flag
	void unSetStopBorder();

protected:

	// get valid 4 or 6 neighbors of a point
	virtual void getNeighbors( TinyVector< int, Dim > &, vector< TinyVector< int, Dim > > & );
	
	// eliminate neighbors that are already ALIVE 
	void checkNeighbors( vector< TinyVector< int, Dim > > & );

	// update point neighbors
	void updateNeighbors( TinyVector< int, Dim > & );

	// update individual point
	// the first argument is the point to update
	// the second argument is the point from where we do the update
	virtual void updatePoint( TinyVector< int, Dim > &, TinyVector< int, Dim > & ) = 0;

	// get valid distance value
	Pixel getDistance( TinyVector< int, Dim > & );

	// get valid distance value
	int getState( TinyVector< int, Dim > & );

	// convert back and fort linear and matrix indexes
	virtual int array2int( TinyVector< int, Dim > & ) = 0;
	virtual TinyVector< int, Dim > int2array( int ) = 0;

	// void force stop
	virtual bool checkForceStop( TinyVector< int, Dim > & position ){ return false; };

	// atributes

	// priority queue
	//vtkPriorityQueue * pqueue;

	// DiscritePriorityQueue< TinyVector< Pixel, 3 > > queueLiron;

	PriorityQueue< TinyVector< Pixel, Dim + 1 >, greaterPD< Pixel, Dim > , vector< TinyVector< Pixel, Dim + 1 > > > pqNB;

	// data dimensions
	TinyVector< int, Dim > dimensions;

	// initial alive set coordinates and distances
	vector< TinyVector< int, Dim > > aliveSetPositions;
	vector< Pixel > aliveSetDistances;

	// point stoping
	bool flagPuntoParada;
	TinyVector< int, Dim > PuntoParada;

	// distance stoping
	bool flagDistanciaParada;
	Pixel DistanciaParada;

	// border stoping
	bool flagBordeParada;
	int BordeParada;

	// force stop
	bool flagForceStop;

	// distances (taken from execute())
	Array< Pixel, Dim > distance;

	// state array = ACCEPTED_SET, FAR_SET, TRIAL_SET
	Array< int, Dim > state;
	
	// weights (taken from constructor)
	Array< Pixel, Dim > weight;

	// arrival times
	Array< int, Dim > times;

	// freezing parameters
	Pixel dmax;
	Pixel dth;

};


template< class Pixel, int const Dim >
nbfFastMarching< Pixel, Dim > :: nbfFastMarching( Array< Pixel, Dim > & weight, TinyVector< int, Dim > & center )
: nbfArray< Pixel, Dim >( weight, center )
//, queueLiron(1000,1)
{
	this->construct( weight );
}

template< class Pixel, int const Dim >
nbfFastMarching< Pixel, Dim > :: nbfFastMarching( Array< Pixel, Dim > & weight )
: nbfArray< Pixel, Dim >( weight )
{
	this->construct( weight );
}

template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: construct( Array< Pixel, Dim > & weight )
{
	this->weight.reference( weight );
	this->dimensions = weight.shape();
	this->state.resize( this->dimensions );

    //this->pqueue = vtkPriorityQueue::New();

	// freezing: thresholding value
	this->dth = 0;

	// reset all stoping flags

	this->PuntoParada = 0;
	this->flagPuntoParada = false;

	this->DistanciaParada = numeric_limits<Pixel>::max();
	this->flagDistanciaParada = false;

	this->BordeParada = 0;
	this->flagBordeParada = false;

	this->flagForceStop = false;
}


template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: setStopPoint( TinyVector< int, Dim > & pp )
{
	flagDistanciaParada = false;
	flagBordeParada = false;

	PuntoParada = pp;
	flagPuntoParada = true;
}

template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: unSetStopPoint()
{
	PuntoParada = Dim;
	flagPuntoParada = false;
}

template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: setStopDistance( Pixel dp )
{
	flagPuntoParada = false;
	flagBordeParada = false;

	DistanciaParada = dp;
	flagDistanciaParada = true;
}

template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: unSetStopDistance()
{
  DistanciaParada = numeric_limits< Pixel > :: max();
	flagDistanciaParada = false;
}

template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: setStopBorder( int db )
{
	flagDistanciaParada = false;
	flagPuntoParada = false;

	BordeParada = db;
	flagBordeParada = true;
}

template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: unSetStopBorder()
{
	BordeParada = 0;
	flagBordeParada = false;
}


template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: setAliveSet( TinyVector< int, Dim > position,
													 Pixel distance )
{
	vector< TinyVector< int, Dim > > vPositions;
	vector< Pixel > vDistances;

	vPositions.push_back( position );
	vDistances.push_back( distance );

	this->setAliveSet( vPositions, vDistances );
}


template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: setWeights( Array< Pixel, Dim > & weight )
{
	this->weight.reference( weight );
}

template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: setAliveSet( vector< TinyVector< int, Dim > > & positions,
													 vector< Pixel > & distances )
{
	// reset attributes

	this->aliveSetPositions.clear();
	this->aliveSetDistances.clear();

	// deep copy both arguments into attributes

	typename vector< TinyVector< int, Dim > > :: iterator iterPositions;
	typename vector< Pixel > :: iterator iterDistances;

	for ( iterPositions = positions.begin(), iterDistances = distances.begin();
		( iterPositions != positions.end() ) && ( iterDistances != distances.end());
		++iterPositions, ++iterDistances )
	{
		// if it is inside image range
		if ( this->weight.isInRange( *iterPositions ) ){
			this->aliveSetPositions.push_back( *iterPositions );
			this->aliveSetDistances.push_back( *iterDistances );
		}
	}
}


template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: inicialize()
{
	// reset priority queue
	//this->pqueue->Reset();
	//this->pqueue->Allocate( this->distance.numElements() );

	this->pqNB.clear();

	// set all distances to infinity
	this->distance = numeric_limits<Pixel>::max();

	if ( this->dth != 0 ){
		this->times.resize( this->dimensions );
		this->times = numeric_limits<Pixel>::max();
	}

	// set all states to FAR
	this->state = FAR_SET;

	// initialize all points in initial alive set
	typename vector< TinyVector< int, Dim > > :: iterator iterPositions;
	typename vector< Pixel > :: iterator iterDistances;

	for ( iterPositions = this->aliveSetPositions.begin(), 
		iterDistances = this->aliveSetDistances.begin(); 
		iterPositions != this->aliveSetPositions.end(),
		iterDistances != this->aliveSetDistances.end();
	++iterPositions, ++iterDistances ){

		// set distance value
		this->distance( *iterPositions ) = *iterDistances;

		// move to ACCEPTED set
		this->state( *iterPositions ) = ACCEPTED_SET;

		if ( this->dth != 0 ){
			this->times( *iterPositions ) = 0;
		}
	}

	// update neighbors of points in the initial ACCEPTED set
	for ( iterPositions = this->aliveSetPositions.begin();
		iterPositions != this->aliveSetPositions.end();
		++iterPositions ){
			this->updateNeighbors( *iterPositions );
		}
}


template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: getNeighbors( TinyVector< int, Dim > & currentPoint,
													vector< TinyVector< int, Dim > > & neighbors )
{
	nbfArray< Pixel, Dim > :: getNeighbors( currentPoint, neighbors );

	this->checkNeighbors( neighbors );

	//vector< TinyVector< int, Dim > > :: iterator iter = neighbors.begin();

	//while ( iter != neighbors.end() ){
	//	if ( this->state( *iter ) == ACCEPTED_SET ){
	//		iter = neighbors.erase( iter );
	//	}
	//	else{
	//		++iter;
	//	}
	//}
}

template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: checkNeighbors( vector< TinyVector< int, Dim > > & neighbors )
{
	typename vector< TinyVector< int, Dim > > :: iterator iter = neighbors.begin();

	while ( iter != neighbors.end() ){
		if ( ( this->state( *iter ) == ACCEPTED_SET ) || ( this->weight( *iter ) == numeric_limits<Pixel>::max() ) ){
			iter = neighbors.erase( iter );
		}
		else{
			++iter;
		}
	}
}

template< class Pixel, int const Dim >
Pixel nbfFastMarching< Pixel, Dim > :: getDistance( TinyVector< int, Dim > & point )
{
	// if inside the array
	if ( this->distance.isInRange( point ) ){
		return ( this->distance( point ) );
	}
	// else infinity
	else{
		return numeric_limits<Pixel>::max();
	}
}


template< class Pixel, int const Dim >
int nbfFastMarching< Pixel, Dim > :: getState( TinyVector< int, Dim > & point )
{
	// if inside the array
	if ( this->state.isInRange( point ) ){
		return ( this->state( point ) );
	}
	// else infinity
	else{
		return FAR_SET;
	}
}


template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: updateNeighbors( TinyVector< int, Dim > & currentPoint )
{
	// get valid neighbors first
	vector< TinyVector< int, Dim > > neighbors;
	this->getNeighbors( currentPoint, neighbors );

	// update each neighbor
	for ( unsigned int i = 0; i < neighbors.size(); i++ ){
		this->updatePoint( neighbors[i], currentPoint );
		
		if ( this->dth != 0 ){
            // freezing: update arrival times
			this->times( neighbors[i] ) = blitz::extrema::min( this->times( currentPoint ) + 1, this->times( neighbors[i] ) );
		}
	}
}

template< class Pixel, int const Dim >
void nbfFastMarching< Pixel, Dim > :: updateQueue( TinyVector< int, Dim > & current )
{
	// get global image ID
	//int currentId = this->array2int(current);

	// if in FAR_SET move to TRIAL_SET
	if ( this->state( current ) == FAR_SET ){
		this->state( current ) = TRIAL_SET;
	}
	//// if already in TRIAL_SET, remove from queue
	//else if ( this->state( current ) == TRIAL_SET ){
	//	float p = this->pqueue->DeleteId( currentId );
	//	if ( p == VTK_LARGE_FLOAT ){
	//		cout << "point not in queue" << endl;
	//	}
	//}

	// insert in queue
	//this->pqueue->Insert( this->distance(current), currentId );
	TinyVector< Pixel, Dim + 1 > next;
	if ( Dim == 2 ){
		next[firstDim] = current[firstDim];
		next[secondDim] = current[secondDim];
		next[thirdDim] = this->distance(current);
	}
	if ( Dim == 3 ){
		next[firstDim] = current[firstDim];
		next[secondDim] = current[secondDim];
		next[thirdDim] = current[thirdDim];
		next[fourthDim] = this->distance(current);
	}

	this->pqNB.push( next );
	//if ( this->distance(current) < numeric_limits< Pixel > :: max() ){
	//	this->queueLiron.insert(this->distance(current)*998+.5,next);
	//}
}

template< class Pixel, int const Dim >
TinyVector< int, Dim > nbfFastMarching< Pixel, Dim > :: execute( Array< Pixel, Dim > & output )
{
	TinyVector< int, Dim > currentPosition;

	// if initial alive set not empty
	if( this->aliveSetPositions.size() != 0 )
	{
		// initialize state and distances
		this->distance.reference( output );
		this->inicialize();

		// freezing: reset furthest point
		this->dmax = 0;

		// stoping criteria
		int stopPosicion;
		bool stopDistancia;
		bool stopBorde = false;

		do {

			TinyVector< Pixel, Dim + 1 > firstOut;
			TinyVector< int, Dim > firstPos;
			do {
				firstOut = pqNB.top();
				pqNB.pop();
				firstPos(0) = firstOut(0);
				firstPos(1) = firstOut(1);
				if ( Dim > 2 ){
					firstPos(2) = firstOut(2);
				}
			} while( ( this->state( firstPos ) == ACCEPTED_SET ) && ( pqNB.size() != 0 ) );

			//cout << pqNB.size() << endl;
			 //queueLiron
			//do {
			//	this->queueLiron.pop(firstOut);
			//	firstPos(0) = floor( firstOut(0) );
			//	firstPos(1) = floor( firstOut(1) );
			//} while( ( this->state( firstPos ) == ACCEPTED_SET ) && ( !this->queueLiron.isEmpty() ) );

			// get point in TRIAL with the smallest distance value

			// remove it from the priority queue
			float vdistance;
			//int id = this->pqueue->Pop(0,vdistance);	 // Pop first value.

			//// move it to the ACCEPTED set
			//currentPosition = this->int2array(id);

			currentPosition = firstPos;
			vdistance = this->distance( currentPosition );

			this->state( currentPosition ) = ACCEPTED_SET;

			//cout << currentPosition << endl;
			//this->flagForceStop = this->checkForceStop( currentPosition );

			//// if freeze
			//if ( this->dth && ( this->times( currentPosition ) < this->dmax - this->dth ) ){
			//	this->times( currentPosition ) = numeric_limits<Pixel>::max();
			//	this->distance( currentPosition ) = numeric_limits<Pixel>::max();
			//}
			//else{
			//	// only if using freezing
			//	if ( this->dth != 0 ){
			//		// update furthest point
			//		this->dmax = blitz::minmax::max( this->times( currentPosition ), this->dmax );			
			//	}

			//	// update neighbours
			//	this->updateNeighbors( currentPosition );
			//}

			// update neighbours
			this->updateNeighbors( currentPosition );

			// check stop conditions
			stopPosicion = sum( abs( currentPosition - PuntoParada ) );
			//if ( !this->distance.isInRange( currentPosition ) ){
			//	cout << "fastmarching error " << currentPosition << endl;
			//}
			stopDistancia = ( vdistance - DistanciaParada ) < 0;
			//if ( this->hasCut ){
			//	stopBorde = !( ( currentPosition(firstDim) > this->center[firstDim] ) & ( currentPosition(secondDim) == ( this->center[secondDim] + 1 ) ) );
			//}
			//else{
			//	stopBorde = ( currentPosition(BordeParada) < output.ubound( BordeParada ) );
			//}
			////stopBorde = ( pd.position(BordeParada) > output.lbound( BordeParada ) );

		//} while ( ( this->pqueue->GetNumberOfItems() != 0 ) &&
		} while ( ( this->pqNB.size() != 0 ) &&
		//} while ( ( !this->queueLiron.isEmpty() ) &&
			( !flagForceStop ) &&
			( ( stopPosicion != 0 ) || !( flagPuntoParada ) ) &&
			( ( stopDistancia != 0 ) || !( flagDistanciaParada ) ) &&
			( ( stopBorde != 0 ) || !( flagBordeParada ) ) );
	}
	else{
		output = numeric_limits<Pixel>::max();
	}
	return currentPosition;
}

#endif /* FILE_nbfFastMarching */
