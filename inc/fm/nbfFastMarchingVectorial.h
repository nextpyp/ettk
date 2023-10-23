#ifndef FILE_nbfFastMarching
#define FILE_nbfFastMarching

// Class nbfFastMarching.
//
// Implements 2D and 3D fast marching method with arbitrary weights.
// 

#include <vtkPriorityQueue.h>

#include <nbfArray.h>
#include <nbfMetric.h>

#include <vector>

template< class Pixel, int const Dim, class Scalar = Pixel >
class nbfFastMarching : public nbfArray< Pixel, Dim >
{
public:

	static const int ACCEPTED_SET = 0, FAR_SET = 1, TRIAL_SET = 2;

	// constructor takes weight array as input
	nbfFastMarching( Array< Pixel, Dim > & );
	nbfFastMarching( Array< Pixel, Dim > &, TinyVector< int, Dim > & );

	~nbfFastMarching(){this->pqueue->Delete();};

	void construct( Array< Pixel, Dim > & );

	// set individual point and corresponding distance as initial set
	void setAliveSet( TinyVector< int, Dim >, Scalar );

	// set multiple points and corresponding distances as initial set
	void setAliveSet( vector< TinyVector< int, Dim > > &,  vector< Scalar > & = 0 );

	// set new set of weights
	void setWeights( Array< Pixel, Dim > & );

	// run fast marching, write result on argument
	virtual TinyVector< int, Dim > execute( Array< Scalar, Dim > & );

	// set freezing threshold value
	void setFreezing( Scalar th ){ this->dth = th; };

	void updateQueue( TinyVector< int, Dim > & );

	// Stopping conditions

	// set stoping point coordinates
	void setStopPoint( TinyVector< int, Dim > & );

	// unset stoping point flag
	void unSetStopPoint();

	// set stoping distance
	void setStopDistance( Scalar );

	// unset stoping distance flag
	void unSetStopDistance();

	// set stoping border at given dimension
	void setStopBorder( int );

	// unset stoping border flag
	void unSetStopBorder();

protected:

	// fast marching initalization
	virtual void initialize();

	// fast marching loop update (if we want to do something inside the loop)
	virtual void loopUpdate( TinyVector< int, Dim > & ){};

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
	Scalar getDistance( TinyVector< int, Dim > & );

	// get valid distance value
	int getState( TinyVector< int, Dim > & );

	// convert back and fort linear and matrix indexes
	virtual int array2int( TinyVector< int, Dim > & ) = 0;
	virtual TinyVector< int, Dim > int2array( int ) = 0;

	// void force stop
	virtual bool checkForceStop( TinyVector< int, Dim > & position ){ return false; };

	// atributes

	// priority queue
	vtkPriorityQueue * pqueue;

	// data dimensions
	TinyVector< int, Dim > dimensions;

	// initial alive set coordinates and distances
	vector< TinyVector< int, Dim > > aliveSetPositions;
	vector< Scalar > aliveSetDistances;

	// point stoping
	bool flagPuntoParada;
	TinyVector< int, Dim > PuntoParada;

	// distance stoping
	bool flagDistanciaParada;
	Scalar DistanciaParada;

	// border stoping
	bool flagBordeParada;
	int BordeParada;

	// force stop
	bool flagForceStop;

	// distances (taken from execute())
	Array< Scalar, Dim > distance;

	// state array = ACCEPTED_SET, FAR_SET, TRIAL_SET
	Array< int, Dim > state;
	
	// weights (taken from constructor)
	Array< Pixel, Dim > weight;

	// arrival times
	Array< int, Dim > times;

	// freezing parameters
	Scalar dmax;
	Scalar dth;

	// \infty value independent of representation
	Pixel infty;
};


template< class Pixel, int const Dim, class Scalar >
nbfFastMarching< Pixel, Dim, Scalar > :: nbfFastMarching( Array< Pixel, Dim > & weight, TinyVector< int, Dim > & center )
: nbfArray< Pixel, Dim >( weight, center )
{
	this->construct( weight );
}

template< class Pixel, int const Dim, class Scalar >
nbfFastMarching< Pixel, Dim, Scalar > :: nbfFastMarching( Array< Pixel, Dim > & weight )
: nbfArray< Pixel, Dim >( weight )
{
	this->construct( weight );
}

template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: construct( Array< Pixel, Dim > & weight )
{
	this->weight.reference( weight );
	this->dimensions = weight.shape();
	this->state.resize( this->dimensions );
	this->times.resize( this->dimensions );

    this->pqueue = vtkPriorityQueue::New();

	// freezing: thresholding value
	this->dth = 0;

	// infty value
	this->infty = numeric_limits<Scalar>::max();

	// reset all stoping flags

	this->PuntoParada = 0;
	this->flagPuntoParada = false;

	this->DistanciaParada = numeric_limits<Scalar>::max();
	this->flagDistanciaParada = false;

	this->BordeParada = 0;
	this->flagBordeParada = false;

	this->flagForceStop = false;
}


template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: setStopPoint( TinyVector< int, Dim > & pp )
{
	flagDistanciaParada = false;
	flagBordeParada = false;

	PuntoParada = pp;
	flagPuntoParada = true;
}

template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: unSetStopPoint()
{
	PuntoParada = dimensiones;
	flagPuntoParada = false;
}

template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: setStopDistance( Scalar dp )
{
	flagPuntoParada = false;
	flagBordeParada = false;

	DistanciaParada = dp;
	flagDistanciaParada = true;
}

template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: unSetStopDistance()
{
	DistanciaParada = infinito;
	flagDistanciaParada = false;
}

template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: setStopBorder( int db )
{
	flagDistanciaParada = false;
	flagPuntoParada = false;

	BordeParada = db;
	flagBordeParada = true;
}

template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: unSetStopBorder()
{
	BordeParada = 0;
	flagBordeParada = false;
}


template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: setAliveSet( TinyVector< int, Dim > position,
													 Scalar distance )
{
	vector< TinyVector< int, Dim > > vPositions;
	vector< Scalar > vDistances;

	vPositions.push_back( position );
	vDistances.push_back( distance );

	this->setAliveSet( vPositions, vDistances );
}


template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: setWeights( Array< Pixel, Dim > & weight )
{
	this->weight.reference( weight );
}

template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: setAliveSet( vector< TinyVector< int, Dim > > & positions,
													 vector< Scalar > & distances )
{
	// reset attributes

	this->aliveSetPositions.clear();
	this->aliveSetDistances.clear();

	// deep copy both arguments into attributes

	vector< TinyVector< int, Dim > > :: iterator iterPositions;
	vector< Scalar > :: iterator iterDistances;

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


template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: initialize()
{
	// reset priority queue
	this->pqueue->Reset();
	this->pqueue->Allocate( this->distance.numElements() );

	// set all distances to infinity
	this->distance = numeric_limits<Scalar>::max();

	this->times = numeric_limits<Scalar>::max();

	// set all states to FAR
	this->state = FAR_SET;

	// initialize all points in initial alive set
	vector< TinyVector< int, Dim > > :: iterator iterPositions;
	vector< Scalar > :: iterator iterDistances;

	for ( iterPositions = this->aliveSetPositions.begin(), 
		iterDistances = this->aliveSetDistances.begin(); 
		iterPositions != this->aliveSetPositions.end(),
		iterDistances != this->aliveSetDistances.end();
	++iterPositions, ++iterDistances ){

		// set distance value
		this->distance( *iterPositions ) = *iterDistances;

		// move to ACCEPTED set
		this->state( *iterPositions ) = ACCEPTED_SET;

		this->times( *iterPositions ) = 0;
	}

	// update neighbors of points in the initial ACCEPTED set
	for ( iterPositions = this->aliveSetPositions.begin();
		iterPositions != this->aliveSetPositions.end();
		++iterPositions ){
			this->updateNeighbors( *iterPositions );
		}
}


template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: getNeighbors( TinyVector< int, Dim > & currentPoint,
													vector< TinyVector< int, Dim > > & neighbors )
{
	nbfArray< Pixel, Dim > :: getNeighbors( currentPoint, neighbors );

	this->checkNeighbors( neighbors );
}

template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: checkNeighbors( vector< TinyVector< int, Dim > > & neighbors )
{
	vector< TinyVector< int, Dim > > :: iterator iter = neighbors.begin();

	while ( iter != neighbors.end() ){
		// workaround for vector data type
		Pixel tmp( this->weight( *iter ) / numeric_limits<Scalar>::max() );

		// check if already accepted or \infty weight
		//if ( ( this->state( *iter ) == ACCEPTED_SET ) | ( pow2(tmp) ==  1 ) ){
		bool isInfty = nbfMetric::mod( this->weight( *iter ) ) == numeric_limits<Scalar>::max();
		if ( ( this->state( *iter ) == ACCEPTED_SET ) | isInfty ){
			iter = neighbors.erase( iter );
		}
		else{
			++iter;
		}
	}
}

template< class Pixel, int const Dim, class Scalar >
Scalar nbfFastMarching< Pixel, Dim, Scalar > :: getDistance( TinyVector< int, Dim > & point )
{
	// if inside the array
	if ( this->distance.isInRange( point ) ){
		return ( this->distance( point ) );
	}
	// else infinity
	else{
		return numeric_limits<Scalar>::max();
	}
}


template< class Pixel, int const Dim, class Scalar >
int nbfFastMarching< Pixel, Dim, Scalar > :: getState( TinyVector< int, Dim > & point )
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


template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: updateNeighbors( TinyVector< int, Dim > & currentPoint )
{
	// get valid neighbors first
	vector< TinyVector< int, Dim > > neighbors;
	this->getNeighbors( currentPoint, neighbors );

	// update each neighbor
	for ( unsigned int i = 0; i < neighbors.size(); i++ ){
		this->updatePoint( neighbors[i], currentPoint );
		
		if ( this->dth != 0 ){
            // freezing: update arrival times
			this->times( neighbors[i] ) = blitz::minmax::min( this->times( currentPoint ) + 1, this->times( neighbors[i] ) );
		}
	}
}

template< class Pixel, int const Dim, class Scalar >
void nbfFastMarching< Pixel, Dim, Scalar > :: updateQueue( TinyVector< int, Dim > & current )
{
	// get global image ID
	int currentId = this->array2int(current);

	// if in FAR_SET move to TRIAL_SET
	if ( this->state( current ) == FAR_SET ){
		this->state( current ) = TRIAL_SET;
	}
	// if already in TRIAL_SET, remove from queue
	else if ( this->state( current ) == TRIAL_SET ){
		float p = this->pqueue->DeleteId( currentId );
		if ( p == VTK_LARGE_FLOAT ){
			cout << "point not in queue" << endl;
		}
	}

	// insert in queue
	this->pqueue->Insert( this->distance(current), currentId );
}

template< class Pixel, int const Dim, class Scalar >
TinyVector< int, Dim > nbfFastMarching< Pixel, Dim, Scalar > :: execute( Array< Scalar, Dim > & output )
{
	TinyVector< int, Dim > currentPosition;

	// if initial alive set not empty
	if( this->aliveSetPositions.size() != 0 )
	{
		// initialize state and distances
		this->distance.reference( output );
		this->initialize();

		// freezing: reset furthest point
		this->dmax = 0;

		// stoping criteria
		int stopPosicion;
		bool stopDistancia;
		bool stopBorde;

		do {
			// get point in TRIAL_SET with the smallest distance value
			double vdistance;
			int id = this->pqueue->Pop(0,vdistance);

			// get coordinates from queue index
			currentPosition = this->int2array(id);

			//cout << currentPosition << ", d = " << vdistance << endl;

			// move it to the ACCEPTED_SET
			this->state( currentPosition ) = ACCEPTED_SET;

			// update inside loop
			this->loopUpdate( currentPosition );

			//this->flagForceStop = this->checkForceStop( currentPosition );

			// if freeze
			if ( this->dth && ( this->times( currentPosition ) < this->dmax - this->dth ) ){
				this->times( currentPosition ) = numeric_limits<Scalar>::max();
				this->distance( currentPosition ) = numeric_limits<Scalar>::max();
			}
			else{
				// only if using freezing
				if ( this->dth != 0 ){
					// update furthest point
					this->dmax = blitz::minmax::max( this->times( currentPosition ), this->dmax );			
				}

				// update neighbours
				this->updateNeighbors( currentPosition );
			}

			// check stop conditions
			stopPosicion = sum( abs( currentPosition - PuntoParada ) );
			if ( !this->distance.isInRange( currentPosition ) ){
				cout << "fastmarching error " << currentPosition << endl;
			}
			stopDistancia = ( this->distance( currentPosition ) - DistanciaParada ) < 0;
			if ( this->hasCut ){
				stopBorde = !( ( currentPosition(firstDim) > this->center[firstDim] ) & ( currentPosition(secondDim) == ( this->center[secondDim] + 1 ) ) );
			}
			else{
				stopBorde = ( currentPosition(BordeParada) < output.ubound( BordeParada ) );
			}
			//stopBorde = ( pd.position(BordeParada) > output.lbound( BordeParada ) );

			//if ( sum( where( this->state == ACCEPTED_SET,1,0) ) == 3500 ){
			//	break;
			//}

		} while ( ( this->pqueue->GetNumberOfItems() != 0 ) &&
			( !flagForceStop ) &&
			( ( stopPosicion != 0 ) || !( flagPuntoParada ) ) &&
			( ( stopDistancia != 0 ) || !( flagDistanciaParada ) ) &&
			( ( stopBorde != 0 ) || !( flagBordeParada ) ) );
	}
	else{
		output = numeric_limits<Scalar>::max();
	}

	return currentPosition;
}

#endif /* FILE_nbfFastMarching */