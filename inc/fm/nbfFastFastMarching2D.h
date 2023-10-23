#pragma once

/** @file nbfFastFastMarching.h
*	2D optimized Fast Marching. Part of FM.
*/

#include <nbfArray.h>

#include <vector>
#include <queue>

#include <fm/DiscritePriorityQueue.h>

#undef QUEUE_LIRON

template< class Pixel >
class FastNode{
public:
	int x, y;
	Pixel d;
	FastNode(int ax,int ay,Pixel ad): x(ax), y(ay), d(ad){};
	FastNode(){};
};

template< class Pixel >
struct fastGreaterPD : public binary_function< FastNode< Pixel >, FastNode< Pixel >, bool >
{
  inline bool operator()( const FastNode< Pixel > & x, const FastNode< Pixel > & y) const 
    {
      return x.d > y.d;
    }
};

//template< class T, class Cmp, class C = std::vector<T> >
//class PriorityQueue : public std::priority_queue< T,C,Cmp>
//{
//public :
//  void clear(){ c.clear();};
//  explicit PriorityQueue() : std::priority_queue<T,C,Cmp>(){};
//  
//};

/** Optimized 2D fast marching.

	Implements 2D fast fast marching method with arbitrary weights.
*/
template< class Pixel >
class nbfFastFastMarching2D
{
public:

	static const int ACCEPTED_SET = 0, FAR_SET = 1, TRIAL_SET = 2;

	// constructor takes weight array as input
	nbfFastFastMarching2D( Array< Pixel, 2 > &, int = 512 );

	~nbfFastFastMarching2D();

	void construct( Array< Pixel, 2 > & );

	// fast marching initalization
	void inicialize();

	// set individual point and corresponding distance as initial set
	void setAliveSet( TinyVector< int, 2 >, Pixel );

	// set multiple points and corresponding distances as initial set
	void setAliveSet( vector< TinyVector< int, 2 > > &,  vector< Pixel > & = 0 );

	// Stopping conditions

	// set stoping point coordinates
	void setStopPoint( TinyVector< int, 2 > & );

	// unset stoping point flag
	void unSetStopPoint();

	// set stoping distance
	void setStopDistance( Pixel );

	// unset stoping distance flag
	void unSetStopDistance();

	// run fast marching, write result on argument
	virtual TinyVector< int, 2 > execute( Array< Pixel, 2 > & );

protected:

	// update point neighbors
	void updateNeighbors();

	// update individual point
	void updatePointXp();
	void updatePointXm();
	void updatePointYp();
	void updatePointYm();

	DiscritePriorityQueue< FastNode< Pixel > > queueLiron;

	PriorityQueue< FastNode< Pixel >, fastGreaterPD< Pixel > , vector< FastNode< Pixel > > > pqNB;

	// data dimensions
	TinyVector< int, 2 > dimensions;
	int maxX, maxY;

	int currentX, currentY;
	int currentAliveX, currentAliveY;

	Pixel currentAliveDistance;

	// point stoping
	bool flagPuntoParada;
	TinyVector< int, 2 > PuntoParada;

	// distance stoping
	bool flagDistanciaParada;
	Pixel DistanciaParada;

	// initial alive set coordinates and distances
	vector< TinyVector< int, 2 > > aliveSetPositions;
	vector< Pixel > aliveSetDistances;

	// distances (taken from execute())
	Array< Pixel, 2 > distance;
	Pixel * pDistance;

	int rows, cols;

	int * pState;
	
	// weights (taken from constructor)
	Array< Pixel, 2 > weight;
	Pixel * pWeight;

};

template< class Pixel >
nbfFastFastMarching2D< Pixel > :: nbfFastFastMarching2D( Array< Pixel, 2 > & weight, int queueSize )
: queueLiron( queueSize, max(weight) )
{
	this->construct( weight );
}

template< class Pixel >
nbfFastFastMarching2D< Pixel > :: ~nbfFastFastMarching2D()
{
	delete [] this->pState;
}

template< class Pixel >
void nbfFastFastMarching2D< Pixel > :: construct( Array< Pixel, 2 > & weight )
{
	this->weight.resize( weight.shape() );
	this->weight = weight;
	this->pWeight = this->weight.dataZero();

	this->dimensions = weight.shape();

	this->rows = weight.rows();
	this->cols = weight.cols();
	this->maxX = this->rows - 1;
	this->maxY = this->cols - 1;

	this->pState = new int[ this->rows * this->cols ];

	this->PuntoParada = 0;
	this->flagPuntoParada = false;

	this->DistanciaParada = numeric_limits<Pixel>::max();
	this->flagDistanciaParada = false;
}

template< class Pixel >
void nbfFastFastMarching2D< Pixel > :: setAliveSet( TinyVector< int, 2 > position,
											        Pixel distance )
{
	vector< TinyVector< int, 2 > > vPositions;
	vector< Pixel > vDistances;

	vPositions.push_back( position );
	vDistances.push_back( distance );

	this->setAliveSet( vPositions, vDistances );
}


template< class Pixel >
void nbfFastFastMarching2D< Pixel > :: setAliveSet( vector< TinyVector< int, 2 > > & positions,
											        vector< Pixel > & distances )
{
	// reset attributes

	this->aliveSetPositions.clear();
	this->aliveSetDistances.clear();

	// deep copy both arguments into attributes

	vector< TinyVector< int, 2 > > :: iterator iterPositions = positions.begin();
	vector< Pixel > :: iterator iterDistances = distances.begin();

	while( iterPositions != positions.end() ){
		if ( weight.isInRange( *iterPositions ) ){
			this->aliveSetPositions.push_back( *iterPositions );
			this->aliveSetDistances.push_back( *iterDistances );
		}
		++iterPositions; 
		++iterDistances;
	}
}


template< class Pixel >
void nbfFastFastMarching2D< Pixel > :: inicialize()
{
#ifndef QUEUE_LIRON
	this->pqNB.clear();
#endif

	// set all distances to infinity and all states to FAR
	int size = this->rows * this->cols;
	for ( int i = 0; i < size; i++ ){
		this->pDistance[ i ] = numeric_limits<Pixel>::max();
		this->pState[ i ] = FAR_SET;
	}

	// initialize all points in initial alive set
	vector< TinyVector< int, 2 > > :: iterator iterPositions;
	vector< Pixel > :: iterator iterDistances;

	for ( iterPositions = this->aliveSetPositions.begin(), 
		iterDistances = this->aliveSetDistances.begin(); 
		iterPositions != this->aliveSetPositions.end(),
		iterDistances != this->aliveSetDistances.end();
		++iterPositions, ++iterDistances ){

		// set distance values
		this->pDistance[ (*iterPositions)[firstDim] * this->cols + (*iterPositions)[secondDim] ] = *iterDistances;

		FastNode< Pixel > next( (*iterPositions)[firstDim], (*iterPositions)[secondDim], (*iterDistances) );
#ifdef QUEUE_LIRON
		this->queueLiron.insertSeed( (*iterDistances), next );
#else
		this->pqNB.push( next );
#endif
	}
}

template< class Pixel >
void nbfFastFastMarching2D< Pixel > :: updateNeighbors()
{
	// if valid insert into neigbors list
	int delta = this->currentAliveX * this->cols + this->currentAliveY;

	int offset;

	if ( this->currentAliveX > 0 ){
		offset = delta - this->cols;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentX = this->currentAliveX - 1;
			this->updatePointXm();		
		}
	}
	if ( this->currentAliveX < this->maxX ){
		offset = delta + this->cols;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentX = this->currentAliveX + 1;
			this->updatePointXp();		
		}
	}
	if ( this->currentAliveY > 0 ){
		offset = delta - 1;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentY = this->currentAliveY - 1;
			this->updatePointYm();		
		}
	}
	if ( this->currentAliveY < this->maxY ){
		offset = delta + 1;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentY = this->currentAliveY + 1;
			this->updatePointYp();		
		}
	}
}

template< class Pixel >
TinyVector< int, 2 > nbfFastFastMarching2D< Pixel > :: execute( Array< Pixel, 2 > & output )
{
	// if initial alive set not empty
	if( this->aliveSetPositions.size() != 0 )
	{
		// initialize state and distances
		this->distance.resize( this->weight.shape() );

		this->pDistance = this->distance.dataZero();

		this->inicialize();

		FastNode< Pixel > firstOut;

		int stopPosicion;
		bool stopDistancia;

		do {
			// get point in TRIAL with the smallest distance value
			
#ifdef QUEUE_LIRON
			do {
				this->queueLiron.pop(firstOut);
			} while( ( this->pState[ firstOut.x * this->cols + firstOut.y ] == ACCEPTED_SET ) && ( !this->queueLiron.isEmpty() ) );
#else
			do {
				firstOut = pqNB.top();
				pqNB.pop();
			} while( ( this->pState[ firstOut.x * this->cols + firstOut.y ] == ACCEPTED_SET ) & ( pqNB.size() != 0 ) );
#endif

			this->currentAliveX = firstOut.x;
			this->currentAliveY = firstOut.y;

			int delta = this->currentAliveX * this->cols + this->currentAliveY;

			this->currentAliveDistance = this->pDistance[ delta ];

			this->pState[ delta ] = ACCEPTED_SET;

			// update neighbours
			this->updateNeighbors();

			stopPosicion = sum( abs( TinyVector<int,2>(this->currentAliveX,this->currentAliveY) - PuntoParada ) );
			stopDistancia = ( this->currentAliveDistance - DistanciaParada ) < 0;

#ifdef QUEUE_LIRON
		} while ( ( !this->queueLiron.isEmpty() ) &&
			      ( ( stopPosicion != 0 ) || !( flagPuntoParada ) ) &&
			      ( ( stopDistancia != 0 ) || !( flagDistanciaParada ) ) );
#else
		} while ( ( this->pqNB.size() != 0 ) &&
			      ( ( stopPosicion != 0 ) || !( flagPuntoParada ) ) &&
			      ( ( stopDistancia != 0 ) || !( flagDistanciaParada ) ) );
#endif
	}
	else{
		//output = numeric_limits<Pixel>::max();
	}
	//output = this->distance( Range(1,this->distance.ubound(firstDim)-1), Range(1,this->distance.ubound(secondDim)-1) );
	output.resize( this->distance.shape() );
	output = this->distance;
	return TinyVector<int,2>(this->currentAliveX,this->currentAliveY);
}

template< class Pixel >
void nbfFastFastMarching2D< Pixel > :: updatePointXp()
{
	Pixel uB2 = numeric_limits<Pixel>::max();

	int delta = this->currentX * this->cols + this->currentAliveY;

	if ( this->currentX < this->maxX ){
		uB2 = this->pDistance[ delta + this->cols ];
	}

	if ( this->currentAliveDistance <= uB2 ){

		Pixel uA1 = numeric_limits<Pixel>::max();
		Pixel uA2 = numeric_limits<Pixel>::max();

		if ( this->currentAliveY < this->maxY ){
			uA1 = this->pDistance[ delta + 1 ];
		}
		if ( this->currentAliveY > 0 ){
			uA2 = this->pDistance[ delta - 1];
		}

		// Pixel uxj = min( uA1, uA2 );
		Pixel uxj = (uA1<uA2)?uA1:uA2;

		// retrieve weight
		Pixel tx = this->pWeight[ delta ];

		// tentative distance value
		Pixel uxjxm = numeric_limits<Pixel>::max();

		if ( uxj <= this->currentAliveDistance ){
			uxjxm = ( uxj + this->currentAliveDistance + sqrt( 2.0 * pow2(tx) - pow2( uxj - this->currentAliveDistance ) ) ) / 2.0;
		}
		else{
			uxjxm = this->currentAliveDistance + tx;
		}

		Pixel * currentDistance = &this->pDistance[ delta ];

		// update distance value
		if ( uxjxm < (*currentDistance) ){
			(*currentDistance) = uxjxm;
			this->pState[ delta ] = TRIAL_SET;
			FastNode< Pixel > next( this->currentX, this->currentAliveY, (*currentDistance) );
#ifdef QUEUE_LIRON
			this->queueLiron.insert( (*currentDistance), next );
#else
			this->pqNB.push( next );
#endif
		}
	}
}

template< class Pixel >
void nbfFastFastMarching2D< Pixel > :: updatePointXm()
{
	Pixel uB2 = numeric_limits<Pixel>::max();

	int delta = this->currentX * this->cols + this->currentAliveY;

	if ( this->currentX > 0 ){
		uB2 = this->pDistance[ delta - this->cols ];
	}

	if ( this->currentAliveDistance <= uB2 ){

		Pixel uA1 = numeric_limits<Pixel>::max();
		Pixel uA2 = numeric_limits<Pixel>::max();

		if ( this->currentAliveY < this->maxY ){
			uA1 = this->pDistance[ delta + 1 ];
		}
		if ( this->currentAliveY > 0 ){
			uA2 = this->pDistance[ delta - 1 ];
		}

		// Pixel uxj = min( uA1, uA2 );
		Pixel uxj = (uA1<uA2)?uA1:uA2;

		// retrieve weight
		Pixel tx = this->pWeight[ delta ];

		// tentative distance value
		Pixel uxjxm = numeric_limits<Pixel>::max();

		if ( uxj <= this->currentAliveDistance ){
			uxjxm = ( uxj + this->currentAliveDistance + sqrt( 2.0 * pow2(tx) - pow2( uxj - this->currentAliveDistance ) ) ) / 2.0;
		}
		else{
			uxjxm = this->currentAliveDistance + tx;
		}

		Pixel * currentDistance = &this->pDistance[ delta ];

		// update distance value
		if ( uxjxm < (*currentDistance) ){
			(*currentDistance) = uxjxm;

			this->pState[ delta ] = TRIAL_SET;

			FastNode< Pixel > next( this->currentX, this->currentAliveY, (*currentDistance) );
#ifdef QUEUE_LIRON
			this->queueLiron.insert( (*currentDistance), next );
#else
			this->pqNB.push( next );
#endif
		}
	}
}

template< class Pixel >
void nbfFastFastMarching2D< Pixel > :: updatePointYp()
{
	Pixel uB2 = numeric_limits<Pixel>::max();

	int delta = this->currentAliveX * this->cols + this->currentY;

	if ( this->currentY < this->maxY ){
		uB2 = this->pDistance[ delta + 1 ];
	}

	if ( this->currentAliveDistance <= uB2 ){

		Pixel uA1 = numeric_limits<Pixel>::max();
		Pixel uA2 = numeric_limits<Pixel>::max();

		if ( this->currentAliveX > 0 ){
			uA1 = this->pDistance[ delta - this->cols ];
		}
		if ( this->currentAliveX < this->maxX ){
			uA2 = this->pDistance[ delta + this->cols ];
		}

		// Pixel uxj = min( uA1, uA2 );
		Pixel uxj = (uA1<uA2)?uA1:uA2;

		// retrieve weight
		Pixel tx = this->pWeight[ delta ];
	
		// tentative distance value
		Pixel uxjxm = numeric_limits<Pixel>::max();

		if ( uxj <= this->currentAliveDistance ){
			uxjxm = ( uxj + this->currentAliveDistance + sqrt( 2.0 * pow2(tx) - pow2( uxj - this->currentAliveDistance ) ) ) / 2.0;
		}
		else{
			uxjxm = this->currentAliveDistance + tx;
		}

		Pixel * currentDistance = &this->pDistance[ delta ];

		// update distance value
		if ( uxjxm < (*currentDistance) ){

			(*currentDistance) = uxjxm;
			this->pState[ delta ] = TRIAL_SET;

			FastNode< Pixel > next( this->currentAliveX, this->currentY, (*currentDistance) );
#ifdef QUEUE_LIRON
			this->queueLiron.insert( (*currentDistance), next );
#else
			this->pqNB.push( next );
#endif
		}
	}
}

template< class Pixel >
void nbfFastFastMarching2D< Pixel > :: updatePointYm()
{
	Pixel uB2 = numeric_limits<Pixel>::max();

	int delta = this->currentAliveX * this->cols + this->currentY;

	if ( this->currentY > 0 ){
		uB2 = this->pDistance[ delta - 1  ];
	}

	if ( this->currentAliveDistance <= uB2 ){

		Pixel uA1 = numeric_limits<Pixel>::max();
		Pixel uA2 = numeric_limits<Pixel>::max();

		if ( this->currentAliveX > 0 ){
			uA1 = this->pDistance[ delta - this->cols ];
		}
		if ( this->currentAliveX < this->maxX ){
			uA2 = this->pDistance[ delta + this->cols ];
		}

		// Pixel uxj = min( uA1, uA2 );
		Pixel uxj = (uA1<uA2)?uA1:uA2;

		// retrieve weight
		Pixel tx = this->pWeight[ delta ];
	
		// tentative distance value
		Pixel uxjxm = numeric_limits<Pixel>::max();

		if ( uxj <= this->currentAliveDistance ){
			uxjxm = ( uxj + this->currentAliveDistance + sqrt( 2.0 * pow2(tx) - pow2( uxj - this->currentAliveDistance ) ) ) / 2.0;
		}
		else{
			uxjxm = this->currentAliveDistance + tx;
		}

		Pixel * currentDistance = &this->pDistance[ delta ];

		// update distance value
		if ( uxjxm < (*currentDistance) ){
			(*currentDistance) = uxjxm;

			this->pState[ delta ] = TRIAL_SET;

			FastNode< Pixel > next( this->currentAliveX, this->currentY, (*currentDistance) );
#ifdef QUEUE_LIRON
			this->queueLiron.insert( (*currentDistance), next );
#else
			this->pqNB.push( next );
#endif
		}
	}
}

template< class Pixel >
void nbfFastFastMarching2D< Pixel > :: setStopPoint( TinyVector< int, 2 > & pp )
{
	flagDistanciaParada = false;

	PuntoParada = pp;
	flagPuntoParada = true;
}

template< class Pixel >
void nbfFastFastMarching2D< Pixel > :: setStopDistance( Pixel dp )
{
	flagPuntoParada = false;

	DistanciaParada = dp;
	flagDistanciaParada = true;
}