#pragma once

/** @file nbfFastFastMarching3D.h
*	3D optimized Fast Marching. Part of FM.
*/

#include <nbfArray.h>

#include <vector>
#include <queue>

#include <fm/DiscritePriorityQueue.h>

#undef QUEUE_LIRON

template< class Pixel >
class FastNode{
public:
	int x, y, z;
	Pixel d;
	FastNode(int ax,int ay, int az, Pixel ad ): x(ax), y(ay), z(az), d(ad){};
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

template< class T, class Cmp, class C = std::vector<T> >
class PriorityQueue : public std::priority_queue< T,C,Cmp>
{
public :
  void clear(){ c.clear();};
  explicit PriorityQueue() : std::priority_queue<T,C,Cmp>(){};
  
};

/** Optimized 3D fast marching.

	Implements 3D fast fast marching method with arbitrary weights.
*/
template< class Pixel >
class nbfFastFastMarching3D
{
public:

	static const int ACCEPTED_SET = 0, FAR_SET = 1, TRIAL_SET = 2;

	// constructor takes weight array as input
	nbfFastFastMarching3D( Array< Pixel, 3 > &, int = 512 );

	~nbfFastFastMarching3D();

	void construct( Array< Pixel, 3 > & );

	// fast marching initalization
	void inicialize();

	// set individual point and corresponding distance as initial set
	void setAliveSet( TinyVector< int, 3 >, Pixel );

	// set multiple points and corresponding distances as initial set
	void setAliveSet( vector< TinyVector< int, 3 > > &,  vector< Pixel > & = 0 );

	// Stopping conditions

	// set stoping point coordinates
	void setStopPoint( TinyVector< int, 3 > & );

	// unset stoping point flag
	void unSetStopPoint();

	// set stoping distance
	void setStopDistance( Pixel );

	// unset stoping distance flag
	void unSetStopDistance();

	// run fast marching, write result on argument
	virtual TinyVector< int, 3 > execute( Array< Pixel, 3 > & );

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
	TinyVector< int, 3 > dimensions;
	int maxX, maxY, maxZ;

	int currentX, currentY, currentZ;
	int currentAliveX, currentAliveY, currentAliveZ;

	Pixel currentAliveDistance;

	// point stoping
	bool flagPuntoParada;
	TinyVector< int, 3 > PuntoParada;

	// distance stoping
	bool flagDistanciaParada;
	Pixel DistanciaParada;

	// initial alive set coordinates and distances
	vector< TinyVector< int, 3 > > aliveSetPositions;
	vector< Pixel > aliveSetDistances;

	// distances (taken from execute())
	Array< Pixel, 3 > distance;
	Pixel * pDistance;

	int rows, cols, depth;

	int * pState;
	
	// weights (taken from constructor)
	Array< Pixel, 3 > weight;
	Pixel * pWeight;

};

template< class Pixel >
nbfFastFastMarching3D< Pixel > :: nbfFastFastMarching3D( Array< Pixel, 3 > & weight, int queueSize )
: queueLiron( queueSize, max(weight) )
{
	this->construct( weight );
}

template< class Pixel >
nbfFastFastMarching3D< Pixel > :: ~nbfFastFastMarching3D()
{
	delete [] this->pState;
}

template< class Pixel >
void nbfFastFastMarching3D< Pixel > :: construct( Array< Pixel, 3 > & weight )
{
	this->weight.resize( weight.shape() );
	this->weight = weight;
	this->pWeight = this->weight.dataZero();

	this->dimensions = weight.shape();

	this->rows = weight.rows();
	this->cols = weight.cols();
	this->depth = weight.depth();
	this->maxX = this->rows - 1;
	this->maxY = this->cols - 1;
	this->maxZ = this->depth - 1;

	this->pState = new int[ this->rows * this->cols * this->depth ];

	this->PuntoParada = 0;
	this->flagPuntoParada = false;

	this->DistanciaParada = numeric_limits<Pixel>::max();
	this->flagDistanciaParada = false;
}

template< class Pixel >
void nbfFastFastMarching3D< Pixel > :: setAliveSet( TinyVector< int, 3 > position,
											        Pixel distance )
{
	vector< TinyVector< int, 3 > > vPositions;
	vector< Pixel > vDistances;

	vPositions.push_back( position );
	vDistances.push_back( distance );

	this->setAliveSet( vPositions, vDistances );
}


template< class Pixel >
void nbfFastFastMarching3D< Pixel > :: setAliveSet( vector< TinyVector< int, 3 > > & positions,
											        vector< Pixel > & distances )
{
	// reset attributes

	this->aliveSetPositions.clear();
	this->aliveSetDistances.clear();

	// deep copy both arguments into attributes

	vector< TinyVector< int, 3 > > :: iterator iterPositions = positions.begin();
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
void nbfFastFastMarching3D< Pixel > :: inicialize()
{
#ifndef QUEUE_LIRON
	this->pqNB.clear();
#endif

	// set all distances to infinity and all states to FAR
	int size = this->rows * this->cols * this->depth;
	for ( int i = 0; i < size; i++ ){
		this->pDistance[ i ] = numeric_limits<Pixel>::max();
		this->pState[ i ] = FAR_SET;
	}

	// initialize all points in initial alive set
	vector< TinyVector< int, 3 > > :: iterator iterPositions;
	vector< Pixel > :: iterator iterDistances;

	for ( iterPositions = this->aliveSetPositions.begin(), 
		iterDistances = this->aliveSetDistances.begin(); 
		iterPositions != this->aliveSetPositions.end(),
		iterDistances != this->aliveSetDistances.end();
		++iterPositions, ++iterDistances ){

		// set distance values
		this->pDistance[ ( (*iterPositions)[firstDim] * this->depth + (*iterPositions)[secondDim] ) * this->cols + (*iterPositions)[thirdDim] ] = *iterDistances;

		FastNode< Pixel > next( (*iterPositions)[firstDim], (*iterPositions)[secondDim], (*iterPositions)[thirdDim], (*iterDistances) );
#ifdef QUEUE_LIRON
		this->queueLiron.insertSeed( (*iterDistances), next );
#else
		this->pqNB.push( next );
#endif
	}
}

template< class Pixel >
void nbfFastFastMarching3D< Pixel > :: updateNeighbors()
{
	// if valid insert into neigbors list
	int delta = ( this->currentAliveX * this->depth + this->currentAliveY ) * this->cols;

	int offset;

	if ( this->currentAliveX > 0 ){
		offset = delta - this->cols * this->depth;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentX = this->currentAliveX - 1;
			this->updatePointXm();		
		}
	}
	if ( this->currentAliveX < this->maxX ){
		offset = delta + this->cols * this->depth;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentX = this->currentAliveX + 1;
			this->updatePointXp();		
		}
	}
	if ( this->currentAliveY > 0 ){
		offset = delta - this->cols;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentY = this->currentAliveY - 1;
			this->updatePointYm();		
		}
	}
	if ( this->currentAliveY < this->maxY ){
		offset = delta + this->cols;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentY = this->currentAliveY + 1;
			this->updatePointYp();		
		}
	}
	if ( this->currentAliveZ > 0 ){
		offset = delta - 1;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentZ = this->currentAliveZ - 1;
			this->updatePointZm();		
		}
	}
	if ( this->currentAliveZ < this->maxZ ){
		offset = delta + 1;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentZ = this->currentAliveZ + 1;
			this->updatePointZp();		
		}
	}
}

template< class Pixel >
TinyVector< int, 3 > nbfFastFastMarching3D< Pixel > :: execute( Array< Pixel, 3 > & output )
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
			} while( ( this->pState[ ( firstOut.x * this->depth + firstOut.y ) * this->cols + firstOut.z ] == ACCEPTED_SET ) && ( !this->queueLiron.isEmpty() ) );
#else
			do {
				firstOut = pqNB.top();
				pqNB.pop();
			} while( ( this->pState[ ( firstOut.x * this->depth + firstOut.y ) * this->cols + firstOut.z ] == ACCEPTED_SET ) & ( pqNB.size() != 0 ) );
#endif

			this->currentAliveX = firstOut.x;
			this->currentAliveY = firstOut.y;
			this->currentAliveZ = firstOut.z;

			int delta = ( this->currentAliveX * this->depth + this->currentAliveY ) * this->cols + this->currentAliveZ;

			this->currentAliveDistance = this->pDistance[ delta ];

			this->pState[ delta ] = ACCEPTED_SET;

			// update neighbours
			this->updateNeighbors();

			stopPosicion = sum( abs( TinyVector<int,3>(this->currentAliveX,this->currentAliveY,this->currentAliveZ) - PuntoParada ) );
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
	return TinyVector<int,3>(this->currentAliveX,this->currentAliveY,this->currentAliveZ);
}

template< class Pixel >
void nbfFastFastMarching3D< Pixel > :: updateNeighbors()
{
	int delta = ( this->currentAliveX * this->depth + this->currentAliveY ) * this->cols + this->currentAliveZ;

	int offset;

	// 1 of 6
	if ( this->currentAliveX > 0 ){
		offset = delta - this->cols * this->depth;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentX = this->currentAliveX - 1;
			this->updatePointXm();		
		}
	}
	// 2 of 6
	if ( this->currentAliveX < this->maxX ){
		offset = delta + this->cols * this->depth;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentX = this->currentAliveX + 1;
			this->updatePointXp();		
		}
	}
	// 3 of 6
	if ( this->currentAliveY > 0 ){
		offset = delta - this->cols;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentY = this->currentAliveY - 1;
			this->updatePointYm();		
		}
	}
	// 4 of 6
	if ( this->currentAliveY < this->maxY ){
		offset = delta + this->cols;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentY = this->currentAliveY + 1;
			this->updatePointYp();		
		}
	}
	// 5 of 6
	if ( this->currentAliveZ > 0 ){
		offset = delta - 1;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentZ = this->currentAliveZ - 1;
			this->updatePointZm();		
		}
	}
	// 6 of 6
	if ( this->currentAliveZ < this->maxZ ){
		offset = delta + 1;
		if ( ( this->pState[ offset ] != ACCEPTED_SET ) && ( this->pWeight[ offset ] < numeric_limits<Pixel>::max() ) ){
			this->currentZ = this->currentAliveZ + 1;
			this->updatePointZp();		
		}
	}
}

template< class Pixel >
void nbfFastFastMarching3D< Pixel > :: updatePointXp()
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
		Pixel uxjxm = solve3D( uA, uB, uC );

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
			FastNode< Pixel > next( this->currentX, this->currentAliveY, this->currentAliveZ, (*currentDistance) );
#ifdef QUEUE_LIRON
			this->queueLiron.insert( (*currentDistance), next );
#else
			this->pqNB.push( next );
#endif
		}
	}
}

template< class Pixel >
void nbfFastFastMarching3D< Pixel > :: updatePoint( TinyVector< int, 3 > & current, TinyVector< int, 3 > & alive )
{
	// coordinates of point we are updating
	int x = current(firstDim);
	int y = current(secondDim);
	int z = current(thirdDim);

	// 6-connected neighbors
	vector< TinyVector< int, 3 > > neighbors;
	TinyVector< int, 3 > offset;

	// find the fixed dimension
	if ( current(firstDim) != alive(firstDim) ){
		offset = current;
		offset(secondDim) = current(secondDim) + 1;
		if ( this->distance.isInRange( offset ) ){
			neighbors.push_back( offset );
		}

		offset = current;
		offset(thirdDim) = current(thirdDim) + 1;
		if ( this->distance.isInRange( offset ) ){
			neighbors.push_back( offset );
		}

		offset = current;
		offset(secondDim) = current(secondDim) - 1;
		if ( this->distance.isInRange( offset ) ){
			neighbors.push_back( offset );
		}

		offset = current;
		offset(thirdDim) = current(thirdDim) - 1;
		if ( this->distance.isInRange( offset ) ){
			neighbors.push_back( offset );
		}

		offset = current;
		offset(secondDim) = current(secondDim) + 1;
		if ( this->distance.isInRange( offset ) ){
			neighbors.push_back( offset );
		}
	}
	else{
		if ( current(secondDim) != alive(secondDim) ){
			offset = current;
			offset(firstDim) = current(firstDim) + 1;
			if ( this->distance.isInRange( offset ) ){
				neighbors.push_back( offset );
			}

			offset = current;
			offset(thirdDim) = current(thirdDim) + 1;
			if ( this->distance.isInRange( offset ) ){
				neighbors.push_back( offset );
			}

			offset = current;
			offset(firstDim) = current(firstDim) - 1;
			if ( this->distance.isInRange( offset ) ){
				neighbors.push_back( offset );
			}

			offset = current;
			offset(thirdDim) = current(thirdDim) - 1;
			if ( this->distance.isInRange( offset ) ){
				neighbors.push_back( offset );
			}

			offset = current;
			offset(firstDim) = current(firstDim) + 1;
			if ( this->distance.isInRange( offset ) ){
				neighbors.push_back( offset );
			}
		}
		else{
			offset = current;
			offset(secondDim) = current(secondDim) + 1;
			if ( this->distance.isInRange( offset ) ){
				neighbors.push_back( offset );
			}

			offset = current;
			offset(firstDim) = current(firstDim) + 1;
			if ( this->distance.isInRange( offset ) ){
				neighbors.push_back( offset );
			}

			offset = current;
			offset(secondDim) = current(secondDim) - 1;
			if ( this->distance.isInRange( offset ) ){
				neighbors.push_back( offset );
			}

			offset = current;
			offset(firstDim) = current(firstDim) - 1;
			if ( this->distance.isInRange( offset ) ){
				neighbors.push_back( offset );
			}

			offset = current;
			offset(secondDim) = current(secondDim) + 1;
			if ( this->distance.isInRange( offset ) ){
				neighbors.push_back( offset );
			}
		}
	}

	// tentative distance value
	Scalar uxm;

	// current distance value
	Scalar * currentD;

	// retrieve weight
	Pixel t = this->weight( current );

	// for each 6-connected neighbor
	for ( int i = 0; i < neighbors.size() - 1; i++ ){
		uxm = this->solve3D( neighbors[i], neighbors[i+1], alive, t );
	
		// update distance value
		//currentD = &( this->distance( current ) );
		//(*currentD) = min( *currentD, uxm );
		this->distance( current ) = min( this->distance( current ), uxm );
	}

	this->updateQueue( current );
}


template< class Pixel >
Scalar nbfFastMarching3D< Pixel > :: solve3D( TinyVector< int, 3 > & A,
											  TinyVector< int, 3 > & B,
										      TinyVector< int, 3 > & C,
										      Pixel & t )
{
	Scalar u, uA, uB, uC;

	uA = this->getDistance(A);
	uB = this->getDistance(B);
	uC = this->getDistance(C);

	int sA = this->state( A );
	int sB = this->state( B );

	if ( ( sA == ACCEPTED_SET ) && ( sB == ACCEPTED_SET ) ){
		u = ( uA + uB + uC + sqrt( 3 * ( t * t - uA * uA - uB * uB - uC * uC ) + ( uA + uB + uC ) * ( uA + uB + uC ) ) ) / 3.0;
	}
	else{
		if ( sA == ACCEPTED_SET ){
			u = ( uA + uC + sqrt( 2 * t * t - ( uA - uC ) * ( uA - uC ) ) ) / 2.0;
		}
		else{
			if ( sB == ACCEPTED_SET ){
				u = ( uB + uC + sqrt( 2 * t * t - ( uB - uC ) * ( uB - uC ) ) ) / 2.0;
			}
			else{
				u = uC + t;
			}
		}
	}
	return u;
}


template< class Pixel >
void nbfFastFastMarching3D< Pixel > :: setStopPoint( TinyVector< int, 3 > & pp )
{
	flagDistanciaParada = false;

	PuntoParada = pp;
	flagPuntoParada = true;
}

template< class Pixel >
void nbfFastFastMarching3D< Pixel > :: setStopDistance( Pixel dp )
{
	flagPuntoParada = false;

	DistanciaParada = dp;
	flagDistanciaParada = true;
}