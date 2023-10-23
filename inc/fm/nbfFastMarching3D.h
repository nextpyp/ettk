#ifndef FILE_nbfFastMarching3D
#define FILE_nbfFastMarching3D

// Class nbfFastMarching3D.
//
// Implements 3D fast marching method with arbitrary weights.
// 

#include <fm/nbfFastMarching.h>

template< class Pixel >
class nbfFastMarching3D : public nbfFastMarching< Pixel, 3 >
{
public:

	// constructor takes weight array as input
	nbfFastMarching3D( Array< Pixel, 3 > & );
	nbfFastMarching3D( Array< Pixel, 3 > &, TinyVector< int, 3 > & );

protected:

	// update individual point
	// the first argument is the point to update
	// the second argument is the point from where we do the update
	void updatePoint( TinyVector< int, 3 > &, TinyVector< int, 3 > & );

	Pixel solve3D( TinyVector< int, 3 > &, TinyVector< int, 3 > &, Pixel &, Pixel & );

	/** Handle weight computation on demand.
		When there is no need to compute distances in the entire domain
		and the weight function is computationally expensive to get, this
		method provides significantly computational savings, by only computing
		weight values as they are needed.
	*/
	virtual Pixel getWeight( TinyVector< int, 3 > & );

	// convert back and fort linear and matrix indexes
	int array2int( TinyVector< int, 3 > & );
	TinyVector< int, 3 > int2array( int );
};


template< class Pixel >
nbfFastMarching3D< Pixel > :: nbfFastMarching3D( Array< Pixel, 3 > & weight )
: nbfFastMarching< Pixel, 3 >( weight )
{
}

template< class Pixel >
nbfFastMarching3D< Pixel > :: nbfFastMarching3D( Array< Pixel, 3 > & weight, TinyVector< int, 3 > & center )
: nbfFastMarching< Pixel, 3 >( weight, center )
{
}

template< class Pixel >
void nbfFastMarching3D< Pixel > :: updatePoint( TinyVector< int, 3 > & current, TinyVector< int, 3 > & alive )
{
	// coordinates of point we are updating
	//int x = current(firstDim);
	//int y = current(secondDim);
	//int z = current(thirdDim);

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
	Pixel uxm;

	// current distance value
	Pixel * currentD = &( this->distance( current ) );

	// retrieve weight
	Pixel t = this->getWeight( current );

	Pixel uC = this->distance(alive);

	// for each 6-connected neighbor
	for ( int i = 0; i < neighbors.size() - 1; i++ ){
		uxm = this->solve3D( neighbors[i], neighbors[i+1], uC, t );
	
		// update distance value
		//currentD = &( this->distance( current ) );
		(*currentD) = min( (*currentD), uxm );
		//this->distance( current ) = min( this->distance( current ), uxm );
	}

	this->updateQueue( current );

	//// if in FAR_SET move to TRIAL_SET and push to queue
	//if ( this->state( current ) == FAR_SET ){
	//	this->state( current ) = TRIAL_SET;
	//	nbfDataPoint< Pixel, 3 > tmp( current );
	//	tmp.distance = &((*this->distance)(current));
	//	pqueue.push( tmp );
	//}
}


template< class Pixel >
Pixel nbfFastMarching3D< Pixel > :: solve3D( TinyVector< int, 3 > & A,
											TinyVector< int, 3 > & B,
											Pixel & uC,
											Pixel & t )
{
	Pixel u, uA, uB;

	uA = this->distance(A);
	uB = this->distance(B);
	//uC = this->distance(C);

	int sA = this->state( A );
	int sB = this->state( B );

	if ( ( sA == this->ACCEPTED_SET ) && ( sB == this->ACCEPTED_SET ) ){
		u = ( uA + uB + uC + sqrt( 3 * ( t*t - uA * uA - uB * uB - uC * uC ) + ( uA + uB + uC ) * ( uA + uB + uC ) ) ) / 3.0;
	}
	else{
		if ( sA == this->ACCEPTED_SET ){
			u = ( uA + uC + sqrt( 2 * t*t - ( uA - uC ) * ( uA - uC ) ) ) / 2.0;
		}
		else{
			if ( sB == this->ACCEPTED_SET ){
				u = ( uB + uC + sqrt( 2 * t*t - ( uB - uC ) * ( uB - uC ) ) ) / 2.0;
			}
			else{
				// WORKAROUND - to allow both the scalar and vector cases
				u = uC + t;
			}
		}
	}
	return u;
}

template< class Pixel >
Pixel nbfFastMarching3D< Pixel > :: getWeight( TinyVector< int, 3 > & p )
{
	return this->weight( p );
}

template< class Pixel >
int nbfFastMarching3D< Pixel > :: array2int( TinyVector< int, 3 > & t ){
	return t(2) * this->distance.rows() * this->distance.cols() +
		   t(1) * this->distance.rows() + t(0);
}

template< class Pixel >
TinyVector< int, 3 > nbfFastMarching3D< Pixel > :: int2array( int p ){
	Pixel k = floor( 1.0 * p / this->distance.rows() / this->distance.cols() );
	Pixel kmod = fmod( 1.0 * p, this->distance.rows() * this->distance.cols() );
	return TinyVector< int, 3 >( fmod( 1.0 * kmod, this->distance.rows() ),
								 floor( 1.0 * kmod / this->distance.rows() ),
								 k );
}
#endif /* FILE_nbfFastMarching3D */