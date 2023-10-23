#ifndef FILE_nbfFastMarching2D
#define FILE_nbfFastMarching2D

#include <fm/nbfFastMarching.h>

// Class nbfFastMarching2D.
//
// Implements 2D fast marching method with arbitrary weights.
// 

template< class Pixel >
class nbfFastMarching2D : public nbfFastMarching< Pixel, 2 >
{
public:

	// constructor takes weight array as input
	nbfFastMarching2D( Array< Pixel, 2 > & );
	nbfFastMarching2D( Array< Pixel, 2 > &, TinyVector< int, 2 > & );
	
	~nbfFastMarching2D(){};

protected:

	// update individual point
	// the first argument is the point to update
	// the second argument is the point from where we do the update
	void updatePoint( TinyVector< int, 2 > &, TinyVector< int, 2 > & );

	// 2D. get valid distance value for point (x,y)
	Pixel getDistance( Pixel, Pixel );

	// convert back and fort linear and matrix indexes
	int array2int( TinyVector< int, 2 > & );
	TinyVector< int, 2 > int2array( int );

};

template< class Pixel >
nbfFastMarching2D< Pixel > :: nbfFastMarching2D( Array< Pixel, 2 > & weight )
: nbfFastMarching< Pixel, 2 >( weight )
{
}

template< class Pixel >
nbfFastMarching2D< Pixel > :: nbfFastMarching2D( Array< Pixel, 2 > & weight, TinyVector< int, 2 > & center )
: nbfFastMarching< Pixel, 2 >( weight, center )
{
}

template< class Pixel >
Pixel nbfFastMarching2D< Pixel > :: getDistance( Pixel x, Pixel y )
{
	TinyVector< int, 2 > point(x,y);

	// if inside the array
	if ( this->distance.isInRange( point ) ){
		return ( this->distance( point ) );
	}
	// else infinity
	else{
		return numeric_limits<Pixel>::max();
	}
}

template< class Pixel >
void nbfFastMarching2D< Pixel > :: updatePoint( TinyVector< int, 2 > & current, TinyVector< int, 2 > & alive )
{
	// coordinates of point we are updating
	int x = current(firstDim);
	int y = current(secondDim);

	int dimension, offset;

	// distance of linked point
	Pixel uxm = this->distance( alive );

	if ( current(firstDim) == alive(firstDim) ){
		dimension = firstDim;
		offset = current(secondDim) - alive(secondDim);
	}
	else{
		dimension = secondDim;
		offset = current(firstDim) - alive(firstDim);
	}

	Pixel uB2, uA1, uA2;
	if ( dimension == firstDim ){
		uB2 = this->getDistance( x, y + offset );
		uA1 = this->getDistance( x - 1, y );
		uA2 = this->getDistance( x + 1, y );
	}
	else{
		uB2 = this->getDistance( x + offset, y );
		uA1 = this->getDistance( x, y + 1 );
		uA2 = this->getDistance( x, y - 1);
	}

	// retrieve weight
	Pixel tx = this->weight( current );
	
	// tentative distance value
	Pixel uxjxm;

	if ( uxm <= uB2 ){
		Pixel uxj = min( uA1, uA2 );
		if ( uxj <= uxm ){
			uxjxm = ( uxj + uxm + sqrt( 2.0 * tx * tx - ( uxj - uxm ) * ( uxj - uxm ) ) ) / 2.0;
		}
		else{
			uxjxm = uxm + tx;
		}
	}
	else{
		uxjxm = numeric_limits<Pixel>::max();
	}

	// update distance value
	this->distance( current ) = min( this->distance( current ), uxjxm );

	this->updateQueue( current );

	//// get global image ID
	//int currentId = this->array2int(current);

	//// if in FAR_SET move to TRIAL_SET
	//if ( this->state( current ) == FAR_SET ){
	//	this->state( current ) = TRIAL_SET;
	//}
	//// if already in TRIAL_SET, remove from queue
	//else if ( this->state( current ) == TRIAL_SET ){
	//	float p = this->pqueue->DeleteId( currentId );
	//	if ( p == VTK_LARGE_FLOAT ){
	//		cout << "point not in queue" << endl;
	//	}
	//}

	//// insert in queue
	//this->pqueue->Insert( this->priority(current), currentId );
}


template< class Pixel >
int nbfFastMarching2D< Pixel > :: array2int( TinyVector< int, 2 > & t ){
	return t(0) * this->distance.cols() + t(1);
}

template< class Pixel >
TinyVector< int, 2 > nbfFastMarching2D< Pixel > :: int2array( int p ){
	return TinyVector< int, 2 >( floor( 1.0 * p / this->distance.cols() ),
								 fmod( 1.0 * p, this->distance.cols() ) );
}


#endif /* FILE_nbfFastMarching2D */