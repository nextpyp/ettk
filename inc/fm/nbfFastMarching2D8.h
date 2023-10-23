#ifndef FILE_nbfFastMarching2D8
#define FILE_nbfFastMarching2D8

#include <fm/nbfFastMarching2D.h>

// Class nbfFastMarching2D8.
//
// Implements 2D fast marching method with arbitrary weights.
// 

template< class Pixel >
class nbfFastMarching2D8 : public nbfFastMarching2D< Pixel >
{
public:

	// constructor takes weight array as input
	nbfFastMarching2D8( Array< Pixel, 2 > & );
	nbfFastMarching2D8( Array< Pixel, 2 > &, TinyVector< int, 2 > & );
	~nbfFastMarching2D8(){};

protected:

	// update individual point
	// the first argument is the point to update
	// the second argument is the point from where we do the update
	void updatePoint( TinyVector< int, 2 > &, TinyVector< int, 2 > & );

	virtual Pixel solve2D( TinyVector< int, 2 > &, TinyVector< int, 2 > &, TinyVector< int, 2 > & );

};

template< class Pixel >
nbfFastMarching2D8< Pixel > :: nbfFastMarching2D8( Array< Pixel, 2 > & weight )
: nbfFastMarching2D< Pixel >( weight )
{
}

template< class Pixel >
nbfFastMarching2D8< Pixel > :: nbfFastMarching2D8( Array< Pixel, 2 > & weight, TinyVector< int, 2 > & center )
: nbfFastMarching2D< Pixel >( weight, center )
{
}

template< class Pixel >
void nbfFastMarching2D8< Pixel > :: updatePoint( TinyVector< int, 2 > & current, TinyVector< int, 2 > & alive )
{
	// coordinates of point we are updating
	int x = alive(firstDim);
	int y = alive(secondDim);

	int dimension;

	if ( current(firstDim) == alive(firstDim) ){
		dimension = firstDim;
	}
	else{
		dimension = secondDim;
	}

	TinyVector< int, 2 > first, second;

	if ( dimension == firstDim ){
		first(firstDim) = x + 1; first(secondDim) = y;
		second(firstDim) = x - 1; second(secondDim) = y;
	}
	else{
		first(firstDim) = x; first(secondDim) = y + 1;
		second(firstDim) = x; second(secondDim) = y - 1;
	}

	Pixel dFirst = this->solve2D(current,alive,first);
	Pixel dSecond = this->solve2D(current,alive,second);

	// update distance value
	this->distance( current ) = min( this->distance( current ), min( dFirst, dSecond ) );

	this->updateQueue( current );

}

template< class Pixel >
Pixel nbfFastMarching2D8< Pixel > :: solve2D( TinyVector< int, 2 > & current, TinyVector< int, 2 > & alive, TinyVector< int, 2 > & extra )
{
	// get distances
	Pixel uA, uB;
	if ( this->hasCut & this->isAtCut( current, extra ) ){
		uA = numeric_limits< Pixel > :: max();
	}
	else
	{
		uA = nbfFastMarching<Pixel,2>::getDistance( extra );    // diagonal
	}

	if ( this->hasCut & this->isAtCut( current, alive ) ){
		uB = numeric_limits< Pixel > :: max();
	}
	else
	{
		uB = nbfFastMarching<Pixel,2>::getDistance( alive );    // 4-connected
	}

	// retrieve weight
	Pixel tABC = this->weight( current );

	// tentative distance value
	Pixel uAB;

	if ( uA < uB ){
		if ( tABC <= sqrt(2.0) * ( uB - uA ) ){
			return uA + sqrt(2.0) * tABC;
		}
		else{
			return uB + sqrt( tABC * tABC - pow2( uB - uA ) );
		}
	}
	else{
		return uB + tABC;
	}
}

#endif /* FILE_nbfFastMarching2D8 */