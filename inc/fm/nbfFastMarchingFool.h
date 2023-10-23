#ifndef FILE_nbfFastMarchingFool2D
#define FILE_nbfFastMarchingFool2D

// Class nbfFastMarchingFool.
//
// Implements 2D and 3D fast marching method with arbitrary weights.
// 

#include <fm/nbfFastMarching2D.h>

template< class Pixel >
class nbfFastMarchingFool2D : public nbfFastMarching2D< Pixel >
{
public:

	// constructor takes weight array as input
	nbfFastMarchingFool2D( Array< Pixel, 2 > & );
	~nbfFastMarchingFool2D(){};

protected:

	// the second argument is the point from where we do the update
	void updatePoint( TinyVector< int, 2 > &, TinyVector< int, 2 > & );

};

template< class Pixel >
nbfFastMarchingFool2D< Pixel > :: nbfFastMarchingFool2D( Array< Pixel, 2 > & weight )
: nbfFastMarching2D< Pixel >( weight )
{}

template< class Pixel >
void nbfFastMarchingFool2D< Pixel > :: updatePoint( TinyVector< int, 2 > & current, TinyVector< int, 2 > & alive )
{
	Pixel newD = this->distance( alive ) + 1;

	// update distance value
	this->distance( current ) = min( this->distance( current ), newD );

	this->updateQueue( current );

	//// if in FAR_SET move to TRIAL_SET and push to queue
	//if ( ( this->state( current ) == FAR_SET ) &&  ( (*this->weight)( current ) < numeric_limist<Pixel>::max() ) ){
	//	this->state( current ) = TRIAL_SET;
	//	nbfDataPoint< Pixel, 2 > tmp( current );
	//	tmp.distance = &((*this->distance)(current));
	//	pqueue.push( tmp );
	//}
}


#endif /* FILE_nbfFastMarchingFool2D */