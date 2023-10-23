#ifndef FILE_nbfFastMarchingFool3D
#define FILE_nbfFastMarchingFool3D

// Class nbfFastMarchingFool3D.
//
// 3D region growing
// 

#include <fm/nbfFastMarching3D26.h>

template< class Pixel >
class nbfFastMarchingFool3D : public nbfFastMarching3D< Pixel >
{
public:

	// constructor takes weight array as input
	nbfFastMarchingFool3D( Array< Pixel, 3 > &, Pixel );
	~nbfFastMarchingFool3D(){};

protected:

	// the second argument is the point from where we do the update
	void updatePoint( TinyVector< int, 3 > &, TinyVector< int, 3 > & );

	Pixel threshold;
};

template< class Pixel >
nbfFastMarchingFool3D< Pixel > :: nbfFastMarchingFool3D( Array< Pixel, 3 > & weight, Pixel th )
: nbfFastMarching3D< Pixel >( weight )
{ this->threshold = th; }

template< class Pixel >
void nbfFastMarchingFool3D< Pixel > :: updatePoint( TinyVector< int, 3 > & current, TinyVector< int, 3 > & alive )
{
	// if inside, assign priority according to gray values
	if ( this->weight( current ) < this->threshold ){
		this->distance( current ) = this->weight( current );
		this->updateQueue( current );
	}
	
	//// if in FAR_SET move to TRIAL_SET and push to queue
	//if ( ( this->state( current ) == FAR_SET ) &&  ( (*this->weight)( current ) < numeric_limist<Pixel>::max() ) ){
	//	this->state( current ) = TRIAL_SET;
	//	nbfDataPoint< Pixel, 2 > tmp( current );
	//	tmp.distance = &((*this->distance)(current));
	//	pqueue.push( tmp );
	//}
}


#endif /* FILE_nbfFastMarchingFool2D */