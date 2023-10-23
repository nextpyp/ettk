#ifndef FILE_nbfFastMarching3D26
#define FILE_nbfFastMarching3D26

#include <fm/nbfFastMarching3D.h>
#include <nbfLinearInterpolator.h>

// Class nbfFastMarching3D26.
//
// Implements 3D fast marching method with arbitrary weights.
// 

template< class Pixel >
class nbfFastMarching3D26 : public nbfFastMarching3D< Pixel >
{
public:

	// constructor takes weight array as input
	nbfFastMarching3D26( Array< Pixel, 3 > & );
	nbfFastMarching3D26( Array< Pixel, 3 > &, TinyVector< int, 3 > & );
	~nbfFastMarching3D26(){};

protected:

	// update individual point
	// the first argument is the point to update
	// the second argument is the point from where we do the update
	void updatePoint( TinyVector< int, 3 > &, TinyVector< int, 3 > & );

	Pixel solve3D( TinyVector< int, 3 > &, TinyVector< int, 3 > &, TinyVector< int, 3 > &, TinyVector< int, 3 > & );

	virtual void get26Neighbors( TinyVector< int, 3 > &, TinyVector< int, 3 > &, vector< TinyVector< int, 3 > > &, vector< TinyVector< int, 3 > > & );

};

template< class Pixel >
nbfFastMarching3D26< Pixel > :: nbfFastMarching3D26( Array< Pixel, 3 > & weight )
: nbfFastMarching3D< Pixel >( weight )
{
}

template< class Pixel >
nbfFastMarching3D26< Pixel > :: nbfFastMarching3D26( Array< Pixel, 3 > & weight, TinyVector< int, 3 > & center )
: nbfFastMarching3D< Pixel >( weight, center )
{
}

template< class Pixel >
void nbfFastMarching3D26< Pixel > :: updatePoint( TinyVector< int, 3 > & current, TinyVector< int, 3 > & alive )
{
	vector< TinyVector< int, 3 > > neighborsB;
	vector< TinyVector< int, 3 > > neighborsC;

	this->get26Neighbors( current, alive, neighborsB, neighborsC );

	// for each 6-connected neighbor
	for ( int i = 0; i < neighborsB.size(); i++ ){
		Pixel uxm = this->solve3D( alive, neighborsB[i], neighborsC[i], current );
	
		// update distance value
		this->distance( current ) = min( this->distance( current ), uxm );
	}

	this->updateQueue( current );
}

using namespace blitz::extrema;

template< class Pixel >
Pixel nbfFastMarching3D26< Pixel > :: solve3D( TinyVector< int, 3 > & A, TinyVector< int, 3 > & B, TinyVector< int, 3 > & C, TinyVector< int, 3 > & current )
{
	// get distances
	Pixel uA, uB, uC;

	uA = nbfFastMarching<Pixel,3>::getDistance( A );	  // closest
	uB = nbfFastMarching<Pixel,3>::getDistance( B );	  // face diagonal
	uC = nbfFastMarching<Pixel,3>::getDistance( C );	  // body diagonal

	// use interpolated weights
#if 0
	// get interpolated weight
	nbfLinearInterpolator< Pixel, 3 > interpolator( this->weight );
	TinyVector< Pixel, 3 > center;
	center(firstDim) = ( (Pixel)C(firstDim) + (Pixel)current(firstDim) ) / 2.0;
	center(secondDim) = ( (Pixel)C(secondDim) + (Pixel)current(secondDim) ) / 2.0;
	center(thirdDim) = ( (Pixel)C(thirdDim) + (Pixel)current(thirdDim) ) / 2.0;

	Pixel t = interpolator.interpolateSingle(center);
#else
	Pixel t = this->weight( current );
#endif

	Pixel k1 = uA - uB;
	Pixel k2 = uB - uC;

	Pixel t1, t3;
	Pixel uD;
	Pixel delta;

	uD = min( uA + t, min( uB + sqrt(2.0) * t, uC + sqrt(3.0) * t ) );

	// ABC
	delta = t*t - (k1*k1+k2*k2);
	if ( delta > 0 ){
		t1 = 1 - k1 / sqrt( delta );
		t3 = k2 / sqrt( delta );

		if ( ( t1 > 0 ) & ( t3 > 0 ) & ( t1 + t3 < 1 ) )
		{
			uD = min( uD, uA + sqrt( delta ) );
		}
	}

	// AB
	delta = t*t - k1*k1;
	if ( delta > 0 ){
		t1 = 1 - k1 / sqrt( delta );
		if ( ( t1 > 0 ) & ( t1 < 1 ) ){
			uD = min( uD, uA + sqrt( delta ) );
		}
	}

	// BC
	delta = t*t -k2*k2;
	if ( delta > 0 ){
		t3 = sqrt(2.0) * k2 / sqrt( delta );
		if ( ( t3 > 0 ) & ( t3 < 1 ) )
		{
			uD = min( uD, uB + sqrt(2.0) * sqrt( delta ) );
		}
	}

	// AC
	delta = t*t - ( uA - uC ) * ( uA - uC ) / 2;
	if ( delta > 0 ){
		t3 = ( uA - uC ) / 2 / sqrt( delta );
		t1 = 1 - t3;
		if ( ( t1 > 0 ) & ( t3 > 0 ) )
		{
			uD = min( uD, uA + sqrt( delta ) );
		}
	}

	return uD;
}

template< class Pixel >
void nbfFastMarching3D26< Pixel > :: get26Neighbors( TinyVector< int, 3 > & current, 
													 TinyVector< int, 3 > & alive, 
													 vector< TinyVector< int, 3 > > & neighborsB, 
													 vector< TinyVector< int, 3 > > & neighborsC )
{
	neighborsB.clear();
	neighborsC.clear();

	if ( current(firstDim) != alive(firstDim) ){
		neighborsB.push_back( alive + TinyVector<int,3>(0,1,0) );
		neighborsC.push_back( alive + TinyVector<int,3>(0,1,1) );

		neighborsB.push_back( alive + TinyVector<int,3>(0,1,0) );
		neighborsC.push_back( alive + TinyVector<int,3>(0,1,-1) );

		neighborsB.push_back( alive + TinyVector<int,3>(0,-1,0) );
		neighborsC.push_back( alive + TinyVector<int,3>(0,-1,1) );

		neighborsB.push_back( alive + TinyVector<int,3>(0,-1,0) );
		neighborsC.push_back( alive + TinyVector<int,3>(0,-1,-1) );

		neighborsB.push_back( alive + TinyVector<int,3>(0,0,1) );
		neighborsC.push_back( alive + TinyVector<int,3>(0,1,1) );

		neighborsB.push_back( alive + TinyVector<int,3>(0,0,1) );
		neighborsC.push_back( alive + TinyVector<int,3>(0,-1,1) );

		neighborsB.push_back( alive + TinyVector<int,3>(0,0,-1) );
		neighborsC.push_back( alive + TinyVector<int,3>(0,1,-1) );

		neighborsB.push_back( alive + TinyVector<int,3>(0,0,-1) );
		neighborsC.push_back( alive + TinyVector<int,3>(0,-1,-1) );

	}
	else{
		if ( current(secondDim) != alive(secondDim) ){
			neighborsB.push_back( alive + TinyVector<int,3>(1,0,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(1,0,1) );

			neighborsB.push_back( alive + TinyVector<int,3>(1,0,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(1,0,-1) );

			neighborsB.push_back( alive + TinyVector<int,3>(-1,0,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(-1,0,1) );

			neighborsB.push_back( alive + TinyVector<int,3>(-1,0,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(-1,0,-1) );

			neighborsB.push_back( alive + TinyVector<int,3>(0,0,1) );
			neighborsC.push_back( alive + TinyVector<int,3>(1,0,1) );

			neighborsB.push_back( alive + TinyVector<int,3>(0,0,1) );
			neighborsC.push_back( alive + TinyVector<int,3>(-1,0,1) );

			neighborsB.push_back( alive + TinyVector<int,3>(0,0,-1) );
			neighborsC.push_back( alive + TinyVector<int,3>(1,0,-1) );

			neighborsB.push_back( alive + TinyVector<int,3>(0,0,-1) );
			neighborsC.push_back( alive + TinyVector<int,3>(-1,0,-1) );
		}
		else{
			neighborsB.push_back( alive + TinyVector<int,3>(0,1,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(1,1,0) );

			neighborsB.push_back( alive + TinyVector<int,3>(0,1,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(-1,1,0) );

			neighborsB.push_back( alive + TinyVector<int,3>(0,-1,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(1,-1,0) );

			neighborsB.push_back( alive + TinyVector<int,3>(0,-1,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(-1,-1,0) );

			neighborsB.push_back( alive + TinyVector<int,3>(1,0,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(1,1,0) );

			neighborsB.push_back( alive + TinyVector<int,3>(1,0,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(1,-1,0) );

			neighborsB.push_back( alive + TinyVector<int,3>(-1,0,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(-1,1,0) );

			neighborsB.push_back( alive + TinyVector<int,3>(-1,0,0) );
			neighborsC.push_back( alive + TinyVector<int,3>(-1,-1,0) );
		}
	}

	if ( this->hasCut )
	{
		// now get rid of the ones that are across the cut
		vector< TinyVector< int, 3 > > :: iterator iterB = neighborsB.begin();
		vector< TinyVector< int, 3 > > :: iterator iterC = neighborsC.begin();

		while ( iterB != neighborsB.end() ){
			// check if not on the cut
			if ( this->isAtCut( current, *iterB ) || this->isAtCut( current, *iterC ) ){
				iterB = neighborsB.erase( iterB );
				iterC = neighborsC.erase( iterC );
			}
			else{
				++iterB;
				++iterC;
			}
		}
	}
}


#endif /* FILE_nbfFastMarching3D26 */