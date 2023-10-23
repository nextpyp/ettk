#ifndef FILE_nbfNarrowBand
#define FILE_nbfNarrowBand

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vector>

class nbfNarrowBand
{

public:

	// build narrow band : take the full array, the list of positions in the NB and the th value.
	template< class Pixel, int const Dim >
	static void build( Array< Pixel, Dim > &, std::vector< TinyVector< int, 2 > > &, Pixel );

};

template< class Pixel, int const Dim >
void nbfNarrowBand :: build( Array< Pixel, Dim > & A, std::vector< TinyVector< int, 2 > > & I, Pixel th ){
	
	I.clear();
	Array< Pixel, Dim > :: iterator iter = A.begin();
	while ( iter != A.end() ){
		if ( abs( (*iter) ) < th ){
			I.push_back( iter.position() );
		}
		++iter;
	}
}

#endif // FILE_nbfNarrowBand
