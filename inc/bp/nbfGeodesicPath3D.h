#ifndef FILE_nbfGeodesicPath3D
#define FILE_nbfGeodesicPath3D

// Class nbfGeodesicPath3D.
//
// Find geodesic paths in 3D images.
// 

#include <nbfArray.h>

#include <vector>
#include <algorithm>

template< class Pixel >
class nbfGeodesicPath3D : public nbfArray< Pixel, 3 >
{
public:

	// takes $g$ array and center point as input
	nbfGeodesicPath3D( Array< Pixel, 3 > & );

	~nbfGeodesicPath3D(){};

	// get geodesic curve given start point
	// the end point is assumed to have zero distance
	// return true if succesful, false otherwise

	// (weights, start, path)
	bool getPath( Array< Pixel, 3 > &, TinyVector< int, 3 > &, vector< TinyVector< int, 3 > > & );
	
};


template< class Pixel >
nbfGeodesicPath3D< Pixel > :: nbfGeodesicPath3D( Array< Pixel, 3 > & weight )
: nbfArray< Pixel, 3 >( weight )
{
}

template< class Pixel >
bool nbfGeodesicPath3D< Pixel> :: getPath( Array< Pixel, 3 > & weights, 
									       TinyVector< int, 3 > & start,
										   vector< TinyVector< int, 3 > > & path )
{	
	// reset if not empty
	path.clear();

	TinyVector< int, 3 > next;
	
	vector< TinyVector< int, 3 > > neighbors;

	// back propagation - minimum distance neighbor
	next = start;
	path.push_back( next );
	
	Pixel minDistance;

	// stop when reach null distance
	while ( weights( next ) != 0 )
	{
		bool hasMultipleMinima = false;

		this->getFullNeighbors( next, neighbors );

		TinyVector< int, 3 > current = next;
		
		if ( neighbors.size() < 1 ){
			break;
		}
		else
		{

			// look for neighbor position with greatest gradient
			next = neighbors[0];
			Pixel minGradient = ( weights(next) - weights(current) ) / sqrt( pow2(next[0]-current[0])+pow2(next[1]-current[1])+pow2(next[2]-current[2]) + 0.0 );
			for ( int i = 1; i < neighbors.size(); i++ ){
				Pixel gradient = ( weights(neighbors[i]) - weights(current) ) / sqrt(pow2(neighbors[i][0]-current[0])+pow2(neighbors[i][1]-current[1])+pow2(neighbors[i][2]-current[2]) + 0.0);
				if ( gradient < minGradient ){
					minGradient = gradient;
					next = neighbors[i];
				}
			}


			//minDistance = weights(neighbors[0]);
			//next = neighbors[0];
			//for ( int i = 1; i < neighbors.size(); i++ ){
			//	if ( weights( neighbors[i] ) < minDistance ){
			//		minDistance = weights( neighbors[i] );
			//		next = neighbors[i];
			//		hasMultipleMinima = false;
			//	}
			//	else if ( weights( neighbors[i] ) == minDistance ){ 
			//		hasMultipleMinima = true;
			//		assert(0);
			//	}
			//}

			path.push_back( next );
			cout << next << " - " << weights(next) << endl;
		}
	}
	return true;
}


#endif /* FILE_nbfGeodesicPath3D */