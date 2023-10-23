#ifndef FILE_nbfFastGeodesicPath3D
#define FILE_nbfFastGeodesicPath3D

// Class nbfFastGeodesicPath.
//
// Find geodesic paths in 3D images.
// 

#include <fm/nbfFastMarchingFool.h>
#include <nbfLinearInterpolator.h>

#include <vector>
#include <algorithm>

template< class Pixel >
class nbfFastGeodesicPath3D
{
public:

	// takes $g$ array and center point as input
	nbfFastGeodesicPath3D( Array< Pixel, 3 > & );

	~nbfFastGeodesicPath3D(){};

	// get geodesic curve given start point
	// the end point is assumed to have zero distance
	// return true if succesful, false otherwise

	bool getPath( TinyVector< int, 3 > &, vector< TinyVector< Pixel, 3 > > & );

	TinyVector< Pixel, 3 > getPointGradient( TinyVector< int, 3 > & );

	// get path as implicit function (only closed paths)
	//void getImplicitPath( vector< TinyVector< int, 3 > > &, Array< Pixel, 3 > & );
	//void getImplicitPath( vector< TinyVector< Pixel, 3 > > &, Array< Pixel, 3 > &, Pixel = numeric_limits<Pixel>::max() );

protected:

	// weight array
	Array< Pixel, 3 > weight;
};


template< class Pixel >
nbfFastGeodesicPath3D< Pixel > :: nbfFastGeodesicPath3D( Array< Pixel, 3 > & weight )
{
	this->weight.reference( weight );
}


template< class Pixel >
TinyVector< Pixel, 3 > nbfFastGeodesicPath3D< Pixel> :: getPointGradient( TinyVector< int, 3 > & position )
{	
	int x = position[firstDim];
	int y = position[secondDim];
	int z = position[thirdDim];

	TinyVector< Pixel, 3 > gradient = 0;

	Pixel d;
	if ( this->weight.isInRange( position ) ){
		d = this->weight( position );
	}
	else{
		return gradient;
	}

	if ( this->weight.isInRange(x-1,y,z) ){
		gradient[firstDim] = min( gradient[firstDim], this->weight(x-1,y,z) - d );
	}
	if ( this->weight.isInRange(x+1,y,z) ){
		gradient[firstDim] = min( gradient[firstDim], this->weight(x+1,y,z) - d );
	}
	if ( this->weight.isInRange(x,y-1,z) ){
		gradient[secondDim] = min( gradient[secondDim], this->weight(x,y-1,z) - d );
	}
	if ( this->weight.isInRange(x,y+1,z) ){
		gradient[secondDim] = min( gradient[secondDim], this->weight(x,y+1,z) - d );
	}
	if ( this->weight.isInRange(x,y,z-1) ){
		gradient[thirdDim] = min( gradient[thirdDim], this->weight(x,y,z-1) - d );
	}
	if ( this->weight.isInRange(x,y,z+1) ){
		gradient[thirdDim] = min( gradient[thirdDim], this->weight(x,y,z+1) - d );
	}
	return gradient;
}


// This back propagation (state-of-the-art) first computes a gradient image
// and then interpolates to get the gradient at intermediate steps.
// The gradient is computed assuming that distances were computed with 8-neighbor fast
// marching, so the gradient direction at each pixel is that from where the distance
// was updated.
// We now bilinearly interpolate the x and y components of the gradient image
// to get the advancing direction at arbitrary point in the image domain.

template< class Pixel >
bool nbfFastGeodesicPath3D< Pixel> :: getPath( TinyVector< int, 3 > & start,
											  vector< TinyVector< Pixel, 3 > > & path )
{	
	// reset path if not empty
	path.clear();

	// store next point in geodesic
	TinyVector< Pixel, 3 > next;
	
	// push first point into path
	next = start;
	path.push_back( next );
	
	// posta
	Array< Pixel, 3 > dx( this->weight.shape() );
	Array< Pixel, 3 > dy( this->weight.shape() );
	Array< Pixel, 3 > dz( this->weight.shape() );

	dx = numeric_limits<Pixel>::max();
	dy = numeric_limits<Pixel>::max();
	dz = numeric_limits<Pixel>::max();
	
	nbfLinearInterpolator< Pixel, 3 > interpX( dx );
	nbfLinearInterpolator< Pixel, 3 > interpY( dy );
	nbfLinearInterpolator< Pixel, 3 > interpZ( dz );

	nbfLinearInterpolator< Pixel, 3 > interpD( this->weight );

	vector< TinyVector< int, 3 > > neighbors;

	while ( this->weight( next ) > 0 ){

		//cout << this->weight( next ) << endl;

		// made last position coordinates handy
		Pixel x = next[firstDim];
		Pixel y = next[secondDim];
		Pixel z = next[thirdDim];

		// get square box where the point is contained
		int lx = floor( x ); int ux = ceil( x );
		int ly = floor( y ); int uy = ceil( y );
		int lz = floor( z ); int uz = ceil( z );

		// Before interpolating the gradient, make sure all adjacent points 
		// have a gradient value (costly operation)
		neighbors.clear();
		neighbors.push_back( TinyVector< int, 3 >( lx, ly, lz ) );
		neighbors.push_back( TinyVector< int, 3 >( lx, uy, lz ) );
		neighbors.push_back( TinyVector< int, 3 >( ux, ly, lz ) );
		neighbors.push_back( TinyVector< int, 3 >( ux, uy, lz ) );
		neighbors.push_back( TinyVector< int, 3 >( lx, ly, uz ) );
		neighbors.push_back( TinyVector< int, 3 >( lx, uy, uz ) );
		neighbors.push_back( TinyVector< int, 3 >( ux, ly, uz ) );
		neighbors.push_back( TinyVector< int, 3 >( ux, uy, uz ) );

		vector< TinyVector< int, 3 > > :: iterator iter = neighbors.begin();
		while ( iter != neighbors.end() ){
			if ( dx.isInRange(*iter) == true ){
				if ( dx( *iter ) == numeric_limits<Pixel>::max() ){
					TinyVector< Pixel, 3 > grad;
					grad = this->getPointGradient( *iter );
					dx( *iter ) = grad[firstDim];
					dy( *iter ) = grad[secondDim];
					dz( *iter ) = grad[thirdDim];
				}
			}
			++iter;
		}

		TinyVector< Pixel, 3 > gradient( interpX.interpolateSingle( next ), 
			                             interpY.interpolateSingle( next ), 
			                             interpZ.interpolateSingle( next ) );

		// re-normalize gradient vector
		Pixel norm = sqrt( pow2( gradient[firstDim] ) + pow2( gradient[secondDim] ) + pow2( gradient[thirdDim] ) );

		if ( norm > 0 ){
			gradient = gradient / norm;
		}
		else{
			cout << "ERROR - " << next << endl;
			cout << gradient << endl;
			for ( int i = 0; i < neighbors.size(); i++ ){
				cout << neighbors[i] << ", [" << dx(neighbors[i]) << "," << dy(neighbors[i]) << "," << dz(neighbors[i]) << "]\n";
			}
			return false;
		}

		// es muy poco sensible al tamaño del paso, achicarlo solo genera discrepancias minimas
		Pixel lambda = .25;

		next = next - lambda * gradient;
		//next[firstDim] = next[firstDim] + lambda * gradient[firstDim];
		//next[secondDim] = next[secondDim] + lambda * gradient[secondDim];

		// keep inside domain
		if ( next[firstDim] < this->weight.lbound(firstDim) ){
			next[firstDim] = this->weight.lbound(firstDim);
		}
		else{
			if ( next[firstDim] > this->weight.ubound(firstDim) ){
				next[firstDim] = this->weight.ubound(firstDim);
			}
		}

		if ( next[secondDim] < this->weight.lbound(secondDim) ){
			next[secondDim] = this->weight.lbound(secondDim);
		}
		else{
			if ( next[secondDim] > this->weight.ubound(secondDim) ){
				next[secondDim] = this->weight.ubound(secondDim);
			}
		}

		if ( next[thirdDim] < this->weight.lbound(thirdDim) ){
			next[thirdDim] = this->weight.lbound(thirdDim);
		}
		else{
			if ( next[thirdDim] > this->weight.ubound(thirdDim) ){
				next[thirdDim] = this->weight.ubound(thirdDim);
			}
		}

		int nfx = floor(next[firstDim]);
		int ncx = ceil(next[firstDim]);
		int nfy = floor(next[secondDim]);
		int ncy = ceil(next[secondDim]);
		int nfz = floor(next[thirdDim]);
		int ncz = ceil(next[thirdDim]);

		if ( this->weight( nfx, nfy, nfz ) == 0 ){
			path.push_back( next );
			next = TinyVector< int, 3 >(nfx,nfy,nfz);
		}
		else{
			if ( this->weight( nfx, ncy, nfz ) == 0 ){
				path.push_back( next );
				next = TinyVector< int, 3 >(nfx,ncy,nfz);
			}
			else{
				if ( this->weight( ncx, nfy, nfz ) == 0 ){
					path.push_back( next );
					next = TinyVector< int, 3 >(ncx,nfy,nfz);
				}
				else{
					if ( this->weight( ncx, ncy, nfz ) == 0 ){
						path.push_back( next );
						next = TinyVector< int, 3 >(ncx,ncy,nfz);
					}
					else{
						if ( this->weight( nfx, nfy, ncz ) == 0 ){
							path.push_back( next );
							next = TinyVector< int, 3 >(nfx,nfy,ncz);
						}
						else{
							if ( this->weight( nfx, ncy, ncz ) == 0 ){
								path.push_back( next );
								next = TinyVector< int, 3 >(nfx,ncy,ncz);
							}
							else{
								if ( this->weight( ncx, nfy, ncz ) == 0 ){
									path.push_back( next );
									next = TinyVector< int, 3 >(ncx,nfy,ncz);
								}
								else{
									if ( this->weight( ncx, ncy, ncz ) == 0 ){
										path.push_back( next );
										next = TinyVector< int, 3 >(ncx,ncy,ncz);
									}
								}
							}
						}
					}
				}
			}
		}

		path.push_back( next );
		cout << next << ", d = " << interpD.interpolateSingle( next ) << endl;
		if ( path.size() > 500 ){
			cout << next << endl;
			return false;
		}
	}
}


#endif /* FILE_nbfFastGeodesicPath */