#ifndef FILE_nbfFastGeodesicPath
#define FILE_nbfFastGeodesicPath

// Class nbfFastGeodesicPath.
//
// Find geodesic paths in 2D images.
// 

#include <fm/nbfFastMarchingFool.h>
#include <nbfLinearInterpolator.h>

#include <vector>
#include <algorithm>

template< class Pixel >
class nbfFastGeodesicPath
{
public:

	// takes $g$ array and center point as input
	nbfFastGeodesicPath( Array< Pixel, 2 > & );

	~nbfFastGeodesicPath(){};

	// get geodesic curve given start point
	// the end point is assumed to have zero distance
	// return true if succesful, false otherwise

	bool getPath( TinyVector< int, 2 > &, vector< TinyVector< Pixel, 2 > > & );

	TinyVector< Pixel, 2 > getPointGradient( TinyVector< int, 2 > & );

	// get path as implicit function (only closed paths)
	void getImplicitPath( vector< TinyVector< int, 2 > > &, Array< Pixel, 2 > & );
	void getImplicitPath( vector< TinyVector< Pixel, 2 > > &, Array< Pixel, 2 > &, Pixel = numeric_limits<Pixel>::max() );

protected:

	// weight array
	Array< Pixel, 2 > weight;
};


template< class Pixel >
nbfFastGeodesicPath< Pixel > :: nbfFastGeodesicPath( Array< Pixel, 2 > & weight )
{
	this->weight.reference( weight );
}


template< class Pixel >
TinyVector< Pixel, 2 > nbfFastGeodesicPath< Pixel> :: getPointGradient( TinyVector< int, 2 > & position )
{	
	int x = position[firstDim];
	int y = position[secondDim];

	TinyVector< Pixel, 2 > gradient = 0;

	Pixel d;
	if ( this->weight.isInRange( position ) ){
		d = this->weight( position );
	}
	else{
		return gradient;
	}

	if ( this->weight.isInRange(x-1,y) ){
		gradient[firstDim] = min( gradient[firstDim], this->weight(x-1,y) - d );
	}
	if ( this->weight.isInRange(x+1,y) ){
		gradient[firstDim] = min( gradient[firstDim], this->weight(x+1,y) - d );
	}
	if ( this->weight.isInRange(x,y-1) ){
		gradient[secondDim] = min( gradient[secondDim], this->weight(x,y-1) - d );
	}
	if ( this->weight.isInRange(x,y+1) ){
		gradient[secondDim] = min( gradient[secondDim], this->weight(x,y+1) - d );
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
bool nbfFastGeodesicPath< Pixel> :: getPath( TinyVector< int, 2 > & start,
								             vector< TinyVector< Pixel, 2 > > & path )
{	
	// reset path if not empty
	path.clear();

	// store next point in geodesic
	TinyVector< Pixel, 2 > next;
	
	// push first point into path
	next = start;
	path.push_back( next );
	
	// posta
	Array< Pixel, 2 > dx( this->weight.shape() );
	Array< Pixel, 2 > dy( this->weight.shape() );

	dx = numeric_limits<Pixel>::max();
	dy = numeric_limits<Pixel>::max();
	
	nbfLinearInterpolator< Pixel, 2 > interpX( dx );
	nbfLinearInterpolator< Pixel, 2 > interpY( dy );

	vector< TinyVector< int, 2 > > neighbors;

	while ( this->weight( next ) > 0 ){

		// made last position coordinates handy
		Pixel x = next[firstDim];
		Pixel y = next[secondDim];

		// get square box where the point is contained
		int lx = floor( x ); int ux = ceil( x );
		int ly = floor( y ); int uy = ceil( y );

		// Before interpolating the gradient, make sure all adjacent points 
		// have a gradient value (costly operation)
		neighbors.clear();
		neighbors.push_back( TinyVector< int, 2 >( lx, ly ) );
		neighbors.push_back( TinyVector< int, 2 >( lx, uy ) );
		neighbors.push_back( TinyVector< int, 2 >( ux, ly ) );
		neighbors.push_back( TinyVector< int, 2 >( ux, uy ) );

		vector< TinyVector< int, 2 > > :: iterator iter = neighbors.begin();
		while ( iter != neighbors.end() ){
			if ( dx.isInRange(*iter) == true ){
				if ( dx( *iter ) == numeric_limits<Pixel>::max() ){
					TinyVector< Pixel, 2 > grad;
					grad = this->getPointGradient( *iter );
					dx( *iter ) = grad[firstDim];
					dy( *iter ) = grad[secondDim];
				}
			}
			++iter;
		}

		TinyVector< Pixel, 2 > gradient( interpX.interpolateSingle( next ), 
			                             interpY.interpolateSingle( next ) );

		// re-normalize gradient vector
		Pixel norm = sqrt( pow2( gradient[firstDim] ) + pow2( gradient[secondDim] ) );

		if ( norm > 0 ){
			gradient = gradient / norm;
		}
		else{
			cout << next << endl;
			cout << gradient << endl;
			for ( int i = 0; i < neighbors.size(); i++ ){
				cout << neighbors[i] << ", [" << dx(neighbors[i]) << "," << dy(neighbors[i]) << "]\n";
			}
		}

		// es muy poco sensible al tamaño del paso, achicarlo solo genera discrepancias minimas
		Pixel lambda = .25;

		next = next + lambda * gradient;
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

		int nfx = floor(next[firstDim]);
		int ncx = ceil(next[firstDim]);
		int nfy = floor(next[secondDim]);
		int ncy = ceil(next[secondDim]);

		if ( this->weight( nfx, nfy ) == 0 ){
			path.push_back( next );
			next = TinyVector< int, 2 >(nfx,nfy);
		}
		else{
			if ( this->weight( nfx, ncy ) == 0 ){
				path.push_back( next );
				next = TinyVector< int, 2 >(nfx,ncy);
			}
			else{
				if ( this->weight( ncx, nfy ) == 0 ){
					path.push_back( next );
					next = TinyVector< int, 2 >(ncx,nfy);
				}
				else{
					if ( this->weight( ncx, ncy ) == 0 ){
						path.push_back( next );
						next = TinyVector< int, 2 >(ncx,ncy);
					}
				}
			}
		}

		path.push_back( next );
		if ( path.size() > 15000 ){
			cout << next << endl;
			return false;
		}
	}
}


template< class Pixel >
void nbfFastGeodesicPath< Pixel > :: getImplicitPath( vector< TinyVector< int, 2 > > & path,
												      Array< Pixel, 2 > & A )
{
	if ( path.size() == 0 ){
		A = numeric_limits<Pixel>::max();
	}
	else{
		// if open path then switch to simple point representation
		if ( ( path[0](firstDim) - path[ path.size() - 1 ](firstDim) > 1 ) |
             ( path[0](secondDim) - path[ path.size() - 1 ](secondDim) > 1 ) )
		{
			A = -1;
			for ( int i = 0; i < path.size(); i++ )
			{
				A(path[i]) = 1;
			}
		}
		else
		{
			Array< Pixel, 2 > D( A.shape() );

			// compute iso-distances to path first
			A = 1;
			nbfFastFastMarching2D< Pixel > fm2d(A);

			vector< TinyVector< int, 2 > > positions;
			vector< Pixel > distances;
			for ( int i = 0; i < path.size(); i++ ){
				positions.push_back( path[i] );
				distances.push_back( 0 );
			}
			fm2d.setAliveSet(positions,distances);
			//fm2d.setStopDistance(10);
			fm2d.execute(D);

			// now compute inside/outside
			A[path] = numeric_limits<Pixel>::max();
			TinyVector< int, 2 > alive( this->center(firstDim), this->center(secondDim) );
			nbfFastMarchingFool2D< Pixel > fm(A);
			fm.setAliveSet( alive, 0 );
			Array< Pixel, 2 > S( A.shape() );
			fm.execute(S);
			// if not sufficient points inside
			if ( sum( S == numeric_limits<Pixel>::max() ) < 5 )
			{
				A = 1;
			}
			else
			{
				A = where( S < numeric_limits<Pixel>::max(), -D, D );
			}
		}
	}
}


template< class Pixel >
void nbfFastGeodesicPath< Pixel > :: getImplicitPath( vector< TinyVector< Pixel, 2 > > & path,
												  Array< Pixel, 2 > & A,
												  Pixel maxDistance )
{
	if ( path.size() == 0 ){
		A = numeric_limits<Pixel>::max();
	}
	else{
		Array< Pixel, 2 > D( A.shape() );
		D = numeric_limits<Pixel>::max();

		Array< Pixel, 2 > S( A.shape() );
		S = 1;

		TinyVector< Pixel, 2 > position;
		
		TinyVector< Pixel, 2 > gradient;

		// build inside/outside regions
		for ( int i = 0; i < path.size(); i++ ){

			Pixel cX = path[i](firstDim);
			Pixel cY = path[i](secondDim);

			int lx = floor( cX ); int ux = ceil( cX );
			int ly = floor( cY ); int uy = ceil( cY );
	
			if ( i < path.size() - 1 ){
				gradient = path[i+1] - path[i];
			}

			if ( D.isInRange(lx,ly) ){
				D( lx, ly ) = min( D( lx, ly ), sqrt( pow2(cX-lx) + pow2(cY-ly) ) );
				// handle this appart
				if ( ( lx == ux ) | ( ly == uy ) ){
					S(lx,ly) = numeric_limits<Pixel>::max();
				}
				else{
					position = TinyVector< Pixel, 2 >(lx,ly) - path[i];
					if ( gradient[1] * position[0] - gradient[0] * position[1] < 0 ){
						S( lx, ly ) = numeric_limits<Pixel>::max();
					}
				}
			}
			if ( D.isInRange(lx,uy) ){
				D( lx, uy ) = min( D( lx, uy ), sqrt( pow2(cX-lx) + pow2(cY-uy) ) );
				position = TinyVector< Pixel, 2 >(lx,uy) - path[i];
				if ( gradient[1] * position[0] - gradient[0] * position[1] < 0){
					S( lx, uy ) = numeric_limits<Pixel>::max();
				}
			}
			if ( D.isInRange(ux,ly) ){
				D( ux, ly ) = min( D( ux, ly ), sqrt( pow2(cX-ux) + pow2(cY-ly) ) );
				position = TinyVector< Pixel, 2 >(ux,ly) - path[i];
				if ( gradient[1] * position[0] - gradient[0] * position[1] < 0 ){
					S( ux, ly ) = numeric_limits<Pixel>::max();
				}
			}
			if ( D.isInRange(ux,uy) ){
				D( ux, uy ) = min( D( ux, uy ), sqrt( pow2(cX-ux) + pow2(cY-uy) ) );
				position = TinyVector< Pixel, 2 >(ux,uy) - path[i];
				if ( gradient[1] * position[0] - gradient[0] * position[1] < 0 ){
					S( ux, uy ) = numeric_limits<Pixel>::max();
				}
			}			
		}

		if ( maxDistance > 1 ){
			// compute iso-distances to path first
			A = 1;
			nbfFastMarching2D< Pixel > fm2d(A);

			vector< TinyVector< int, 2 > > positions;
			vector< Pixel > distances;

			Array< Pixel, 2 > :: iterator iter = D.begin();
			while( iter != D.end() ){
				if ( (*iter) < numeric_limits<Pixel>::max() ){
					positions.push_back( iter.position() );
					distances.push_back( *iter );
				}
				++iter;
			}
			fm2d.setAliveSet(positions,distances);
			fm2d.setStopDistance(maxDistance);
			fm2d.execute(D);
		}

		// now compute inside/outside

		TinyVector< int, 2 > alive( 0, 0 );
		nbfFastMarchingFool2D< Pixel > fm(S);
		fm.setAliveSet( alive, 0 );
		Array< Pixel, 2 > sign( A.shape() );
		fm.execute(sign);

		A = where( sign < numeric_limits<Pixel>::max(), D, -D );
	}	
}


#endif /* FILE_nbfFastGeodesicPath */