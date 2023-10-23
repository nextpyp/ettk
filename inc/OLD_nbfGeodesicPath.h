#ifndef FILE_nbfGeodesicPath
#define FILE_nbfGeodesicPath

// Class nbfGeodesicPath.
//
// Find geodesic paths in 2D images.
// 

#include <nbfArray.h>

template< class Pixel >
class nbfGeodesicPath : nbfArray< Pixel, 2 >
{
private:

	Array< Pixel, 2 > * weight;
	int orientation;

public:

	// constructor takes weight array as input
	nbfGeodesicPath( Array< Pixel, 2 > &, int = firstDim );

	~nbfGeodesicPath(){};

	void setOrientation( int o ){ this->orientation = o; }

	// get geodesic curve given end point 
	// (pre-computed distances must be provided)
	// the end point is assumed to have zero distance
	void getPath( Array< Pixel, 2 > &, TinyVector< int, 2 >, TinyVector< int, 2 >, vector< TinyVector< int, 2 > > &, bool = false );
	void getPath( Array< Pixel, 2 > &, TinyVector< int, 2 >, TinyVector< int, 2 >, Array< Pixel, 2 > & );
	void getPath( Array< Pixel, 2 > &, TinyVector< int, 2 >, TinyVector< int, 2 >, vector< TinyVector< int, 2 > > &, Array< Pixel, 2 > & );

	// get minimal circular path
	void getCircularPath( vector< TinyVector< int, 2 > > & );
	void getCircularPath( Array< Pixel, 2 > & );
	void getCircularPath( vector< TinyVector< int, 2 > > &, Array< Pixel, 2 > & );

	void getForwardPath( Array< Pixel, 2 > &, TinyVector< int, 2 >, vector< TinyVector< int, 2 > > & );
	void getForwardCircularPath( vector< TinyVector< int, 2 > > &, Array< Pixel, 2 > & );

	// get path as implicit function
	void getImplicitPath( vector< TinyVector< int, 2 > > &, Array< Pixel, 2 > & );

};

template< class Pixel >
nbfGeodesicPath< Pixel > :: nbfGeodesicPath( Array< Pixel, 2 > & weight,
												  int orientation )
												  : nbfCheckedArray< Pixel, 2 >( weight )
{
	this->weight = &weight;
	this->orientation = orientation;
}


template< class Pixel >
void nbfGeodesicPath< Pixel > :: getCircularPath( vector< TinyVector< int, 2 > > & path )
{
	Array< Pixel, 2 > ipath( this->weight->shape() );
	this->getCircularPath( path, ipath );
}

template< class Pixel >
void nbfGeodesicPath< Pixel > :: getCircularPath( Array< Pixel, 2 > & ipath )
{
	vector< TinyVector< int, 2 > > path;
	this->getCircularPath( path, ipath );
}

template< class Pixel >
void nbfGeodesicPath< Pixel > :: getCircularPath( vector< TinyVector< int, 2 > > & path, Array< Pixel, 2 > & ipath )
{

	// find a circular path between parallel sides of a 2d image.
	// The circular constraint is enforced by restricting the start and
	// end points of the path to be circular neighbors, i.e. to have the
	// same coordinate in the appropriate dimension.
	nbfFastMarching< Pixel, 2 > fastMarching( *this->weight );

	// get infinity value from input image
//	float infty = max( *this->weight );
//	fastMarching.setRestrictDomain( infty );

	vector< TinyVector< int, 2 > > positions;
	vector< Pixel > distances;

	// build alive sets depending on desired orientation
	switch ( this->orientation ){
			case firstDim:
				for ( int j = 0; j < this->weight->cols(); j++ ){
					if ( (*this->weight)( this->weight->lbound(firstDim), j ) < numeric_limits<Pixel>::max() ){
						TinyVector< int, 2 > last( this->weight->lbound(firstDim), j );
						positions.push_back( last );
						distances.push_back( 0 );
					}
				}
				break;

			case secondDim:

				for ( int i = 0; i < this->weight->rows(); i++ ){
					if ( (*this->weight)( i, this->weight->lbound(secondDim) ) < numeric_limits<Pixel>::max() ){
						TinyVector< int, 2 > last( i, this->weight->lbound(secondDim) );
						positions.push_back( last );
						distances.push_back( 0 );
					}
				}
	}

	// set alive points
	fastMarching.setAliveSet( positions, distances );

	Array< Pixel, 2 > distance( this->weight->shape() );

	// run fast marching (no stopping condition)
	fastMarching.execute( distance );

	TinyVector< int, 2 > start, end;

    // compute first arriving point to the oposite side (start)
	// and corresponding circular neighbor (end)
	switch ( this->orientation ){
			case firstDim:
				start(firstDim) = this->weight->ubound(firstDim);
				start(secondDim) = ( minIndex( distance( this->weight->ubound(firstDim), Range::all() ) ))(firstDim);
				end(firstDim) = this->weight->lbound(firstDim);
				end(secondDim) = start(secondDim);
				break;
			case secondDim:
				start(firstDim) = ( minIndex( distance( Range::all(), this->weight->ubound(secondDim) ) ))(firstDim);
				start(secondDim) = this->weight->ubound(secondDim);
				end(firstDim) = start(firstDim);
				end(secondDim) = this->weight->lbound(secondDim);
	}

	// back propagation from start until we reach a point on the other side
	// (not neccesarily end).
	this->getPath( distance, start, end, path, true );

	// if last point coincides with end, we are done
	if ( sum( fabs( path[ path.size() - 1 ] - end ) ) == 0 ){
		// compute implicit path representation
		this->getImplicitPath( path, ipath );
	}
	else{
		switch ( this->orientation ){
			case firstDim:
				start(secondDim) = path[ path.size() - 1 ](secondDim);
				end(secondDim) = start(secondDim);
				break;
			case secondDim:
				start(firstDim) = path[ path.size() - 1 ](firstDim);
				end(firstDim) = start(firstDim);
		}
		this->getPath( distance, start, end, path, true );
		this->getImplicitPath( path, ipath );
	}
	//else{
	//	// we need a strictly circular path

	//	// re-set start point, so it is the first to be reached by the geodesic
	//	// re-set end point, so the path is circular
	//	switch ( this->orientation ){
	//		case firstDim:
	//			start(secondDim) = path[ path.size() - 1 ](secondDim);
	//			end(secondDim) = start(secondDim);
	//			break;
	//		case secondDim:
	//			start(firstDim) = path[ path.size() - 1 ](firstDim);
	//			end(firstDim) = start(firstDim);
	//	}

	//	// re-compute distances so we make sure the geodesic goes from start to end
	//	positions.clear();
	//	distances.clear();

	//	// use start as alive set
	//	positions.push_back( start );
	//	distances.push_back( 0 );
	//	fastMarching.setAliveSet( positions, distances );

	//	// only compute until we reach end
	//	fastMarching.setStopPoint( end );

	//	fastMarching.execute(distance);

	//	// back propagation from end to start
	//	this->getPath( distance, end, start, path, ipath );
	//}
}

template< class Pixel >
void nbfGeodesicPath< Pixel > :: getForwardCircularPath( vector< TinyVector< int, 2 > > & path, Array< Pixel, 2 > & ipath )
{
	// find a circular path between parallel sides of a 2d image.
	// The circular constraint is enforced by restricting the start and
	// end points of the path to be circular neighbors, i.e. to have the
	// same coordinate in the appropriate dimension.

	TinyVector< int, 2 > start;

    // compute first arriving point to the opposite side (start)
	switch ( this->orientation ){
			case firstDim:
				start(firstDim) = this->weight->ubound(firstDim);
				start(secondDim) = ( minIndex( (*this->weight)( this->weight->ubound(firstDim), Range::all() ) ))(firstDim);
				break;
			case secondDim:
				start(firstDim) = ( minIndex( (*this->weight)( Range::all(), this->weight->ubound(secondDim) ) ))(firstDim);
				start(secondDim) = this->weight->ubound(secondDim);
	}

	if ( ( start(secondDim) < 0 ) || ( start(firstDim) < 0 ) ){
		cout << "WARNING: cannot find starting point." << start << endl;
		return;
	}

	// back propagation from start until we reach a point on the other side
	// (not neccesarily end).
	this->getForwardPath( *this->weight, start, path );

	switch ( this->orientation ){
			case firstDim:
				// if last point "coincides" with end, we are done
				if ( abs( path[ path.size() - 1 ](secondDim) - start(secondDim) ) > 1 ){
					start(secondDim) = path[ path.size() - 1 ](secondDim);
					this->getForwardPath( (*this->weight), start, path );
				}
				break;
			case secondDim:
				if ( abs( path[ path.size() - 1 ](firstDim) - start(firstDim) ) > 1 ){
					start(firstDim) = path[ path.size() - 1 ](firstDim);
					this->getForwardPath( (*this->weight), start, path );
				}
	}
	this->getImplicitPath( path, ipath );
}

template< class Pixel >
void nbfGeodesicPath< Pixel> :: getPath( Array< Pixel, 2 > & distance,
										   TinyVector< int, 2 > start, 
										   TinyVector< int, 2 > end, 
										   vector< TinyVector< int, 2 > > & path,
										   bool reachOtherSide )
{	
	this->weight = &distance;

	path.clear();

	TinyVector< int, 2 > next;
	vector< TinyVector< int, 2 > > neighbors;

	// back propagation - minimum distance neighbor
	next = start;
	path.push_back( next );

	Pixel minDistance;

	while ( ( sum( fabs( next - end ) ) != 0 ) && 
		( ( !reachOtherSide ) || ( next(this->orientation) != end(this->orientation) ) ) ){	

			bool hasMultipleMinima;
			hasMultipleMinima = false;

			this->getFullNeighbors( next, neighbors );

			TinyVector< int, 2 > current = next;

			//cout << next << endl;
			//cout << distance( Range( next(firstDim) - 1, next(firstDim) + 1 ),
			//	Range( next(secondDim) - 1, next(secondDim) + 1 ) ) << endl;

			if ( neighbors.size() > 0 ){
				minDistance = distance(neighbors[0]);
				next = neighbors[0];
				for ( int i = 0; i < neighbors.size(); i++ ){
					if ( distance( neighbors[i] ) < minDistance ){
						minDistance = distance( neighbors[i] );
						next = neighbors[i];
						hasMultipleMinima = false;
					}
					else if ( distance( neighbors[i] ) == minDistance ){ 
						hasMultipleMinima = true;
					}
				}

				if ( hasMultipleMinima == true ){
					Pixel minD = numeric_limits<Pixel>::max();
					Pixel minIndex = 0;
					for ( int i = 0; i < neighbors.size(); i++ ){
						Pixel d = pow2( neighbors[i](firstDim) - current(firstDim) ) + pow2( neighbors[i](secondDim) - current(secondDim) );
						d = ( distance( neighbors[i] ) - distance( current ) ) / sqrt( d );
						if ( d < minD ){
							minD = d;
							minIndex = i;
						}
					}
					next = neighbors[minIndex];
				}

				//cout << "multiple minima " << hasMultipleMinima << endl;
				//cout << "d = " << distance( next ) << endl;

				// detect self intersections
				if ( ( path.size() > 1000 ) || ( sum( abs( path[ path.size() - 1 ] - next ) ) == 0 ) ){
					cout << "back propagation stalled. Path size = " << path.size() << endl;
					cout << path[ path.size() - 1 ] << ", d = " << distance( path[ path.size() - 1 ] ) << endl;
					cout << next << ", d = " << distance( next ) << endl;
					//cout << distance( Range( next(firstDim) - 1, next(firstDim) + 1 ), Range( next(secondDim) - 1, next(secondDim) + 1 ) )<< endl;
					break;
				}
				path.push_back( next );
				//distance( next ) = max( distance );
			}
			else{
				break;
			}
		}

	// add last point
	if ( sum( fabs( next - end ) ) == 0 )
		path.push_back( next );
}


template< class Pixel >
void nbfGeodesicPath< Pixel> :: getForwardPath( Array< Pixel, 2 > & distance,
												TinyVector< int, 2 > start,
												vector< TinyVector< int, 2 > > & path )
{	
	path.clear();

	TinyVector< int, 2 > next;

	// back propagation - minimum distance neighbor
	next = start;
	path.push_back( next );

	while ( ( ( this->orientation == firstDim ) && ( next(firstDim) != 0 ) ) ||
		( ( this->orientation == secondDim ) && ( next(secondDim) != 0 ) ) ){

			Pixel t1, t2, t3;

			if ( this->orientation == firstDim ){

				t1 = distance( next(firstDim) - 1, next(secondDim) - 1 );
				t2 = distance( next(firstDim) - 1, next(secondDim) );
				t3 = distance( next(firstDim) - 1, next(secondDim) + 1 );

				if ( ( t1 < t2 ) & ( t1 < t3 ) ){
					next(firstDim) = next(firstDim) - 1;
					next(secondDim) = next(secondDim) - 1;
				}
				else if ( ( t2 < t1 ) & ( t2 < t3 ) ){
					next(firstDim) = next(firstDim) - 1;
					next(secondDim) = next(secondDim);
				}
				else if ( ( t3 < t1 ) & ( t3 < t2 ) ){
					next(firstDim) = next(firstDim) - 1;
					next(secondDim) = next(secondDim) + 1;
				}
				else{
					next(firstDim) = next(firstDim) - 1;
					next(secondDim) = next(secondDim);
				}
			}
			else{
				t1 = distance( next(firstDim) - 1, next(secondDim) - 1 );
				t2 = distance( next(firstDim), next(secondDim) - 1 );
				t3 = distance( next(firstDim) + 1, next(secondDim) - 1 );

				if ( ( t1 < t2 ) & ( t1 < t3 ) ){
					next(firstDim) = next(firstDim) - 1;
					next(secondDim) = next(secondDim) - 1;
				}
				else if ( ( t2 < t1 ) & ( t2 < t3 ) ){
					next(firstDim) = next(firstDim);
					next(secondDim) = next(secondDim) - 1;
				}
				else if ( ( t3 < t1 ) & ( t3 < t2 ) ){
					next(firstDim) = next(firstDim) + 1;
					next(secondDim) = next(secondDim) - 1;
				}
				else{
					next(firstDim) = next(firstDim);
					next(secondDim) = next(secondDim) - 1;
				}

			}
			path.push_back( next );
		}

		// add last point
		path.push_back( next );
}


template< class Pixel >
void nbfGeodesicPath< Pixel > :: getPath( Array< Pixel, 2 > & distance,
											TinyVector< int, 2 > start, 
											TinyVector< int, 2 > end, 
											Array< Pixel, 2 > & ipath
											)
{
	vector< TinyVector< int, 2 > > path;	
    this->getPath( distance, start, end, path );
	this->getImplicitPath( path, ipath );
}

template< class Pixel >
void nbfGeodesicPath< Pixel > :: getPath( Array< Pixel, 2 > & distance,
											TinyVector< int, 2 > start, 
											TinyVector< int, 2 > end, 
											vector< TinyVector< int, 2 > > & path,
											Array< Pixel, 2 > & ipath
											)
{
	path.clear();	
    this->getPath( distance, start, end, path );
	this->getImplicitPath( path, ipath );
}
													 
template< class Pixel >
void nbfGeodesicPath< Pixel > :: getImplicitPath( vector< TinyVector< int, 2 > > & path,
													 Array< Pixel, 2 > & A )
{
	if ( sum( fabs( A.shape() - this->weight->shape() ) ) != 0 ){
		A.resize( this->weight->shape() );
	}
	A = 10 * blitz::minmax::max( A.rows(), A.cols() );

	switch ( this->orientation ){

		case firstDim:
			for ( int i = 0; i < path.size(); i++ ){
				Pixel x = path[i](firstDim);
				Pixel y = path[i](secondDim);

				firstIndex si;
				A( x, Range( fromStart, y ) ) = where( A( x, Range( fromStart, y ) ) > si - y, si - y, A( x, Range( fromStart, y ) ) );
				A( x, Range( y + 1, toEnd ) )= where( A( x, Range( y + 1, toEnd ) ) > si, si + 1, A( x, Range( y + 1, toEnd ) ) );
			}
			break;

		case secondDim:
			for ( int i = 0; i < path.size(); i++ ){
				Pixel x = path[i](firstDim);
				Pixel y = path[i](secondDim);

				firstIndex fi;
				A( Range( fromStart, x), y ) = where( A( Range( fromStart, x), y ) > fi - x, fi - x, A( Range( fromStart, x), y ) );
				A( Range( x + 1, toEnd ), y ) = where( A( Range( x + 1, toEnd ), y ) > fi, fi + 1, A( Range( x + 1, toEnd ), y ) );
			}
	}
}

#endif /* FILE_nbfGeodesicPath */