#ifndef FILE_nbfArray
#define FILE_nbfArray

// Class nbfArray.

#include <vector>
#include <vtkMath.h>

// The idea is to provide full class specializations for the different dimensions.
// As most of the code is generic, we place in a macro all the common code and
// add the specialized code at the end.

// This is the generic double template class (dummy).
template< class Pixel, int const Dim >
class nbfArray
{
};


// define macro for specializations with only common code

#define NBF_CLASS_ARRAY_MACRO(Dim) \
template< class Pixel > \
class nbfArray< Pixel, Dim > \
{ \
public: \
	nbfArray( Array< Pixel, Dim > & ); \
	nbfArray( Array< Pixel, Dim > &, TinyVector< int, Dim > & ); \
	~nbfArray(){}; \
	void getPatchNeighbors( TinyVector< int, 2 >, int, vector< TinyVector< int, 2 > > & ); \
    void getNeighbors( TinyVector< int, Dim > &, vector< TinyVector< int, Dim > > & ); \
	void getFullNeighbors( TinyVector< int, Dim > &, vector< TinyVector< int, Dim > > & ); \
	void restrictNeighbors( TinyVector< int, Dim > &, vector< TinyVector< int, Dim > > & ); \
	bool isAtCut( TinyVector< int, Dim > &, TinyVector< int, Dim > & ); \
	void getRho( Array< Pixel, Dim > & ); \
	void getRho3D( Array< Pixel, Dim > & ); \
	void getPhi( Array< Pixel, 3 > & ); \
	void getTheta( Array< Pixel, 2 > & ); \
	void setCutDimension( int d ){ this->cutDimension = d; } \
	int getCutDimension(){ return this->cutDimension; } \
	void setCenter( TinyVector< int, Dim > & c ){ this->center = c; } \
 \
protected: \
 \
	Array< Pixel, Dim > data; \
	bool hasCut; \
	TinyVector< int, Dim > center; \
	int cutDimension; \
 \
}; \
 \
template< class Pixel > \
nbfArray< Pixel, Dim > :: nbfArray( Array< Pixel, Dim > & data ) \
{ \
	this->data.reference( data ); \
	this->hasCut = false; \
	this->cutDimension = firstDim; \
} \
 \
 \
template< class Pixel > \
nbfArray< Pixel, Dim > :: nbfArray( Array< Pixel, Dim > & data, TinyVector< int, Dim > & center ) \
{ \
	this->data.reference( data ); \
	this->hasCut = true; \
	this->center = center; \
	this->cutDimension = firstDim; \
} \
 \
 \
template< class Pixel > \
void nbfArray< Pixel, Dim > :: restrictNeighbors( TinyVector< int, Dim > & currentPoint, vector< TinyVector< int, Dim > > & neighbors ) \
{ \
	vector< TinyVector< int, Dim > > :: iterator iter = neighbors.begin(); \
 \
	while ( iter != neighbors.end() ){ \
		if ( this->isAtCut( currentPoint, *iter ) ){ \
			iter = neighbors.erase( iter ); \
		} \
		else{ \
			++iter; \
		} \
	} \
} \
 \
 \
template< class Pixel > \
void nbfArray< Pixel, Dim > :: getPatchNeighbors( TinyVector< int, 2 > currentPoint, \
												  int patchSize, \
												  vector< TinyVector< int, 2 > > & neighbors ) \
{ \
	neighbors.clear(); \
 \
	int x = currentPoint(firstDim); \
	int y = currentPoint(secondDim); \
 \
	patchSize += 1; \
 \
	for ( int i = x - patchSize; i <= x + patchSize; i++ ){ \
		TinyVector< int, 2 > top( i, y - patchSize); \
		TinyVector< int, 2 > bottom( i, y + patchSize); \
		if ( this->data.isInRange( top ) ){ \
			neighbors.push_back( top ); \
		} \
		if ( this->data.isInRange( bottom ) ){ \
			neighbors.push_back( bottom ); \
		} \
	} \
 \
	for ( int j = y - patchSize; j <= y + patchSize; j++ ){ \
		TinyVector< int, 2 > left( x - patchSize, j); \
		TinyVector< int, 2 > right( x + patchSize, j); \
		if ( this->data.isInRange( left ) ){ \
			neighbors.push_back( left ); \
		} \
		if ( this->data.isInRange( right ) ){ \
			neighbors.push_back( right ); \
		} \
 \
	} \
} \
 \
 \
template< class Pixel > \
bool nbfArray< Pixel, Dim > :: isAtCut( TinyVector< int, Dim > & p1, TinyVector< int, Dim > & p2 ) \
{ \
	if ( this->cutDimension == firstDim ) { \
		if ( p1( secondDim ) < p2( secondDim ) ){ \
			if ( ( p1( secondDim ) == this->center( secondDim ) ) && \
				( p1( firstDim ) >= this->center( firstDim ) ) ){ \
				return true; \
			} \
			else { \
				return false; \
			} \
		} \
\
		if ( p2( secondDim ) < p1( secondDim ) ){ \
			if ( ( p2( secondDim ) == this->center( secondDim ) ) && \
				( p2( firstDim ) >= this->center( firstDim ) ) ){ \
				return true; \
			} \
			else { \
				return false; \
			} \
		} \
 \
		return false; \
	} \
	else \
	{ \
		if ( p1( firstDim ) < p2( firstDim ) ){ \
			if ( ( p1( firstDim ) == this->center( firstDim ) ) && \
				( p1( secondDim ) >= this->center( secondDim ) ) ){ \
				return true; \
			} \
			else { \
				return false; \
			} \
		} \
\
		if ( p2( firstDim ) < p1( firstDim ) ){ \
			if ( ( p2( firstDim ) == this->center( firstDim ) ) && \
				( p2( secondDim ) >= this->center( secondDim ) ) ){ \
				return true; \
			} \
			else { \
				return false; \
			} \
		} \
 \
		return false; \
	} \
} \
 \
 \
template< class Pixel > \
void nbfArray< Pixel, Dim > :: getRho( Array< Pixel, Dim > & rho ) \
{ \
	if ( this->hasCut ) \
	{ \
		rho.resize( this->data.shape() ); \
		firstIndex i; \
		secondIndex j; \
		rho = sqrt( pow2( i - this->center(firstDim) ) +  \
			pow2( j - this->center(secondDim) ) + 0.0 ); \
	} \
} \
 \
template< class Pixel > \
void nbfArray< Pixel, Dim > :: getRho3D( Array< Pixel, Dim > & rho ) \
{ \
	if ( this->hasCut ) \
	{ \
		rho.resize( this->data.shape() ); \
		firstIndex i; \
		secondIndex j; \
		thirdIndex k; \
		rho = sqrt( pow2( i - this->center(firstDim) ) +  \
			pow2( j - this->center(secondDim) ) + \
			pow2( k - this->center(thirdDim) ) + 0.0 ); \
	} \
} \
 \
template< class Pixel > \
void nbfArray< Pixel, Dim > :: getTheta( Array< Pixel, 2 > & theta ) \
{ \
	if ( this->hasCut ) \
	{ \
		firstIndex i; \
		secondIndex j; \
 \
		theta.resize( this->data.rows(), this->data.cols() ); \
		theta = atan2( j - this->center(secondDim) + 0.0, i - this->center(firstDim)); \
	} \
} \
 \
template< class Pixel > \
void nbfArray< Pixel, Dim > :: getPhi( Array< Pixel, 3 > & phi ) \
{ \
	if ( this->hasCut ) \
	{ \
 \
		firstIndex i; \
		secondIndex j; \
		thirdIndex k; \
 \
		phi.resize( this->data.shape() ); \
		Array< Pixel, Dim > rho( phi.shape() ); \
		this->getRho3D( rho ); \
		phi = acos( ( k - this->center(thirdDim) ) /   \
			sqrt( pow2( i - this->center(firstDim) ) +  \
			pow2( j - this->center(secondDim) ) + \
			pow2( k - this->center(thirdDim) ) + 0.0 ) ); \
	} \
} \

// expand macros for 2D and 3D

NBF_CLASS_ARRAY_MACRO(2)
NBF_CLASS_ARRAY_MACRO(3)

// expand specializations

template< class Pixel >
void nbfArray< Pixel, 2 > :: getNeighbors( TinyVector< int, 2 > & currentPoint,
		  								   vector< TinyVector< int, 2 > > & neighbors )
{
	// clear neighbors list
	neighbors.clear();

	// get current point coordinates
	int x = currentPoint(firstDim);
	int y = currentPoint(secondDim);

	// build 4-connected potential neighbors
	TinyVector< int, 2 > nfn(x-1,y);
	TinyVector< int, 2 > nfp(x+1,y);
	TinyVector< int, 2 > nsn(x,y-1);
	TinyVector< int, 2 > nsp(x,y+1);

	// if valid insert into neigbors list
	if ( this->data.isInRange( nfn ) )
		neighbors.push_back( nfn );
	if ( this->data.isInRange( nfp ) )
		neighbors.push_back( nfp );
	if ( this->data.isInRange( nsn ) )
		neighbors.push_back( nsn );
	if ( this->data.isInRange( nsp ) )
		neighbors.push_back( nsp );

	if ( this->hasCut )
	{
		this->restrictNeighbors( currentPoint, neighbors );
	}
}


template< class Pixel >
void nbfArray< Pixel, 3 > :: getNeighbors( TinyVector< int, 3 > & currentPoint,
										   vector< TinyVector< int, 3 > > & neighbors )
{
	// clear neighbors list
	neighbors.clear();

	// get current point coordinates
	int x = currentPoint(firstDim);
	int y = currentPoint(secondDim);
	int z = currentPoint(thirdDim);

	// build 6-connected potential neighbors
	TinyVector< int, 3 > nfn(x-1,y,z);
	TinyVector< int, 3 > nfp(x+1,y,z);
	TinyVector< int, 3 > nsn(x,y-1,z);
	TinyVector< int, 3 > nsp(x,y+1,z);
	TinyVector< int, 3 > ntn(x,y,z-1);
	TinyVector< int, 3 > ntp(x,y,z+1);

	// if valid coordinates, insert into neigbors list
	if ( this->data.isInRange( nfn(0), nfn(1), nfn(2) ) )
		neighbors.push_back( nfn );
	if ( this->data.isInRange( nfp(0), nfp(1), nfp(2) ) )
		neighbors.push_back( nfp );
	if ( this->data.isInRange( nsn(0), nsn(1), nsn(2) ) )
		neighbors.push_back( nsn );
	if ( this->data.isInRange( nsp(0), nsp(1), nsp(2) ) )
		neighbors.push_back( nsp );
	if ( this->data.isInRange( ntn(0), ntn(1), ntn(2) ) )
		neighbors.push_back( ntn );
	if ( this->data.isInRange( ntp(0), ntp(1), ntp(2) ) )
		neighbors.push_back( ntp );

	if ( this->hasCut )
	{
		this->restrictNeighbors( currentPoint, neighbors );
	}

}


template< class Pixel >
void nbfArray< Pixel, 2 > :: getFullNeighbors( TinyVector< int, 2 > & currentPoint,
											   vector< TinyVector< int, 2 > > & neighbors )
{
	// clear neighbors list
	neighbors.clear();

	// get current point coordinates
	int x = currentPoint(firstDim);
	int y = currentPoint(secondDim);

	// build 8-connected potential neighbors
	TinyVector< int, 2 > nN(x-1,y);
	TinyVector< int, 2 > nNE(x-1,y+1);
	TinyVector< int, 2 > nE(x,y+1);
	TinyVector< int, 2 > nSE(x+1,y+1);
	TinyVector< int, 2 > nS(x+1,y);
	TinyVector< int, 2 > nSW(x+1,y-1);
	TinyVector< int, 2 > nW(x,y-1);
	TinyVector< int, 2 > nNW(x-1,y-1);

    // TODO //
	// OPTIMIZE LOOKING AT THE CURRENT POINT COORDINATES

	// if valid coordinates insert into neigbors list
	if ( this->data.isInRange( nN ) )
		neighbors.push_back( nN );
	if ( this->data.isInRange( nNE ) )
		neighbors.push_back( nNE );
	if ( this->data.isInRange( nE ) )
		neighbors.push_back( nE );
	if ( this->data.isInRange( nSE ) )
		neighbors.push_back( nSE );
	if ( this->data.isInRange( nS ) )
		neighbors.push_back( nS );
	if ( this->data.isInRange( nSW ) )
		neighbors.push_back( nSW );
	if ( this->data.isInRange( nW ) )
		neighbors.push_back( nW );
	if ( this->data.isInRange( nNW ) )
		neighbors.push_back( nNW );

	if ( this->hasCut )
	{
		this->restrictNeighbors( currentPoint, neighbors );
	}

}


template< class Pixel >
void nbfArray< Pixel, 3 > :: getFullNeighbors( TinyVector< int, 3 > & currentPoint,
											   vector< TinyVector< int, 3 > > & neighbors )
{
	// clear neighbors list
	neighbors.clear();

	// get current point coordinates
	int x = currentPoint(firstDim);
	int y = currentPoint(secondDim);
	int z = currentPoint(thirdDim);

	// build 8-connected potential neighbors
	TinyVector< int, 3 > n1(x,y,z-1);
	TinyVector< int, 3 > n2(x,y,z+1);

	TinyVector< int, 3 > n3(x,y+1,z);
	TinyVector< int, 3 > n4(x,y+1,z-1);
	TinyVector< int, 3 > n5(x,y+1,z+1);

	TinyVector< int, 3 > n6(x,y-1,z);
	TinyVector< int, 3 > n7(x,y-1,z-1);
	TinyVector< int, 3 > n8(x,y-1,z+1);

	TinyVector< int, 3 > n9(x-1,y,z);
	TinyVector< int, 3 > n10(x-1,y,z-1);
	TinyVector< int, 3 > n11(x-1,y,z+1);

	TinyVector< int, 3 > n12(x-1,y+1,z);
	TinyVector< int, 3 > n13(x-1,y+1,z-1);
	TinyVector< int, 3 > n14(x-1,y+1,z+1);

	TinyVector< int, 3 > n15(x-1,y-1,z);
	TinyVector< int, 3 > n16(x-1,y-1,z-1);
	TinyVector< int, 3 > n17(x-1,y-1,z+1);

	TinyVector< int, 3 > n18(x+1,y,z);
	TinyVector< int, 3 > n19(x+1,y,z-1);
	TinyVector< int, 3 > n20(x+1,y,z+1);

	TinyVector< int, 3 > n21(x+1,y+1,z);
	TinyVector< int, 3 > n22(x+1,y+1,z-1);
	TinyVector< int, 3 > n23(x+1,y+1,z+1);

	TinyVector< int, 3 > n24(x+1,y-1,z);
	TinyVector< int, 3 > n25(x+1,y-1,z-1);
	TinyVector< int, 3 > n26(x+1,y-1,z+1);

    // TODO //
	// OPTIMIZE LOOKING AT THE CURRENT POINT COORDINATES

	// if valid coordinates insert into neigbors list
	if ( this->data.isInRange( n1 ) )
		neighbors.push_back( n1 );
	if ( this->data.isInRange( n2 ) )
		neighbors.push_back( n2 );
	if ( this->data.isInRange( n3 ) )
		neighbors.push_back( n3 );
	if ( this->data.isInRange( n4 ) )
		neighbors.push_back( n4 );
	if ( this->data.isInRange( n5 ) )
		neighbors.push_back( n5 );
	if ( this->data.isInRange( n6 ) )
		neighbors.push_back( n6 );
	if ( this->data.isInRange( n7 ) )
		neighbors.push_back( n7 );
	if ( this->data.isInRange( n8 ) )
		neighbors.push_back( n8 );
	if ( this->data.isInRange( n9 ) )
		neighbors.push_back( n9 );
	if ( this->data.isInRange( n10 ) )
		neighbors.push_back( n10 );
	if ( this->data.isInRange( n11 ) )
		neighbors.push_back( n11 );
	if ( this->data.isInRange( n12 ) )
		neighbors.push_back( n12 );
	if ( this->data.isInRange( n13 ) )
		neighbors.push_back( n13 );
	if ( this->data.isInRange( n14 ) )
		neighbors.push_back( n14 );
	if ( this->data.isInRange( n15 ) )
		neighbors.push_back( n15 );
	if ( this->data.isInRange( n16 ) )
		neighbors.push_back( n16 );
	if ( this->data.isInRange( n17 ) )
		neighbors.push_back( n17 );
	if ( this->data.isInRange( n18 ) )
		neighbors.push_back( n18 );
	if ( this->data.isInRange( n19 ) )
		neighbors.push_back( n19 );
	if ( this->data.isInRange( n20 ) )
		neighbors.push_back( n20 );
	if ( this->data.isInRange( n21 ) )
		neighbors.push_back( n21 );
	if ( this->data.isInRange( n22 ) )
		neighbors.push_back( n22 );
	if ( this->data.isInRange( n23 ) )
		neighbors.push_back( n23 );
	if ( this->data.isInRange( n24 ) )
		neighbors.push_back( n24 );
	if ( this->data.isInRange( n25 ) )
		neighbors.push_back( n25 );
	if ( this->data.isInRange( n26 ) )
		neighbors.push_back( n26 );

	if ( this->hasCut )
	{
		this->restrictNeighbors( currentPoint, neighbors );
	}

}


#endif /* FILE_nbfArray */