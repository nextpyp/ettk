#ifndef FILE_nbfSliceDomain
#define FILE_nbfSliceDomain

// Class nbfSliceDomain.
//
// Implements 3D slicing transformation.
// 
// example:
//
//		Array<float,3> A,P,Q;
//		nbfSliceDomain< float, 3 > slice;
//		slice.setCenter( TinyVector<float,3>(2,3,4) );
//
//		// compute Slice transformation
//		slice.cartesian2Slice(A,P);
//
//		// get back cartesian coordinates
//		slice.slice2cartesian(P,Q);

#include <nbfLinearInterpolator.h>

template< class Pixel, int const Dim >
class nbfSliceDomain
{
public:

	// constructor
	nbfSliceDomain();

	// destructor
	~nbfSliceDomain(){};

	// set center of Slice coordinates
	void setCenter( TinyVector< int, Dim > & );

	// restrict computation to: |rho| < M
	void setMaxRho( Pixel p ){ this->maxRho = p; }
	
	// restrict to |rho| > m
	void setMinRho( Pixel p ){ this->minRho = p; }

	// set resolution in rho ( default = 1/2 )
	void setResRho( Pixel p ){ this->resRho = ( this->maxRho - this->minRho ) / p; }

	// set angular resolution
	void setResTheta( Pixel p )
	{ 
		// use full resolution in theta [ 0, 2*PI ]
		this->resTheta = ( this->maxTheta - this->minTheta ) / p; 
		// use half resolution in phi [ 0, PI ]
		this->resPhi = ( this->maxPhi - this->minPhi ) / ( p / 2 ); 
	}

	// restrict computation to: |phi| < M
	void setMaxPhi( Pixel p ){ this->maxPhi = p; }
	
	// restrict to |phi| > m
	void setMinPhi( Pixel p ){ this->minPhi = p; }

	// cartesian to slice conversion
	// parameters: input, output, domain
	void slice( Array< Pixel, 3 > &, Array< Pixel, 3 > &, int, Array< bool, 3 > & );

	// Slice to cartesian conversion (3D)
	void unslice( Array< Pixel, 3 > &, Array< Pixel, 3 > &, int, Array< bool, 3 > &  );

private:

	// store shape of source array
	TinyVector< int, Dim > sourceShape;

	// store shape of Slice array
	TinyVector< int, Dim > sliceShape;

	// rho parameters
	Pixel minRho, maxRho, resRho;

	// theta parameters
	Pixel minTheta, maxTheta, resTheta;

	// phi parameters
	Pixel minPhi, maxPhi, resPhi;

	// Slice center coordinates
	TinyVector< Pixel, Dim > center;
};

template< class Pixel, int const Dim >
nbfSliceDomain< Pixel, Dim > :: nbfSliceDomain()
{
	// default values

	this->minRho = 0;
	this->maxRho = 100;
	this->resRho = 1;

	this->minTheta = 0;
	this->maxTheta = 2 * vtkMath::Pi();
	this->resTheta = ( this->maxTheta - this->minTheta ) / 100;

	this->minPhi = 0;
	this->maxPhi = vtkMath::Pi();
	this->resPhi = ( this->maxPhi - this->minPhi ) / ( 100 / 2 );
}

template< class Pixel, int const Dim >
void nbfSliceDomain< Pixel, Dim > :: setCenter( TinyVector< int, Dim > & c )
{
	this->center = c;
}

template< class Pixel, int const Dim >
void nbfSliceDomain< Pixel, Dim > :: slice( Array< Pixel, 3 > & carte,
									    	Array< Pixel, 3 > & slice,
											int sliceDimension,
											Array< bool,  3 > & inside )
{
	// compute shape of sliced array

	int sizeRho;
	int sizeTheta = 90;
	int sizeHeight = carte.extent(sliceDimension);

	Pixel d1, d2, d3, d4;

	switch ( sliceDimension ){
	case firstDim:
		// distance to four corners
		d1 = pow2( center[secondDim] - carte.lbound(secondDim) ) + pow2( center[thirdDim] - carte.lbound(thirdDim) );
		d2 = pow2( center[secondDim] - carte.lbound(secondDim) ) + pow2( center[thirdDim] - carte.ubound(thirdDim) );
		d3 = pow2( center[secondDim] - carte.ubound(secondDim) ) + pow2( center[thirdDim] - carte.lbound(thirdDim) );
		d4 = pow2( center[secondDim] - carte.ubound(secondDim) ) + pow2( center[thirdDim] - carte.ubound(thirdDim) );
		sizeRho = 2 * ceil( sqrt( max( d1, max( d2, max( d3, d4 ) ) ) ) ) + 1;
		break;
	case secondDim:
		// distance to four corners
		d1 = pow2( center[firstDim] - carte.lbound(firstDim) ) + pow2( center[thirdDim] - carte.lbound(thirdDim) );
		d2 = pow2( center[firstDim] - carte.lbound(firstDim) ) + pow2( center[thirdDim] - carte.ubound(thirdDim) );
		d3 = pow2( center[firstDim] - carte.ubound(firstDim) ) + pow2( center[thirdDim] - carte.lbound(thirdDim) );
		d4 = pow2( center[firstDim] - carte.ubound(firstDim) ) + pow2( center[thirdDim] - carte.ubound(thirdDim) );
		sizeRho = 2 * ceil( sqrt( max( d1, max( d2, max( d3, d4 ) ) ) ) ) + 1;
		break;
	case thirdDim:
		// distance to four corners
		d1 = pow2( center[firstDim] - carte.lbound(firstDim) ) + pow2( center[secondDim] - carte.lbound(secondDim) );
		d2 = pow2( center[firstDim] - carte.lbound(firstDim) ) + pow2( center[secondDim] - carte.ubound(secondDim) );
		d3 = pow2( center[firstDim] - carte.ubound(firstDim) ) + pow2( center[secondDim] - carte.lbound(secondDim) );
		d4 = pow2( center[firstDim] - carte.ubound(firstDim) ) + pow2( center[secondDim] - carte.ubound(secondDim) );
		sizeRho = 2 * ceil( sqrt( max( d1, max( d2, max( d3, d4 ) ) ) ) ) + 1;
		break;
	}

	// store in attributes
	this->sliceShape = TinyVector<int,3>(sizeRho,sizeTheta,sizeHeight);
	this->sourceShape = carte.shape();

	// temp rhos, and theta arrays
	Array< Pixel, 3 > rho( this->sliceShape );
	Array< Pixel, 3 > theta( this->sliceShape );
	Array< Pixel, 3 > height( this->sliceShape );

	firstIndex i;
	secondIndex j;
	thirdIndex k;

	rho = ( sizeRho - 1.0 ) / 2.0 - i;
	theta = j / ( sizeTheta - 1.0 ) * vtkMath::Pi();
	height = k;

	// build interpolation coordinates

	Array< Pixel, 3 > X( this->sliceShape );
	Array< Pixel, 3 > Y( this->sliceShape );
	Array< Pixel, 3 > Z( this->sliceShape );

	switch ( sliceDimension ){
	case firstDim:
		X = height;
		Y = rho * cos(theta) + this->center(secondDim);
		Z = rho * sin(theta) + this->center(thirdDim);
		break;
	case secondDim:
		X = rho * cos(theta) + this->center(firstDim);
		Y = height;
		Z = rho * sin(theta) + this->center(thirdDim);
		break;
	case thirdDim:
		X = rho * cos(theta) + this->center(firstDim);
		Y = rho * sin(theta) + this->center(secondDim);
		Z = height;
	}

	// interpolator
	nbfLinearInterpolator< Pixel, Dim > interpolator( carte );

	// interpolate, store in argument
	interpolator.interpolate( X, Y,	Z, slice, inside );
}

template< class Pixel, int const Dim >
void nbfSliceDomain< Pixel, Dim > :: unslice( Array< Pixel, 3 > & slice,
											  Array< Pixel, 3 > & carte,
											  int sliceDimension,
											  Array< bool, 3 > & inside )
{
	// interpolator
	nbfLinearInterpolator< Pixel, Dim > interpolator( slice );

	firstIndex i;
	secondIndex j;
	thirdIndex k;

	// build output interpolation coordinate rho
	Array< Pixel, 3 > rhoF( this->sourceShape );
	Array< Pixel, 3 > thetaF( this->sourceShape );
	Array< Pixel, 3 > heightF( this->sourceShape );

	switch ( sliceDimension ){
	case firstDim:
		thetaF = atan2( k - this->center(thirdDim), j - this->center(secondDim) );
		rhoF = sqrt( pow2( this->center(secondDim) - j ) +
		             pow2( this->center(thirdDim)  - k ) ) * where( thetaF < 0, 1, -1 );
		rhoF = rhoF + ( slice.rows() - 1 ) / 2.0;
		thetaF = where( thetaF < 0, thetaF + vtkMath::Pi(), thetaF ) / vtkMath::Pi() * ( slice.cols() - 1 );
		heightF = i;
		break;
	case secondDim:
		thetaF = atan2( k - this->center(thirdDim), i - this->center(firstDim) );
		rhoF = sqrt( pow2( this->center(firstDim) - i ) +
		             pow2( this->center(thirdDim)  - k ) ) * where( thetaF < 0, 1, -1 );
		rhoF = rhoF + ( slice.rows() - 1 ) / 2.0;
		thetaF = where( thetaF < 0, thetaF + vtkMath::Pi(), thetaF ) / vtkMath::Pi() * ( slice.cols() - 1 );
		heightF = j;
		break;
	case thirdDim:
		thetaF = atan2( j - this->center(secondDim), i - this->center(firstDim) );
		rhoF = sqrt( pow2( this->center(firstDim) - i ) +
		             pow2( this->center(secondDim) - j ) ) * where( thetaF < 0, 1, -1 );
		rhoF = rhoF + ( slice.rows() - 1 ) / 2;
		thetaF = where( thetaF < 0, thetaF + vtkMath::Pi(), thetaF ) / vtkMath::Pi() * ( slice.cols() - 1 );
		heightF = k;
	}

	// interpolate, store in argument
	interpolator.interpolate( rhoF, thetaF,	heightF, carte, inside );
}

#endif /* FILE_nbfSliceDomain */
