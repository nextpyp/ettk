#ifndef FILE_nbfPolarDomain
#define FILE_nbfPolarDomain

#include <vtkMath.h>

// Class nbfPolarDomain.
//
// Implements 2D and 3D coordinate polar transformation.
// 
// example:
//
//		Array<float,3> A,P,Q;
//		nbfPolarDomain< float, 3 > polar;
//		polar.setCenter( TinyVector<float,3>(2,3,4) );
//
//		// compute polar transformation
//		polar.cartesian2polar(A,P);
//
//		// get back cartesian coordinates
//		polar.polar2cartesian(P,Q);

#include <nbfLinearInterpolator.h>

template< class Pixel, int const Dim >
class nbfPolarDomain
{
public:

	// constructor
	nbfPolarDomain();

	// destructor
	~nbfPolarDomain(){};

	// set center of polar coordinates
	void setCenter( TinyVector< Pixel, Dim > & );

	// restrict computation to: |rho| < M
	void setMaxRho( Pixel p ){ this->maxRho = p; }
	Pixel getMaxRho(){ return this->maxRho; }
	
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
	Pixel getResTheta(){ return ( this->maxTheta - this->minTheta ) / this->resTheta; }

	// restrict computation to: |phi| < M
	void setMaxPhi( Pixel p ){ this->maxPhi = p; }
	
	// restrict to |phi| > m
	void setMinPhi( Pixel p ){ this->minPhi = p; }

	// cartesian to polar conversion (2D)
	//void cartesian2polar( Array< Pixel, 2 > &, Array< Pixel, 2 > & );
	void cartesian2polar( Array< Pixel, 2 > &, Array< Pixel, 2 > &, Array< bool, 2 > & );

	// cartesian to polar conversion (3D)
	void cartesian2polar( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< bool, 3 > & );

	// polar to cartesian conversion (2D)
	void polar2cartesian( Array< Pixel, 2 > &, Array< Pixel, 2 > &, Array< bool, 2 > & );

	// polar to cartesian conversion (3D)
	void polar2cartesian( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< bool, 3 > &  );

	// set boundary condition
	void setBoundaryPoint( Array< Pixel, 2 > &, Pixel &, Pixel & );

private:

	// store shape of source array
	TinyVector< int, Dim > sourceShape;

	// store shape of polar array
	TinyVector< int, Dim > polarShape;

	// rho parameters
	Pixel minRho, maxRho, resRho;

	// theta parameters
	Pixel minTheta, maxTheta, resTheta;

	// phi parameters
	Pixel minPhi, maxPhi, resPhi;

	// polar center coordinates
	TinyVector< Pixel, Dim > center;
};

template< class Pixel, int const Dim >
nbfPolarDomain< Pixel, Dim > :: nbfPolarDomain()
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
void nbfPolarDomain< Pixel, Dim > :: setCenter( TinyVector< Pixel, Dim > & c )
{
	this->center = c;
}

template< class Pixel, int const Dim >
void nbfPolarDomain< Pixel, Dim > :: cartesian2polar( Array< Pixel, 2 > & carte,
													  Array< Pixel, 2 > & polar,
													  Array< bool, 2 >  & inside )
{
	// compute size of polar arrays
	// int rhoSize = floor( this->maxRho / this->resRho );
	int rhoSize = floor( ( this->maxRho - this->minRho ) / this->resRho );
	int thetaSize = floor( ( this->maxTheta - this->minTheta ) / this->resTheta );

	// store in attributes
	this->polarShape = TinyVector<int,2>( rhoSize, thetaSize );
	this->sourceShape = carte.shape();

	// temp rho and theta arrays
	Array< Pixel, 2 > rho( this->polarShape );
	Array< Pixel, 2 > theta( this->polarShape );

	// build rho array
	firstIndex i;
	// rho = this->maxRho / ( rhoSize - 1 ) * i;
	rho = ( this->maxRho - this->minRho ) / ( rhoSize - 1 ) * i + this->minRho;

	// build theta array
	secondIndex j;
	theta = ( this->maxTheta - this->minTheta ) / ( thetaSize - 1 ) * j + this->minTheta;

	// build interpolation coordinate X
	Array< Pixel, 2 > X( this->polarShape );
	X = this->center(firstDim) + rho * cos(theta);

	// build interpolation coordinate Y
	Array< Pixel, 2 > Y( this->polarShape );
	Y = this->center(secondDim) + rho * sin(theta);

	// interpolator
	nbfLinearInterpolator< Pixel, Dim > interpolator( carte );

	// interpolate, store in argument
	interpolator.interpolate( X, Y,	polar, inside );
}

template< class Pixel, int const Dim >
void nbfPolarDomain< Pixel, Dim > :: setBoundaryPoint( Array< Pixel, 2 > & P, Pixel & x, Pixel & y )
{
	// compute size of polar arrays
	Pixel currentRho = sqrt( ( this->center(firstDim) - x ) * ( this->center(firstDim) - x ) +
		( this->center(secondDim) - y ) * ( this->center(secondDim) - y ) );
	currentRho = currentRho / this->maxRho * ( this->polarShape(firstDim) - 1 );
	Pixel currentTheta = atan2( ( y - this->center(secondDim) ) , ( x - this->center(firstDim) ) );
	if ( currentTheta < 0 ){
		currentTheta = currentTheta + 2.0 * vtkMath::Pi();
	}
	currentTheta = currentTheta / this->maxTheta * ( this->polarShape(secondDim) - 1 );
	// check if inside domain
	if ( currentRho + 2 <= P.ubound(firstDim) ){
		P( Range(fromStart,currentRho-2), Range(floor(currentTheta),ceil(currentTheta)) ) = numeric_limits<Pixel>::max();
		P( Range(currentRho+2,toEnd), Range(floor(currentTheta),ceil(currentTheta) ) ) = numeric_limits<Pixel>::max();
	}
}

template< class Pixel, int const Dim >
void nbfPolarDomain< Pixel, Dim > :: polar2cartesian( Array< Pixel, 2 > & polar,
													  Array< Pixel, 2 > & carte,
													  Array< bool, 2 > & inside )
{
	// interpolator
	nbfLinearInterpolator< Pixel, Dim > interpolator( polar );

	// build output interpolation coordinate rho
	Array< Pixel, 2 > rhoF( this->sourceShape );

	firstIndex i;
	secondIndex j;

	// distances to center
	rhoF = sqrt( ( this->center(firstDim) - i ) * ( this->center(firstDim) - i ) +
		( this->center(secondDim) - j ) * ( this->center(secondDim) - j ) );

	// scale to maxRho
	rhoF = ( rhoF - this->minRho ) / ( this->maxRho - this->minRho ) * ( this->polarShape(firstDim) - 1 );

	// build output interpolation coordinate theta
	Array< Pixel, 2 > thetaF( this->sourceShape );

	// theta values [0,2*PI]
	thetaF = atan2( ( j - this->center(secondDim) ) , ( i - this->center(firstDim) ) );
	thetaF = where( thetaF < 0, thetaF + 2.0 * vtkMath::Pi(), thetaF );

	// scale to maxRho and maxTheta
	rhoF = rhoF / this->maxRho * ( this->polarShape(firstDim) - 1 );
	thetaF = thetaF / this->maxTheta * ( this->polarShape(secondDim) - 1 );

	// interpolate, store in argument
	interpolator.interpolate( rhoF, thetaF,	carte, inside );
}


template< class Pixel, int const Dim >
void nbfPolarDomain< Pixel, Dim > :: cartesian2polar( Array< Pixel, 3 > & carte,
													  Array< Pixel, 3 > & polar,
													  Array< bool, 3 > & inside )
{
	// compute size of polar arrays
	int rhoSize = floor( ( this->maxRho - this->minRho ) / this->resRho );
	int thetaSize = floor( ( this->maxTheta - this->minTheta ) / this->resTheta );
	int phiSize = floor( ( this->maxPhi - this->minPhi ) / this->resPhi );

	// store in attributes
	this->polarShape = TinyVector<int,3>(rhoSize,thetaSize,phiSize);
	this->sourceShape = carte.shape();

	// temp rho, theta and phi arrays
	Array< Pixel, 3 > rho( this->polarShape );
	Array< Pixel, 3 > theta( this->polarShape );
	Array< Pixel, 3 > phi( this->polarShape );

	// build rho array
	firstIndex i;
	rho = ( this->maxRho - this->minRho ) / ( rhoSize - 1 ) * i + this->minRho;

	// build theta array
	secondIndex j;
	theta = ( this->maxTheta - this->minTheta ) / ( thetaSize - 1 ) * j + this->minTheta;

	// build phi array
	thirdIndex k;
	phi = ( this->maxPhi - this->minPhi ) / ( phiSize - 1 ) * k + this->minPhi;

	// build interpolation coordinate X
	Array< Pixel, 3 > X( rho.shape() );
	X = this->center(firstDim) + rho * cos(theta) * sin(phi);

	// build interpolation coordinate Y
	Array< Pixel, 3 > Y( rho.shape() );
	Y = this->center(secondDim) + rho * sin(theta) * sin(phi);

	// build interpolation coordinate Z
	Array< Pixel, 3 > Z( rho.shape() );
	Z = this->center(thirdDim) + rho * cos(phi);

	// interpolator
	nbfLinearInterpolator< Pixel, Dim > interpolator( carte );

	// interpolate, store in argument
	interpolator.interpolate( X, Y,	Z, polar, inside );
}

template< class Pixel, int const Dim >
void nbfPolarDomain< Pixel, Dim > :: polar2cartesian( Array< Pixel, 3 > & polar,
													  Array< Pixel, 3 > & carte,
													  Array< bool, 3 > & inside )
{
	// interpolator
	nbfLinearInterpolator< Pixel, Dim > interpolator( polar );

	firstIndex i;
	secondIndex j;
	thirdIndex k;

	// build output interpolation coordinate rho
	Array< Pixel, 3 > rhoF( this->sourceShape );

	// distances to center
	rhoF = sqrt( ( this->center(firstDim) - i ) * ( this->center(firstDim) - i ) +
		( this->center(secondDim) - j ) * ( this->center(secondDim) - j ) +
		( this->center(thirdDim) - k ) * ( this->center(thirdDim) - k ) );

	// scale to maxRho
	rhoF = ( rhoF - this->minRho ) / ( this->maxRho - this->minRho ) * ( this->polarShape(firstDim) - 1 );
	
	// build output interpolation coordinate theta
	Array< Pixel, 3 > thetaF( this->sourceShape );

	// theta values [0,2*PI]
	thetaF = atan2( ( j - this->center(secondDim) ) , ( i - this->center(firstDim) ) );
	thetaF = where( thetaF < 0, thetaF + 2.0 * vtkMath::Pi(), thetaF );

	// scale to maxTheta
	thetaF = thetaF / this->maxTheta * ( this->polarShape(secondDim) - 1 );

	// build output interpolation coordinate phi
	Array< Pixel, 3 > phiF( this->sourceShape );

	// theta values [0,PI]
	phiF = sqrt( ( i - this->center(firstDim) ) * ( i - this->center(firstDim) )
		+ ( j - this->center(secondDim) ) * ( j - this->center(secondDim) ) );
	phiF = atan2( phiF , ( k - this->center(thirdDim) ) );

	// scale to maxPhi
	phiF = ( phiF - this->minPhi ) / ( this->maxPhi - this->minPhi ) * ( this->polarShape(thirdDim) - 1 );

	// interpolate, store in argument
	interpolator.interpolate( rhoF, thetaF,	phiF, carte, inside );
}

#endif /* FILE_nbfPolarDomain */
