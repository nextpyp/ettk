#ifndef FILE_nbfPolarDomain3
#define FILE_nbfPolarDomain3

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

template< class Pixel >
class nbfPolarDomain3
{
public:

	// constructor
	nbfPolarDomain3();

	// destructor
	~nbfPolarDomain3(){};

	// set center of polar coordinates
	void setCenter( TinyVector< Pixel, 3 > & );

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

	// set z scale
	void setZScale( Pixel p )
	{ 
		this->scale = p; 
	}

	// restrict computation to: |phi| < M
	void setMaxPhi( Pixel p ){ this->maxPhi = p; }
	
	// restrict to |phi| > m
	void setMinPhi( Pixel p ){ this->minPhi = p; }

	// cartesian to polar conversion (3D)
	void cartesian2polar( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< bool, 3 > & );
	void cartesian2polarFast( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< bool, 3 > &, int );

	// polar to cartesian conversion (3D)
	void polar2cartesian( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< bool, 3 > &  );

	void setBoundaryPoint( Array< Pixel, 3 > &, TinyVector< int, 3 > & );

private:

	// store shape of source array
	TinyVector< int, 3 > sourceShape;

	// store shape of polar array
	TinyVector< int, 3 > polarShape;

	// rho parameters
	Pixel minRho, maxRho, resRho;

	// theta parameters
	Pixel minTheta, maxTheta, resTheta;

	// phi parameters
	Pixel minPhi, maxPhi, resPhi;

	// z scale
	Pixel scale;

	// polar center coordinates
	TinyVector< Pixel, 3 > center;
};

template< class Pixel >
nbfPolarDomain3< Pixel > :: nbfPolarDomain3()
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

	this->scale = 1.25;
}

template< class Pixel >
void nbfPolarDomain3< Pixel > :: setCenter( TinyVector< Pixel, 3 > & c )
{
	this->center = c;
}

template< class Pixel >
void nbfPolarDomain3< Pixel > :: cartesian2polar( Array< Pixel, 3 > & carte,
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
	//Y = this->center(secondDim) + rho * cos(phi);
	Y = this->center(secondDim) + rho * sin(theta) * sin(phi) * this->scale;

	// build interpolation coordinate Z
	Array< Pixel, 3 > Z( rho.shape() );
	Z = this->center(thirdDim) + rho * cos(phi);
	//Z = this->center(thirdDim) + rho * sin(theta) * sin(phi) * this->scale;

	// interpolator
	nbfLinearInterpolator< Pixel, 3 > interpolator( carte );

	// interpolate, store in argument
	interpolator.interpolate( X, Y,	Z, polar, inside );
}

template< class Pixel >
void nbfPolarDomain3< Pixel > :: polar2cartesian( Array< Pixel, 3 > & polar,
												  Array< Pixel, 3 > & carte,
												  Array< bool, 3 > & inside )
{
	// interpolator
	nbfLinearInterpolator< Pixel, 3 > interpolator( polar );

	firstIndex i;
	secondIndex j;
	thirdIndex k;

	// build output interpolation coordinate rho
	Array< Pixel, 3 > rhoF( this->sourceShape );

	// distances to center
	rhoF = sqrt( ( this->center(firstDim) - i ) * ( this->center(firstDim) - i ) +
		( this->center(secondDim) - j ) / this->scale * ( this->center(secondDim) - j ) / this->scale +
		( this->center(thirdDim) - k ) * ( this->center(thirdDim) - k ) );

	// scale to maxRho
	rhoF = ( rhoF - this->minRho ) / ( this->maxRho - this->minRho ) * ( this->polarShape(firstDim) - 1 );
	
	// build output interpolation coordinate theta
	Array< Pixel, 3 > thetaF( this->sourceShape );

	// theta values [0,2*PI]
	thetaF = atan2( ( ( j - this->center(secondDim) ) / this->scale ), ( i - this->center(firstDim) ) );
	thetaF = where( thetaF < 0, thetaF + 2.0 * vtkMath::Pi(), thetaF );

	// scale to maxTheta
	thetaF = thetaF / this->maxTheta * ( this->polarShape(secondDim) - 1 );

	// build output interpolation coordinate phi
	Array< Pixel, 3 > phiF( this->sourceShape );

	// phi values [0,PI]
	phiF = sqrt( ( i - this->center(firstDim) ) * ( i - this->center(firstDim) )
		+ ( j - this->center(secondDim) ) / this->scale * ( j - this->center(secondDim) ) / this->scale );
	phiF = atan2( phiF , ( k - this->center(thirdDim) ) );

	// scale to maxPhi
	phiF = ( phiF - this->minPhi ) / ( this->maxPhi - this->minPhi ) * ( this->polarShape(thirdDim) - 1 );

	// interpolate, store in argument
	interpolator.interpolate( rhoF, thetaF,	phiF, carte, inside );
}


template< class Pixel >
void nbfPolarDomain3< Pixel > :: cartesian2polarFast( Array< Pixel, 3 > & carte,
												      Array< Pixel, 3 > & polar,
												      Array< bool, 3 > & inside,
													  int minRho )
{
	// compute size of polar arrays
	int rhoSize = floor( ( this->maxRho - this->minRho ) / this->resRho );
	int thetaSize = floor( ( this->maxTheta - this->minTheta ) / this->resTheta );
	int phiSize = floor( ( this->maxPhi - this->minPhi ) / this->resPhi );

	// store in attributes
	this->polarShape = TinyVector<int,3>(rhoSize,thetaSize,phiSize);
	this->sourceShape = carte.shape();

	inside.resize( this->polarShape );
	polar.resize( this->polarShape );

	// interpolator
	nbfLinearInterpolator< Pixel, 3 > interpolator( carte );
	Pixel rho, theta, phi, x, y, z;
	Pixel rhoC = ( this->maxRho - this->minRho ) / ( rhoSize - 1 );
	Pixel thetaC = ( this->maxTheta - this->minTheta ) / ( thetaSize - 1 );
	Pixel phiC = ( this->maxPhi - this->minPhi ) / ( phiSize - 1 );
	for ( int i = 0; i < rhoSize; i++ ){
		for ( int j = 0; j < thetaSize; j++ ){
			for ( int k = 0; k < phiSize; k++ ){
				rho = rhoC * (float)i + this->minRho;
				if ( rho < minRho ){
					polar(i,j,k) = numeric_limits< Pixel > :: max();
				}
				else{
					theta = thetaC * (float)j + this->minTheta;
					phi = phiC * (float)k + this->minPhi;
					x = this->center(firstDim) + rho * cos(theta) * sin(phi);
					y = this->center(secondDim) + rho * sin(theta) * sin(phi) * this->scale;
					z = this->center(thirdDim) + rho * cos(phi);
					polar(i,j,k) = interpolator.interpolateSingle(x,y,z,inside(i,j,k));
				}
			}
		}
	}
}

template< class Pixel >
void nbfPolarDomain3< Pixel > :: setBoundaryPoint( Array< Pixel, 3 > & P, TinyVector< int, 3 > & p )
{
	// compute size of polar arrays
	Pixel currentRho = sqrt( ( this->center(firstDim) - p[0] ) * ( this->center(firstDim) - p[0] ) +
		( this->center(secondDim) - p[1] ) * ( this->center(secondDim) - p[1] ) +
		( this->center(thirdDim) - p[2] ) * ( this->center(thirdDim) - p[2] ) );
	currentRho = ( currentRho - this->minRho ) / ( this->maxRho - this->minRho ) * ( this->polarShape(firstDim) - 1 );
	
	Pixel currentTheta = atan2( ( p[1] - this->center(secondDim) ) , ( p[0] - this->center(firstDim) ) );
	if ( currentTheta < 0 ){
		currentTheta = currentTheta + 2.0 * vtkMath::Pi();
	}
	currentTheta = currentTheta / this->maxTheta * ( this->polarShape(secondDim) - 1 );

	Pixel currentPhi = sqrt( ( p[0] - this->center(firstDim) ) * ( p[0] - this->center(firstDim) )
		+ ( p[1] - this->center(secondDim) ) * ( p[1] - this->center(secondDim) ) );
	currentPhi = atan2( currentPhi , ( p[2] - this->center(thirdDim) ) );
	currentPhi = ( currentPhi - this->minPhi ) / ( this->maxPhi - this->minPhi ) * ( this->polarShape(thirdDim) - 1 );
	
	Range J( floor(currentTheta), ceil(currentTheta) );
	Range K( floor(currentPhi), ceil(currentPhi) );
	P( Range(fromStart,currentRho-2), J, K ) = numeric_limits<Pixel>::max();
	P( Range(currentRho+2,toEnd), J, K ) = numeric_limits<Pixel>::max();
}


#endif /* FILE_nbfPolarDomain3 */