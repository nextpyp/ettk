#ifndef FILE_nbfCylindricalDomain3
#define FILE_nbfCylindricalDomain3

// Class nbfcylindricalDomain.
//
// Implements 2D and 3D coordinate cylindrical transformation.
// 
// example:
//
//		Array<float,3> A,P,Q;
//		nbfcylindricalDomain< float, 3 > cylindrical;
//		cylindrical.setCenter( TinyVector<float,3>(2,3,4) );
//
//		// compute cylindrical transformation
//		cylindrical.cartesian2cylindrical(A,P);
//
//		// get back cartesian coordinates
//		cylindrical.cylindrical2cartesian(P,Q);

#include <nbfLinearInterpolator.h>

template< class Pixel >
class nbfCylindricalDomain3
{
public:

	// constructor
	nbfCylindricalDomain3();

	// destructor
	~nbfCylindricalDomain3(){};

	// set center of cylindrical coordinates
	void setCenter( TinyVector< int, 3 > & );

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

	// cartesian to cylindrical conversion (3D)
	void cartesian2cylindrical( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< bool, 3 > & );
	void cartesian2cylindricalFast( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< bool, 3 > &, int );

	// cylindrical to cartesian conversion (3D)
	void cylindrical2cartesian( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< bool, 3 > &  );

	void setBoundaryPoint( Array< Pixel, 3 > &, TinyVector< int, 3 > & );

private:

	// store shape of source array
	TinyVector< int, 3 > sourceShape;

	// store shape of cylindrical array
	TinyVector< int, 3 > cylindricalShape;

	// rho parameters
	Pixel minRho, maxRho, resRho;

	// theta parameters
	Pixel minTheta, maxTheta, resTheta;

	// phi parameters
	Pixel minPhi, maxPhi, resPhi;

	// z scale
	Pixel scale;

	// cylindrical center coordinates
	TinyVector< Pixel, 3 > center;
};

template< class Pixel >
nbfCylindricalDomain3< Pixel > :: nbfCylindricalDomain3()
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
void nbfCylindricalDomain3< Pixel > :: setCenter( TinyVector< int, 3 > & c )
{
	this->center = c;
}

template< class Pixel >
void nbfCylindricalDomain3< Pixel > :: cartesian2cylindrical( Array< Pixel, 3 > & carte,
												  Array< Pixel, 3 > & cylindrical,
												  Array< bool, 3 > & inside )
{
	// compute size of cylindrical arrays
	int rhoSize = floor( ( this->maxRho - this->minRho ) / this->resRho );
	int thetaSize = floor( ( this->maxTheta - this->minTheta ) / this->resTheta );
	int phiSize = floor( ( this->maxPhi - this->minPhi ) / this->resPhi );
	phiSize = carte.depth();

	// store in attributes
	this->cylindricalShape = TinyVector<int,3>(rhoSize,thetaSize,phiSize);
	this->sourceShape = carte.shape();

	// temp rho, theta and phi arrays
	Array< Pixel, 3 > rho( this->cylindricalShape );
	Array< Pixel, 3 > theta( this->cylindricalShape );
	//Array< Pixel, 3 > phi( this->cylindricalShape );

	// build rho array
	firstIndex i;
	rho = ( this->maxRho - this->minRho ) / ( rhoSize - 1 ) * i + this->minRho;

	// build theta array
	secondIndex j;
	theta = ( this->maxTheta - this->minTheta ) / ( thetaSize - 1 ) * j + this->minTheta;

	// build phi array
	thirdIndex k;
	//phi = ( this->maxPhi - this->minPhi ) / ( phiSize - 1 ) * k + this->minPhi;

	// build interpolation coordinate X
	Array< Pixel, 3 > X( rho.shape() );
	X = this->center(firstDim) + rho * cos(theta);

	// build interpolation coordinate Y
	Array< Pixel, 3 > Y( rho.shape() );
	//Y = this->center(secondDim) + rho * cos(phi);
	Y = this->center(secondDim) + rho * sin(theta);

	// build interpolation coordinate Z
	Array< Pixel, 3 > Z( rho.shape() );
	Z = k;
	//Z = this->center(thirdDim) + rho * sin(theta) * sin(phi) * this->scale;

	//cout << min(theta) << ", " << max(theta) << endl;
	//cout << min(rho) << ", " << max(rho) << endl;
	//cout << min(Z) << ", " << max(Z) << endl;

	// interpolator
	nbfLinearInterpolator< Pixel, 3 > interpolator( carte );

	//nbfMatlabWriter mwriter;
	//mwriter.setFileName("ipath");
	//mwriter.write(X);
	//mwriter.write(Y);
	//mwriter.write(Z);

	// interpolate, store in argument
	interpolator.interpolate( X, Y,	Z, cylindrical, inside );
}

template< class Pixel >
void nbfCylindricalDomain3< Pixel > :: cylindrical2cartesian( Array< Pixel, 3 > & cylindrical,
												  Array< Pixel, 3 > & carte,
												  Array< bool, 3 > & inside )
{
	// interpolator
	nbfLinearInterpolator< Pixel, 3 > interpolator( cylindrical );

	firstIndex i;
	secondIndex j;
	thirdIndex k;

	// build output interpolation coordinate rho
	Array< Pixel, 3 > rhoF( this->sourceShape );

	// distances to center
	rhoF = sqrt( ( this->center(firstDim) - i ) * ( this->center(firstDim) - i ) +
		         ( this->center(secondDim) - j ) * ( this->center(secondDim) - j ) );

	// scale to maxRho
	rhoF = ( rhoF - this->minRho ) / ( this->maxRho - this->minRho ) * ( this->cylindricalShape(firstDim) - 1 );
	
	// build output interpolation coordinate theta
	Array< Pixel, 3 > thetaF( this->sourceShape );

	// theta values [0,2*PI]
	thetaF = atan2( ( j - this->center(secondDim) ), ( i - this->center(firstDim) ) );
	thetaF = where( thetaF < 0, thetaF + 2.0 * vtkMath::Pi(), thetaF );

	// scale to maxTheta
	thetaF = thetaF / this->maxTheta * ( this->cylindricalShape(secondDim) - 1 );

	// build output interpolation coordinate phi
	Array< Pixel, 3 > phiF( this->sourceShape );

	// phi values [0,PI]
	phiF = k;

	// scale to maxPhi
	// phiF = ( phiF - this->minPhi ) / ( this->maxPhi - this->minPhi ) * ( this->cylindricalShape(thirdDim) - 1 );

	// w.write(phiF);

	// interpolate, store in argument
	interpolator.interpolate( rhoF, thetaF,	phiF, carte, inside );
}


#endif /* FILE_nbfcylindricalDomain3 */