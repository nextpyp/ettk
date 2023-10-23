#pragma once

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkTransform.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkImageFourierCenter.h>

#include <io/nbfVTKInterface.h>
#include <em/nbfFourierFilter.h>
#include <em/nbfSubVolume.h>

using namespace blitz;

/** Fourier Filter 3D. For efficiency, filtering is conceived as an all-in-one filter. 
	Normally we apply several different filters in cascade, which corresponds to multiplication
	in reciprocal space. This class allows activation and deactivation of several filters to be
	applied simultaneously. Supported filters are low and high pass, and for dealing with the wedge.
	Bandpass also supported.
	Two different wedges can be applied simultaneously wedge1 and wedge2.
*/
template< class Pixel >
class nbfWedgeFilter3D : public nbfFourierFilter< Pixel, 3 >
{
public:

	nbfWedgeFilter3D();

	~nbfWedgeFilter3D();

	/// Set to use wedge 1 and/or 2 with bounds [+/-w] and transform t.
	void wedge1On( Pixel w, vtkTransform * t = NULL ){ this->wedge1On(-w,w,t); }
	void wedge2On( Pixel w, vtkTransform * t = NULL ){ this->wedge2On(-w,w,t); }
	
	/// Set to use wedge 1 and/or 2 with asymmetric wedge bounds [lwedge,uwedge] and transform t.
	void wedge1On( Pixel, Pixel, vtkTransform * = NULL );
	void wedge2On( Pixel, Pixel, vtkTransform * = NULL );

	void wedge1On( TinyVector< Pixel, 2 >, vtkTransform * = NULL );
	void wedge2On( TinyVector< Pixel, 2 >, vtkTransform * = NULL );

	void wedgeOn( nbfSubVolume< Pixel > );

	/// Specify additional transform for wedge 2. 
	/// This is provided for efficiency so wedge2 image is created when needed.
	void wedge2Transform( vtkTransform * = NULL );

	/// Switch off use of wedge1 and/or 2.
	void wedge1Off(){ this->useWedge1 = false; }
	void wedge2Off(){ this->useWedge2 = false; }

	/// Apply wedge2 *only* to reciprocal space input data. Arguments: (input,output)
	void applyComplexWedge2( Array< Pixel, 3 > &, Array< Pixel, 3 > & );

	Pixel getOverlap();

	/// Export wedge limits
	Pixel lwedge1, uwedge1;
	Pixel lwedge2, uwedge2;

	/// Export wedge transforms
	vtkTransform * t1;
	vtkTransform * t2;

	/// Build wedge image given wedge vectors and transform. Arguments: (vector1, vector2, wedge image).
	void buildSphericalWedge1( Array< Pixel, 2 > & );
	void buildSphericalWedge2( Array< Pixel, 2 > & );
	void buildSphericalWedge( double [3], double [3], vtkTransform *, Array< Pixel, 2 > & );

protected:

	// redefine from parent
	void updateCommonFilter();
	void updateNotCommonFilter();

	/// Build canonical wedge vectors given wedge bounds. Arguments: (lwedge, uwedge, result vector1, result vector2).
	void buildWedgeVectors( Pixel, Pixel, double [3], double [3] );
	
	/// Build wedge image given wedge vectors and transform. Arguments: (vector1, vector2, wedge image).
	void buildWedge( double [3], double [3], vtkTransform *, Array< Pixel, 3 > & );

	/// internal wedge representaions
	Array< Pixel, 3 > wedge1;
	Array< Pixel, 3 > wedge2;

	// store canonical wedge position as plane normals
	double v11[3];
	double v12[3];

	double v21[3];
	double v22[3];

	// store multiple wedges for composition
	vector< TinyVector< double, 3 > > wedge1vectors;
	vector< TinyVector< double, 3 > > wedge2vectors;

	// store multiple tranforms for composition
	vector< vtkTransform * > wedge1transforms;
	vector< vtkTransform * > wedge2transforms;

	// control execution of wedge2
	bool useWedge1, useWedge2, wedge2Updated;
};


template< class Pixel >
nbfWedgeFilter3D< Pixel > :: nbfWedgeFilter3D()
: nbfFourierFilter< Pixel, 3 >(),
  useWedge1(false), useWedge2(false), wedge2Updated(false), lwedge1(-90), uwedge1(90), lwedge2(-90), uwedge2(90)
{ 
	this->t1 = vtkTransform::New();
	this->t2 = vtkTransform::New();
}

template< class Pixel >
nbfWedgeFilter3D< Pixel > :: ~nbfWedgeFilter3D()
{ 
	this->t1->Delete();
	this->t2->Delete();
}


template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: updateCommonFilter()
{
	// call parent's
	nbfFourierFilter< Pixel, 3 > :: updateCommonFilter();

	// BYPASS
	//// apply wedge1
	//if ( this->useWedge1 == true ){
	//	this->buildWedge( this->v11, this->v12, this->t1, this->wedge1 );
	//	this->commonFilter *= this->wedge1;
	//}
}


template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: updateNotCommonFilter()
{
	nbfFourierFilter< Pixel, 3 > :: updateNotCommonFilter();

	// BYPASS
	//// apply wedge2
	//if ( this->useWedge2 == true ){
	//	// update wedge 2 if neccessary
	//	if ( this->wedge2Updated == false ){
	//		this->wedge2Transform();
	//	}
	//	this->notCommonFilter *= this->wedge2;
	//}
}

template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: buildWedgeVectors( Pixel lwedge, Pixel uwedge, double v1[3], double v2[3] )
{
	// check for valid wedge range
	if ( lwedge < -90 ){
		lwedge = -90;
	}
	if ( uwedge > 90 ){
		uwedge = 90.0;
	}

	// store canonical wedge position as plane normals

	v1[0] =   cos( ( 90.0 - uwedge ) * vtkMath::DegreesToRadians() );
	v1[1] =   0;
	v1[2] = - sin( ( 90.0 - uwedge ) * vtkMath::DegreesToRadians() );

	v2[0] =   cos( ( 90.0 - lwedge ) * vtkMath::DegreesToRadians() );
	v2[1] =   0;
	v2[2] = - sin( ( 90.0 - lwedge ) * vtkMath::DegreesToRadians() );
}


template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: wedge1On( Pixel lwedge, Pixel uwedge, vtkTransform * t )
{
	if ( ( lwedge != this->lwedge1 ) || ( uwedge != this->uwedge1 ) ){

		this->lwedge1 = lwedge;
		this->uwedge1 = uwedge;

		// build vectors first
		this->buildWedgeVectors( lwedge, uwedge, this->v11, this->v12 );
	}

	//this->wedge1vectors.clear();
	//this->wedge1vectors1.push_back( TinyVector< double, 3 >( this->v11[0], this->v11[1], this->v11[2] );
	//this->wedge1vectors2.push_back( TinyVector< double, 3 >( this->v12[0], this->v12[1], this->v12[2] );

	if ( t != NULL ){
		//this->t1->SetMatrix( t->GetMatrix() );
		vtkTransform * tmp = vtkTransform::New();
		this->t1->DeepCopy( tmp );
		tmp->Delete();
		this->t1->SetInput( t );
	}

	this->useWedge1 = true;
	this->commonFilterUpdated = false;
}

template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: wedge1On( TinyVector< Pixel, 2 > wedge, vtkTransform * t )
{
	this->wedge1On( wedge[0], wedge[1], t );
}

template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: wedge2On( Pixel lwedge, Pixel uwedge, vtkTransform * t )
{
	if ( ( lwedge != this->lwedge2 ) || ( uwedge != this->uwedge2 ) ){

		this->lwedge2 = lwedge;
		this->uwedge2 = uwedge;

		// build vectors
		this->buildWedgeVectors( lwedge, uwedge, this->v21, this->v22 );
	}

	// store transform
	if ( t != NULL ){
		//this->t2->SetMatrix( t->GetMatrix() );
		vtkTransform * tmp = vtkTransform::New();
		this->t2->DeepCopy( tmp );
		tmp->Delete();
		this->t2->SetInput( t );
	}

	// reset state
	this->useWedge2 = true;
	this->wedge2Updated = false;
	this->notCommonFilterUpdated = false;
}

template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: wedge2On( TinyVector< Pixel, 2 > wedge, vtkTransform * t )
{
	this->wedge2On( wedge[0], wedge[1], t );
}

template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: wedge2Transform( vtkTransform * t )
{
	// concatenate transforms
	vtkTransform * concat = vtkTransform::New();
	if ( t != NULL ){
		//concat->SetMatrix( t->GetMatrix() );
		vtkTransform * tmp = vtkTransform::New();
		concat->DeepCopy( tmp );
		tmp->Delete();
		concat->SetInput( t );
	}
	concat->Concatenate( this->t2 );

	// compute wedge2 image
	this->buildWedge( this->v21, this->v22, concat, this->wedge2 );
	concat->Delete();

	this->wedge2Updated = true;
	this->notCommonFilterUpdated = false;
}

template< class Pixel >
Pixel nbfWedgeFilter3D< Pixel > :: getOverlap()
{
	if ( this->useWedge1 && this->useWedge2 ){
		return sum( this->wedge1 * this->wedge2 ) / product( this->dimensions );
	}
	else{
		return 0;
	}
}

template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: buildWedge( double v1[3], double v2[3], vtkTransform * t, Array< Pixel, 3 > & W )
{
	// now generate filter in reciprocal space
	firstIndex i; secondIndex j; thirdIndex k;

	// apply t to canonical vectors
	double vT1[3], vT2[3];
	if ( t != NULL ){
		t->TransformVector(v1,vT1);
		t->TransformVector(v2,vT2);
	}
	else{
		vT1[0] = v1[0];
		vT1[1] = v1[1];
		vT1[2] = v1[2];
		vT2[0] = v2[0];
		vT2[1] = v2[1];
		vT2[2] = v2[2];
	}

	W.resize( dimensions[0], dimensions[1], dimensions[2] );

	W = ( ( vT1[0] * (i - this->center[0] ) + 
	        vT1[1] * (j - this->center[1] ) + 
	        vT1[2] * (k - this->center[2] ) ) *
	      ( vT2[0] * (i - this->center[0] ) + 
	        vT2[1] * (j - this->center[1] ) + 
	        vT2[2] * (k - this->center[2] ) ) <= 0.0 );
}


template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: buildSphericalWedge1( Array< Pixel, 2 > & W )
{
	if ( this->useWedge1 == true ){
		this->buildSphericalWedge( this->v11, this->v12, this->t1, W );
	}
	else{
		W = 1;
	}
}

template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: buildSphericalWedge2( Array< Pixel, 2 > & W )
{
	if ( this->useWedge2 == true ){
		this->buildSphericalWedge( this->v21, this->v22, this->t2, W );
	}
	else{
		W = 1;
	}
}

template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: buildSphericalWedge( double v1[3], double v2[3], vtkTransform * t, Array< Pixel, 2 > & W )
{
	// now generate filter in reciprocal space
	firstIndex i; secondIndex j;

	// apply t to canonical vectors
	double vT1[3], vT2[3];
	if ( t != NULL ){
		t->TransformVector(v1,vT1);
		t->TransformVector(v2,vT2);
	}
	else{
		vT1[0] = v1[0];
		vT1[1] = v1[1];
		vT1[2] = v1[2];
		vT2[0] = v2[0];
		vT2[1] = v2[1];
		vT2[2] = v2[2];
	}

	Pixel phiK = vtkMath::Pi() / W.ubound(firstDim);
	Pixel thetaK = 2.0 * vtkMath::Pi() / W.ubound(secondDim);

	W = ( ( vT1[0] * cos( j * thetaK ) * sin( i * phiK ) + vT1[1] * sin( j * thetaK ) * sin( i * phiK ) + vT1[2] * cos( i * phiK ) ) *
		  ( vT2[0] * cos( j * thetaK ) * sin( i * phiK ) + vT2[1] * sin( j * thetaK ) * sin( i * phiK ) + vT2[2] * cos( i * phiK ) ) <= 0.0 );

}

template< class Pixel >
void nbfWedgeFilter3D< Pixel > :: applyComplexWedge2( Array< Pixel, 3 > & real, Array< Pixel, 3 > & imag )
{
	if ( this->useWedge2 == false ){
		return;
	}
	else{
		// shift filter to match FFT geometry
		nbfVTKInterface::blitzToVtk( this->wedge2, this->filterData );

		centerF->SetInput( this->filterData );
		centerF->Modified();
		centerF->Update();
		nbfVTKInterface::vtkToBlitzReference( centerF->GetOutput(), this->shiftedFilter );
	}

	// apply filter
	real *= this->shiftedFilter;
	imag *= this->shiftedFilter;
}