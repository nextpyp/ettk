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

template< class Pixel > class nbfProjectionRotationMetric3D;

using namespace blitz;

/** Wedge 3D. For efficiency, filtering is conceived as an all-in-one filter. 
	Normally we apply several different filters in cascade, which corresponds to multiplication
	in reciprocal space. This class allows activation and deactivation of several filters to be
	applied simultaneously. Supported filters are low and high pass, and for dealing with the wedge.
	Bandpass also supported.
	Two different wedges can be applied simultaneously wedge1 and wedge2.
*/
template< class Pixel >
class nbfWedge3D
{
public:

	nbfWedge3D();

	nbfWedge3D(const nbfWedge3D&);
	nbfWedge3D< Pixel > & operator= ( const nbfWedge3D< Pixel > & );
	
	virtual ~nbfWedge3D();

	/// Set wedge bounds [lwedge,uwedge] and transform t.
	void set( Pixel, vtkTransform * = NULL );
	void set( Pixel, Pixel, vtkTransform * = NULL );
	void set( TinyVector< Pixel, 2 > &, vtkTransform * = NULL );

	/// Reset wedge bounds to [-90,90].
	void reset();

	/// Get Fourier image of wedge after applying given transformation.
	void getImage( Array< Pixel, 3 > &, vtkTransform * = NULL );
	void getImage( vtkImageData *, vtkTransform * = NULL  );

	// Specify wedge by arbitrary array
	void setImage( Array< Pixel, 3 > & );

	// Specify wedge in spherical coordinates by arbitrary array
	void setSphericalImage( Array< Pixel, 2 > & );

	void getImageHalf( Array< Pixel, 3 > &, vtkTransform * = NULL );
	void getImageHalf( vtkImageData *, vtkTransform * = NULL  );

	/// Get Spherical Fourier image of wedge after applying given transformation.
	void getSphericalImage( Array< Pixel, 2 > &, vtkTransform * = NULL, nbfProjectionRotationMetric3D< Pixel > * = NULL );

	/// wedge limits
	Pixel lwedge, uwedge;

	Array< Pixel, 3 > wedge;
	Array< Pixel, 2 > sphericalWedge;

	static void smoothWedgeImage( Array< Pixel, 3 > & );

	// indicates if wedge is less than [-90,90]
	bool isEffective(){ return ( ( this->lwedge != -90 ) || ( this->uwedge != 90 ) ); }

protected:

	/// wedge transform
	vtkTransform * trans;

	// store canonical wedge geometry as plane normals
	double normal1[3];
	double normal2[3];
};


template< class Pixel >
nbfWedge3D< Pixel > :: nbfWedge3D()
: lwedge(-90), uwedge(90)
{ 
	this->trans = vtkTransform::New();
}

template< class Pixel >
nbfWedge3D< Pixel > :: nbfWedge3D( const nbfWedge3D & c )
{ 
	this->trans = vtkTransform::New();
	(*this) = c;
}

template< class Pixel >
nbfWedge3D< Pixel > :: ~nbfWedge3D()
{ 
	this->trans->Delete();
}

template< class Pixel >
nbfWedge3D< Pixel > & nbfWedge3D< Pixel > :: operator= ( const nbfWedge3D< Pixel > & param )
{
	this->lwedge = param.lwedge;
	this->uwedge = param.uwedge;
	for ( int i = 0; i < 3; i++ ){
		this->normal1[i] = param.normal1[i];
		this->normal2[i] = param.normal2[i];
	}

	this->trans->SetMatrix( param.trans->GetMatrix() );

	this->wedge.resize( param.wedge.shape() );
	this->wedge = param.wedge;

	this->sphericalWedge.resize( param.sphericalWedge.shape() );
	this->sphericalWedge = param.sphericalWedge;

	return *this;
}

template< class Pixel >
void nbfWedge3D< Pixel > :: set( Pixel wedge, vtkTransform * t )
{
	this->set( -wedge, wedge, t );
}

template< class Pixel >
void nbfWedge3D< Pixel > :: set( TinyVector< Pixel, 2 > & p, vtkTransform * t )
{
	this->set( p[0], p[1], t );
}

template< class Pixel >
void nbfWedge3D< Pixel > :: reset()
{
	this->set( -90, 90 );
}

template< class Pixel >
void nbfWedge3D< Pixel > :: set( Pixel lwedge, Pixel uwedge, vtkTransform * t )
{
	// check for valid wedge range
	if ( lwedge < -90 ){
		this->lwedge = -90;
	}
	else{
		this->lwedge = lwedge;
	}
	if ( uwedge > 90 ){
		this->uwedge = 90.0;
	}
	else{
		this->uwedge = uwedge;
	}

	// store canonical wedge position as plane normals
    /*
	this->normal1[0] =  vtkMath::RadiansFromDegrees( cos( - ( 90.0 - this->uwedge ) ) );
	this->normal1[1] =   0;
	this->normal1[2] =  vtkMath::RadiansFromDegrees( - sin( - ( 90.0 - this->uwedge ) ) );

	this->normal2[0] =  vtkMath::RadiansFromDegrees( cos( - ( 90.0 - this->lwedge ) ) );
	this->normal2[1] =   0;
	this->normal2[2] =  vtkMath::RadiansFromDegrees( - sin( - ( 90.0 - this->lwedge ) ) );
    */
	this->normal1[0] =  cos( - vtkMath::RadiansFromDegrees( 90.0 - this->uwedge ) );
	this->normal1[1] =   0;
	this->normal1[2] =  - sin( - vtkMath::RadiansFromDegrees( 90.0 - this->uwedge ) );

	this->normal2[0] =  cos( - vtkMath::RadiansFromDegrees( 90.0 - this->lwedge ) );
	this->normal2[1] =   0;
	this->normal2[2] =  - sin( - vtkMath::RadiansFromDegrees( 90.0 - this->lwedge ) );
	
    // store transform
	if ( t != NULL ){
		this->trans->SetMatrix( t->GetMatrix() );
	}

	this->wedge.free();
	this->sphericalWedge.free();
}

template< class Pixel >
void nbfWedge3D< Pixel > :: getImage( vtkImageData * data, vtkTransform * t )
{
	Array< Pixel, 3 > A;
	nbfVTKInterface :: vtkToBlitzReference( data, A );
	this->getImage( t, A );
}

template< class Pixel >
void nbfWedge3D< Pixel > :: getImage( Array< Pixel, 3 > & W, vtkTransform * t )
{
	if ( this->wedge.size() > 0 ){
        if ( sum( abs( W.shape() - this->wedge.shape() ) ) > 0 ){
			cerr << "ERROR - Specified wedge size does not match. In " __FILE__ "," << __LINE__ << endl;
		} else {

			// rotate wedge image (if neccessary)
			if ( t != NULL ){

                // build new transform with rotation component only
				vtkTransform * rotation = vtkTransform :: New();
				rotation->SetInput( t );
				double position[3];
				t->GetPosition(position);
				position[0] = - position[0];
				position[1] = - position[1];
				position[2] = - position[2];
				rotation->PostMultiply();
				rotation->Translate(position);

				vtkImageData * data = vtkImageData :: New();
				nbfVTKInterface :: blitzToVtk( this->wedge, data );

				vtkImageChangeInformation * change = vtkImageChangeInformation :: New();
				change->SetInput( data );
				change->SetOriginTranslation( - this->wedge.rows() / 2.0 , - this->wedge.cols() / 2.0, - this->wedge.depth() / 2.0 );

				vtkImageReslice * reslice = vtkImageReslice :: New();
				reslice->SetInput( change->GetOutput() );
				reslice->SetBackgroundLevel( max( this->wedge ) );
				reslice->SetInterpolationModeToCubic();
				reslice->SetResliceTransform( rotation );
				reslice->Update();

				nbfVTKInterface :: vtkToBlitz( reslice->GetOutput(), W );
				reslice->Delete();
				change->Delete();
				data->Delete();
				rotation->Delete();

			} else {
				W = this->wedge;
			}
		}
	} else if ( this->uwedge - this->lwedge >= 180 ){
		W = 1;
	} else {
        // apply rotational part of t to canonical vectors
		double vT1[3], vT2[3];

		// remove translation component
		vtkTransform * rotation = vtkTransform::New();
		rotation->SetInput( t );
		//vtkMatrix4x4 * mat = vtkMatrix4x4 :: New();
		//cout << *t->GetMatrix() << endl;
		//mat->DeepCopy( t->GetMatrix() );
		//rotation->SetMatrix( mat );
		//mat->Delete();
		double position[3];
		t->GetPosition(position);
		position[0] = - position[0];
		position[1] = - position[1];
		position[2] = - position[2];
		rotation->PostMultiply();
		rotation->Translate(position);

		rotation->TransformVector(this->normal1,vT1);
		rotation->TransformVector(this->normal2,vT2);
		rotation->Delete();

		// now generate filter in reciprocal space

		TinyVector< int, 3 > center = floor( W.shape() / 2.0 );
		firstIndex i; secondIndex j; thirdIndex k;

		W = ( ( vT1[0] * ( i - center[0] ) + 
			vT1[1] * ( j - center[1] ) + 
			vT1[2] * ( k - center[2] ) ) *
			( vT2[0] * ( i - center[0] ) + 
			vT2[1] * ( j - center[1] ) + 
			vT2[2] * ( k - center[2] ) ) <= 0.0 );

		this->smoothWedgeImage( W );
	}
}

template< class Pixel >
void nbfWedge3D< Pixel > :: setImage( Array< Pixel, 3 > & W )
{
	if ( min(W) == 0 ){
		this->wedge.resize( W.shape() );
		this->wedge = W;
	} else {
		this->reset();
	}
}

template< class Pixel >
void nbfWedge3D< Pixel > :: setSphericalImage( Array< Pixel, 2 > & W )
{
	if ( W.size() > 0 ){
		this->sphericalWedge.resize( W.shape() );
		this->sphericalWedge = W;
	} else {
		this->reset();
	}
}

template< class Pixel >
void nbfWedge3D< Pixel > :: getImageHalf( vtkImageData * data, vtkTransform * t )
{
	Array< Pixel, 3 > A;
	nbfVTKInterface :: vtkToBlitzReference( data, A );
	this->getImageHalf( t, A );
}

template< class Pixel >
void nbfWedge3D< Pixel > :: getImageHalf( Array< Pixel, 3 > & W, vtkTransform * t )
{
	W.resize( W.rows(), W.cols(), W.depth() / 2 + 1 );

	if ( this->uwedge - this->lwedge >= 180 ){
		W = 1;
	} else {
		// apply rotational part of t to canonical vectors
		double vT1[3], vT2[3];

		// remove translation component
		vtkTransform * rotation = vtkTransform::New();
		rotation->SetInput( t );
		double position[3];
		t->GetPosition(position);
		position[0] = - position[0];
		position[1] = - position[1];
		position[2] = - position[2];
		rotation->PostMultiply();
		rotation->Translate(position);

		rotation->TransformVector(this->normal1,vT1);
		rotation->TransformVector(this->normal2,vT2);
		rotation->Delete();

		// now generate filter in reciprocal space

		TinyVector< int, 3 > center( W.rows() / 2, W.cols() / 2, W.ubound(thirdDim) );
		firstIndex i; secondIndex j; thirdIndex k;

        /*
        cout << "Center = " << center[0] << ", " << center[1] << "," << center[2] << endl;
        cout << "normal1 = " << normal1[0] << ", " << normal1[1] << "," << normal1[2] << endl;
        cout << "normal2 = " << normal2[0] << ", " << normal2[1] << "," << normal2[2] << endl;
        cout << "vT1 = " << vT1[0] << ", " << vT1[1] << "," << vT1[2] << endl;
        cout << "vT2 = " << vT2[0] << ", " << vT2[1] << "," << vT2[2] << endl;
        */
		
        W = ( ( vT1[0] * ( i - center[0] ) + 
			    vT1[1] * ( j - center[1] ) + 
			    vT1[2] * ( k - center[2] ) ) *
			  ( vT2[0] * ( i - center[0] ) + 
			    vT2[1] * ( j - center[1] ) + 
			    vT2[2] * ( k - center[2] ) ) <= 0.0 );

		this->smoothWedgeImage( W );
	}
}


template< class Pixel >
void nbfWedge3D< Pixel > :: getSphericalImage( Array< Pixel, 2 > & W, vtkTransform * t, nbfProjectionRotationMetric3D< Pixel > * metric  )
{
	if ( this->sphericalWedge.size() > 0 ){
		if ( sum( abs( W.shape() - this->sphericalWedge.shape() ) ) > 0 ){
			cerr << "ERROR - Specified wedge size does not match. In " __FILE__ "," << __LINE__ << endl;
		} else {
			if ( t == NULL ){
				W = this->sphericalWedge;
			} else if ( metric == NULL ){
				cerr << "ERROR - Must specify metric object in order to do rotation: " << __FILE__ << "," << __LINE__ << endl;
			} else {
				metric->harmonics( this->sphericalWedge, t, W );
				//nbfMatlabWriter w;
				//w.setFileName("p.matlab");
				//w.write(W);
				// binarize to avoid spherical harmonics artifacts
				W = where( W < .5, 0, 1 );
				// metric->conditionSphericalWindow(W);
			}
		}
	} else if ( this->uwedge - this->lwedge >= 180 ){
		W = 1;
	} else {

		// apply t to canonical vectors
		double vT1[3], vT2[3];
		//vtkTransform * concat = vtkTransform::New();
		//if ( t != NULL ){
		//	concat->SetMatrix( t->GetMatrix() );
		//}
		//if ( this->trans != NULL ){
		//	concat->Concatenate( this->trans );
		//}

		// remove translation component
		vtkTransform * rotation = vtkTransform::New();
		rotation->SetInput( t );
		double position[3];
		t->GetPosition(position);
		position[0] = - position[0];
		position[1] = - position[1];
		position[2] = - position[2];
		rotation->PostMultiply();
		rotation->Translate(position);

		rotation->TransformVector(this->normal1,vT1);
		rotation->TransformVector(this->normal2,vT2);
		rotation->Delete();

		// now generate filter in reciprocal space
		firstIndex i; secondIndex j;

		Pixel phiK = vtkMath::Pi() / W.ubound(firstDim);
		Pixel thetaK = 2.0 * vtkMath::Pi() / W.ubound(secondDim);

		W = ( ( vT1[0] * cos( j * thetaK ) * sin( i * phiK ) + vT1[1] * sin( j * thetaK ) * sin( i * phiK ) + vT1[2] * cos( i * phiK ) ) *
			( vT2[0] * cos( j * thetaK ) * sin( i * phiK ) + vT2[1] * sin( j * thetaK ) * sin( i * phiK ) + vT2[2] * cos( i * phiK ) ) <= 0.0 );

		//double axis[3];
		//vtkMath::Cross(vT1,vT2,axis);
		//vtkMath::Normalize(axis);
		//W = where( abs( axis[0] * cos( j * thetaK ) * sin( i * phiK ) + 
		//			    axis[1] * sin( j * thetaK ) * sin( i * phiK ) + 
		//				axis[2] * cos( i * phiK ) ) > .9, 0, W ); 
	}
}

template< class Pixel >
void nbfWedge3D< Pixel > :: smoothWedgeImage( Array< Pixel, 3 > & W )
{
	vtkImageData * window = vtkImageData::New();
	vtkImageGaussianSmooth * filter = vtkImageGaussianSmooth::New();

	// filter down the window functions (to improve spherical harmonics representation)
	nbfVTKInterface::blitzToVtk( W, window );
	filter->SetInput( window );
	filter->SetRadiusFactors(.5,.5,.5);
	filter->Update();
	nbfVTKInterface::vtkToBlitz( filter->GetOutput(), W );
	filter->Delete();
	window->Delete();
}
