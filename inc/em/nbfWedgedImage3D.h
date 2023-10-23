#pragma once

#include <vtkImageData.h>
#include <vtkTransform.h>
#include <vtkImageResample.h>
#include <vtkImageChangeInformation.h>

#include <em/nbfWedge3D.h>

#include <sstream>

#define NBF_WEDGED_SUB_IMAGE_3D     0
#define NBF_WEDGED_AVERAGE_IMAGE_3D 1

template< class Pixel > class nbfProjectionRotationMetric3D;

/** Interface for VTK-like input-output pipeline filters.
	Update state is kept internally so execution is only done when needed.
	User is responsible for changing the state when doing changes that affect the filter's output.
*/
template< class Pixel >
class nbfWedgedImage3D
{
public:

	nbfWedgedImage3D();
	virtual ~nbfWedgedImage3D();

	/// Get represented volume
	virtual void getImage( Array< Pixel, 3 > &, vtkTransform * = NULL, bool = false ) = 0;
	virtual void getImage( vtkImageData *, vtkTransform * = NULL, bool = false ) = 0;

	nbfWedge3D< Pixel > * getWedge(){ return &this->wedge; }

	virtual void getWedgeImage( Array< Pixel, 3 > &, vtkTransform * = NULL ) = 0;
	virtual void getWedgeImage( Array< Pixel, 2 > &, vtkTransform * = NULL ){};

	virtual void getWedgeImageHalf( Array< Pixel, 3 > &, vtkTransform * = NULL ) = 0;
	virtual void getWedgeImageHalf( Array< Pixel, 2 > &, vtkTransform * = NULL ){};

	virtual void getSphericalWedgeImage( Array< Pixel, 2 > &, vtkTransform * = NULL, nbfProjectionRotationMetric3D< Pixel > * = NULL ) = 0;

	virtual TinyVector< Pixel, 3 > getDimensions() = 0;

	// return object type
	virtual int getTypeId() = 0;

	// methods for serializing/unserializing object data
	virtual void unserialize( stringstream & ) = 0;
	virtual void serialize( stringstream & ) = 0;

	virtual bool isWedgeEffective(){ return this->wedge.isEffective(); }

	// wedge object
	nbfWedge3D< Pixel > wedge;

protected:

};


template< class Pixel >
nbfWedgedImage3D< Pixel > :: nbfWedgedImage3D()
{
}

template< class Pixel >
nbfWedgedImage3D< Pixel > :: ~nbfWedgedImage3D()
{
}