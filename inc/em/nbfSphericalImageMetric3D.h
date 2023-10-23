#pragma once

#include <io/nbfVTKInterface.h>
#include <em/nbfCorrelationImageMetric3D.h>
#include <em/nbfCorrelationImageMetric2D.h>
#include <em/nbfFourierFilter.h>
#include <nbfPolarDomain.h>

#include <nbfTimer.h>

using namespace blitz;

/** 3D volume alignment based on cross-correlation. Given two input volumes, a brute force search 
	is performed with rotated versions of the second image.
	Accounting for the wedge is supported by specifying a nbfWedge object with setting for both volumes.
	Search domain bounds and sampling rate are specified with nbfProjectionRotationMetric3D::setAngleSearch#.
	Search can be done on downsized volumes by specifying a scale factor.
	@todo Use more sophisticated search strategies to speed up brute force approach (manage multiple scales, 
	different sampling rates, etc).
	@see nbfCorrelationMetricFFTW
*/
template< class Pixel >
class nbfSphericalImageMetric3D : public nbfCorrelationImageMetric< Pixel, 3 >
{
public:

	/// Constructor.
	nbfSphericalImageMetric3D( nbfWedgeFilter3D< Pixel > * = NULL );

	~nbfSphericalImageMetric3D();

	/// Redefine from nbfImageMetric
	void execute();

protected:

	nbfCorrelationImageMetric3D< Pixel > * ccc;
};

template< class Pixel >
nbfSphericalImageMetric3D< Pixel > :: nbfSphericalImageMetric3D( nbfWedgeFilter3D< Pixel > * w )
: nbfCorrelationImageMetric< Pixel, 3 >( w )
{
}

template< class Pixel >
nbfSphericalImageMetric3D< Pixel > :: ~nbfSphericalImageMetric3D()
{
}

template< class Pixel >
void nbfSphericalImageMetric3D< Pixel > :: execute()
{
	this->initialize();

	// apply combined wedge to both images (in-place)
	this->filter->applyReal( this->input1, this->input1 );
	this->filter->applyReal( this->input2, this->input2 );

	// create blitz references to image data
	Array< Pixel, 3 > in1, in2;
	nbfVTKInterface::vtkToBlitzReference( this->input1, in1 );
	nbfVTKInterface::vtkToBlitzReference( this->input2, in2 );

	nbfMatlabWriter w;
	w.setFileName("p.matlab");
	w.write(in1);
	w.write(in2);

	// transform to spherical coordinates

	Array< Pixel, 3 > polar1, polar2;
	Array< bool, 3 > B;

	nbfPolarDomain< Pixel, 3 > polar;
	TinyVector< Pixel, 3 > center = ( in1.shape() - 1 ) / 2.0;

	polar.setCenter( center );
	int resRho = max( center[0], max( center[1], center[2] ) );
	// make size even for FFTW
	if ( resRho % 2 != 0 ){
		resRho++;
	}
	polar.setMaxRho( resRho );
	polar.setResRho( resRho );
	polar.setResTheta( 180 );
	
	polar.cartesian2polar( in1, polar1, B );
	polar.cartesian2polar( in2, polar2, B );
	
	w.write(polar1);
	w.write(polar2);

	//P( Range(0,this->innerRadius), Range::all() ) = numeric_limits<float>::max();

	nbfFourierFilter< Pixel, 3 > filter;
	//filter.bandPassOn(.02,.5,.01);

	nbfCorrelationImageMetric3D< Pixel > ccc( &filter );
	ccc.windowOff();
	// ccc.windowPolarOn();

	vtkImageData * data1 = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( polar1, data1 );

	vtkImageData * data2 = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( polar2, data2 );

	ccc.setInput1( data1 );
	ccc.setInput2( data2 );
	ccc.execute();

	this->correlationPeak = ccc.getCorrelationPeak();
	
	Pixel pos[3];
	ccc.getTransform()->GetPosition(pos);

	cout << "theta = " << pos[1]*2  << ", phi " << pos[2]*2 << endl;

	this->transform->SetMatrix( ccc.getTransform()->GetMatrix() );

	data1->Delete();
	data2->Delete();
}