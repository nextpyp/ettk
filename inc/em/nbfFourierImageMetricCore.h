#pragma once

#include <em/nbfFourierImageMetric.h>

using namespace blitz;

/** Compute normalized cross-correlation between two images in Fourier space.
	Input images may be in real or already in reciprocal representation.

	Correlation peak is located with sub-pixel accuracy, using parabolic fit (Frank p. 82) in each axis.

	@todo Improve peak search computation, i.e. compute two most significant peaks and 
	use difference for validation purposes. 
*/
template< class Pixel, int const Dim >
class nbfFourierImageMetricCore : public nbfFourierImageMetric< Pixel, Dim >
{
public:

	/// Constructor with wedge object argument (neccesary to carry out FFT computations).
	nbfFourierImageMetricCore( nbfImageFilter< Pixel, Dim > * = NULL, nbfFourierFilter< Pixel, Dim > * = NULL );

	virtual int getId(){ return NBF_IMAGE_METRIC_CORE; }

	/// Redefine from parent
	void execute();
};

template< class Pixel, int const Dim  >
nbfFourierImageMetricCore< Pixel, Dim > :: nbfFourierImageMetricCore( nbfImageFilter< Pixel, Dim > * i, nbfFourierFilter< Pixel, Dim  > * f )
: nbfFourierImageMetric< Pixel, Dim >(i,f)
{}

template< class Pixel, int const Dim >
void nbfFourierImageMetricCore< Pixel, Dim > :: execute()
{
	vtkTransform * t = vtkTransform::New();
	this->transform->DeepCopy(t);
	nbfFourierImageMetric< Pixel, Dim > :: executeFourierNewHalf(t);

	if ( this->isTransformValid( t ) == true ){
		this->transform->SetInput(t);
	}

	t->Delete();

	this->candidateTransforms.clear();
	this->candidateCorrelationPeaks.clear();
	this->candidateCorrelationScales.clear();
	this->candidateWedgeOverlaps.clear();

	Array< double, 1 > currentMatrix(16);

	double matrix[16];
	vtkMatrix4x4 :: DeepCopy( matrix, this->transform->GetMatrix() );

	for ( int i = 0; i < 16; i++ ){
		currentMatrix(i) = matrix[i];
	}

	this->candidateTransforms.push_back( currentMatrix );
	this->candidateCorrelationPeaks.push_back( this->correlationPeak );
	this->candidateCorrelationScales.push_back( this->correlationScale );
	//this->candidateWedgeOverlaps.push_back( this->getWedgeOverlap( t ) );
	this->candidateWedgeOverlaps.push_back( 1.0 );
}
