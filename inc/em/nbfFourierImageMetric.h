#pragma once

#include <vtkMath.h>
#include <vtkImageData.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfImageWriter.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfFourierFilter.h>

#include <nbfLinearInterpolator.h>
#include <nbfPolarDomain.h>

#include <fftw3.h>

#include <vnl/algo/vnl_powell.h>
#include <vnl/algo/vnl_amoeba.h>
#include <vnl/vnl_cost_function.h>

using namespace blitz;

/** Compute normalized cross-correlation between two images in Fourier space.
	Input images may be in real or already in reciprocal representation.

	Correlation peak is located with sub-pixel accuracy, using parabolic fit (Frank p. 82) in each axis.

	@todo Improve peak search computation, i.e. compute two most significant peaks and 
	use difference for validation purposes. 
*/
template< class Pixel, int const Dim >
class nbfFourierImageMetric : public nbfCorrelationImageMetric< Pixel, Dim >, public vnl_cost_function
{
public:

	/// Constructor with wedge object argument (neccesary to carry out FFT computations).
	nbfFourierImageMetric( nbfImageFilter< Pixel, Dim > * = NULL, nbfFourierFilter< Pixel, Dim > * = NULL );

	virtual ~nbfFourierImageMetric();

	void reinitialize( TinyVector< int, Dim > & );
	void reinitializeHalf( TinyVector< int, Dim > & );

	/// Redefine from parent
	void execute();
	//void execute(){	nbfCorrelationImageMetric< Pixel, Dim > :: execute(); }

	void execute( vtkTransform * );
	//void executeRobust( vtkTransform * );
	// Pixel getDistance(){ return this->getCorrelationPeak(); }
	void executeRotation( vtkTransform * );
	void executeAffine( vtkTransform * );
	void executeFourier( vtkTransform * );
	void executeFourierNew( vtkTransform * );
	void executeFourierNewHalf( vtkTransform * );
	void executeFourierNewHalf2D( vtkTransform * );
	void executeFourierNewWindows( vtkTransform * );
	void executeFourierNewWindowsReal( vtkTransform * );
	void executeFourierComplex( vtkTransform * );

	/// Redefine from vnl_cost_function
	double 	f ( vnl_vector< double > const & x );

	void setSeedTransform( vtkTransform * = NULL );

	void refine( vtkTransform *, Pixel = 1e-2 );

	void setInput1( nbfWedgedImage3D< Pixel > * i1 ){ nbfImageMetric< Pixel, Dim > :: setInput1(i1); this->input1ChangedLocally = true; }
	void setInput1( vtkImageData * i1 ){ nbfImageMetric< Pixel, Dim > :: setInput1(i1); this->input1ChangedLocally = true; }

	void setInput2( nbfWedgedImage3D< Pixel > * i2 ){ nbfImageMetric< Pixel, Dim > :: setInput2(i2); this->input2ChangedLocally = true; }
	void setInput2( vtkImageData * i2 ){ nbfImageMetric< Pixel, Dim > :: setInput2(i2); this->input2ChangedLocally = true; }

	// redefine
	Pixel getCorrelationPeakNoUpdate(){
		return this->correlationPeak; 
	}

	/// Distance value between inputs (function of the correlation value). Redefine from parent.
	virtual Pixel getDistance(){ Pixel p = nbfImageMetric< Pixel, Dim > :: getCorrelationPeak(); return this->correlationPeak; }

	// Compute feature vector as non-zero fourier components after aplication of wedge and fourier filter functions.
	void getLowDimensionRepresentation( Array< complex< Pixel >, 1 > &, Array< Pixel, 1 > & );
	void getLowDimensionRepresentationHalf( Array< complex< Pixel >, 1 > &, Array< Pixel, 1 > &, int = 0 );
	void getLowDimensionRepresentationReal( Array< Pixel, 1 > &, Pixel = 1.0, int = 0 );

	void putLowDimensionRepresentationReal( Array< Pixel, 1 > &, Pixel = 1.0 );

	// for CTF estimation
	void get2DimensionRepresentationHalf( Array< complex< Pixel >, 1 > &, Array< Pixel, 2 > &, int, bool = true );
	
protected:

	Pixel lastX, lastY, lastZ, lastDistance;
	Pixel oneBeforeLastX, oneBeforeLastY, oneBeforeLastZ, oneBeforeLastDistance;

	bool input1ChangedLocally, input2ChangedLocally;

	bool useWedge;

	vtkTransform * seedTransform;

	Array< complex< double >, Dim > FFT1re, FFT2re;
	// Array< complex< double >, Dim > FFTmaskreal;
	// Array< double, Dim > FFT1savedreal, FFT1sqrreal, FFT2real, maskDen;

	// Array< double, Dim > modSqrFFT1;
	Array< Pixel, Dim > wedge1, wedge2;

	Array< Pixel, Dim > mask, referenceMask, maskfft;
};

template< class Pixel, int const Dim  >
nbfFourierImageMetric< Pixel, Dim > :: nbfFourierImageMetric( nbfImageFilter< Pixel, Dim > * i, nbfFourierFilter< Pixel, Dim  > * f )
: nbfCorrelationImageMetric< Pixel, Dim >(i,f), vnl_cost_function(3)
{
	this->useWedge = true;
	this->seedTransform = vtkTransform :: New();
	this->seedTransform->Identity();
}


template< class Pixel, int const Dim  >
nbfFourierImageMetric< Pixel, Dim > :: ~nbfFourierImageMetric()
{
	this->seedTransform->Delete();
}

template< class Pixel, int const Dim >
double nbfFourierImageMetric< Pixel, Dim > :: f( vnl_vector< double > const & rotations )
{
	//if ( fabs(this->lastX-rotations[0]) + fabs(this->lastY-rotations[1]) + fabs(this->lastZ-rotations[2]) < .01 ){
	//	return this->lastDistance;
	//}

	//if ( fabs(this->oneBeforeLastX-rotations[0]) + fabs(this->oneBeforeLastY-rotations[1]) + fabs(this->oneBeforeLastZ-rotations[2]) < .01 ){
	//	return this->oneBeforeLastDistance;
	//}

	//// constrain search by making the cost function huge if outside desired range
	//if ( ( fabs( rotations[0] ) > this->restrictAngularSearch ) || ( fabs( rotations[1] ) > this->restrictAngularSearch ) || ( fabs( rotations[2] ) > this->restrictAngularSearch ) ){
	//	return numeric_limits< Pixel > :: max();
	//}

	double result;

	vtkTransform * transform = vtkTransform::New();
	transform->SetInput( this->seedTransform );

	// build trial transform on top of seed
	transform->RotateX( rotations[0] );
	transform->RotateY( rotations[1] );
	transform->RotateZ( rotations[2] );

	// check if trial transform is within allowed range
	if ( this->isTransformValid( transform ) == false ){
		result = numeric_limits< Pixel > :: max();
	} else {

		//this->execute( transform );
		this->executeFourierNewHalf( transform );

		if ( this->isTransformValid( transform ) == false ){
			cerr << "WARNING - Probing invalid transformation in 3D alignment routine. In " << __FILE__ ":" << __LINE__ << endl;
			cerr << "\t R=[" << rotations[0] << "," << rotations[1] << "," << rotations[2] << "], t=[" << transform->GetPosition()[0] << "," << transform->GetPosition()[1] << "," << transform->GetPosition()[2] << "]" << endl;
			result = numeric_limits< Pixel > :: max();
		} else {
			result = this->getCorrelationPeak();
		}
	
	}

	transform->Delete();

#if 0
	cout << "f (" << rotations[0] <<"," << rotations[1] << "," << rotations[2] << ") = " << result << endl;
	//cout << "f (" << rotations[0] <<"," << rotations[1] << "," << rotations[2] << ") = " << this->getDistance() << endl;
#endif

	//this->oneBeforeLastX = this->lastX;
	//this->oneBeforeLastY = this->lastY;
	//this->oneBeforeLastZ = this->lastZ;
	//this->oneBeforeLastDistance = this->lastDistance;

	//this->lastX = rotations[0];
	//this->lastY = rotations[1];
	//this->lastZ = rotations[2];
	//this->lastDistance = this->getCorrelationPeak();

	return result;
}


template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: setSeedTransform( vtkTransform * t )
{
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

	this->seedTransform->DeepCopy( rotation );

	//// this->seedTransform->DeepCopy( t );
	//vtkTransform * tmp = vtkTransform::New();
	//tmp->DeepCopy( t );
	//this->seedTransform->DeepCopy( tmp );
	//tmp->Delete();
	//this->seedTransform->SetInput( t );
}


template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: executeRotation( vtkTransform * myTransform )
{
	// get rotated image first
	this->wedgedInput2->getImage( this->input2, myTransform );
	this->setInput2( this->input2 );

	// compute corresponding translation
	nbfCorrelationImageMetric< Pixel, 3 > :: execute();
	
	vtkTransform * t = vtkTransform::New();
	//t->SetMatrix( myTransform->GetMatrix() );
	t->SetInput( myTransform );
	float translation[3];
	this->transform->GetPosition(translation);
	translation[0] *= -1;
	translation[1] *= -1;
	translation[2] *= -1;
	t->Translate(translation);
	
	// cout << *t << endl;

	// store translation in result
	vtkTransform * tmp = vtkTransform::New();
	this->transform->DeepCopy( tmp );
	tmp->Delete();
	//this->transform->SetMatrix( t->GetMatrix() );
	this->transform->SetInput( t );

	// now compute distance
	this->executeAffine( t );
	t->Delete();
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: executeAffine( vtkTransform * transform )
{
	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");

	// apply transform to second image
	if ( transform != NULL ){
		this->wedgedInput2->getImage( this->input2, transform );
		this->setInput2( this->input2 );
	}

	if ( this->input1Changed ){
		// build wedge filter
		this->wedge1.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
		this->wedgedInput1->getWedgeImage( this->wedge1 );
	}
	if ( this->input2Changed ){
		this->wedge2.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput2->getDimensions() );
		this->wedgedInput2->getWedgeImage( this->wedge2, transform );

		this->wedge2 = this->wedge1 * this->wedge2;

		this->fourierFilter->wedgeOn( this->wedge2 );
	}

	Pixel normalization = sum( this->wedge2 );
	
	Array< Pixel, Dim > A;

	if ( this->input1ChangedLocally == true ){
		
		// window data in real space
		this->imageFilter->execute( this->input1 );

		nbfVTKInterface::vtkToBlitz( this->input1, A );
		this->fourierFilter->initializeFFTW( TinyVector< int, 3 >( A.shape() ) );
		real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
		imag( this->fourierFilter->blitzFFT1 ) = 0;

		fftw_execute( this->fourierFilter->fftplan );

		this->modSqrFFT1.resize( this->fourierFilter->blitzFFT2.shape() );
		this->modSqrFFT1 = this->fourierFilter->blitzFFT2;

		this->input1ChangedLocally = false;
	} else {
		this->fourierFilter->blitzFFT2 = this->modSqrFFT1;
	}

	// apply filter
	this->fourierFilter->execute( this->fourierFilter->blitzFFT2 );
	// this->normalizeFourier( this->fourierFilter->blitzFFT2 );

	// Fourier modulo of first image
	this->modSqrFFT2.resize( this->fourierFilter->blitzFFT2.shape() );
	real(this->modSqrFFT2) = log( 1 + sqrt( real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ) ) );

	// window data2 in real space
	if ( this->input2Changed ){
		this->imageFilter->execute( this->input2 );
		this->input2Changed = false;
	}

	nbfVTKInterface::vtkToBlitz( this->input2, A );
	this->fourierFilter->initializeFFTW( TinyVector< int, 3 >( A.shape() ) );

	real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
	imag( this->fourierFilter->blitzFFT1 ) = 0;

	fftw_execute( this->fourierFilter->fftplan );

	// apply filter
	this->fourierFilter->execute( this->fourierFilter->blitzFFT2 );
	// this->normalizeFourier( this->fourierFilter->blitzFFT2 );

	// Fourier modulo of second image
	real(this->fourierFilter->blitzFFT2) = log( 1 + sqrt( real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ) ) );

	//this->fourierFilter->updateFilter();
	//this->wedge2 = this->wedge1 * this->wedge2 * this->fourierFilter->filter;

	// The L1 norm gives a much smoother energy function for minimization
	// Without normalization the peak is better localized.
	this->correlationPeak = sum( abs( real(this->modSqrFFT2) - real(this->fourierFilter->blitzFFT2 ) ) ) / normalization; 
	//this->correlationPeak = sum( abs( wedge1 - wedge2 ) / ( 1.0 * this->modSqrFFT1.size() ) ); // / normalization; 
	//this->correlationPeak = sum( abs( ( wedge1 - wedge2 ) / ( 1.0 + wedge2 ) ) ); // / normalization; 

	//this->correlationPeak = sum( abs( this->modSqrFFT1 / this->modSqrFFT2  - 1 ) * this->wedge2 ) / normalization;

	//this->modSqrFFT1 = this->modSqrFFT1 - min(this->modSqrFFT1); 
	//this->modSqrFFT1 = this->modSqrFFT1 / max(this->modSqrFFT1); 

	//this->modSqrFFT2 = this->modSqrFFT2 - min(this->modSqrFFT2); 
	//this->modSqrFFT2 = this->modSqrFFT2 / max(this->modSqrFFT2); 

	// The L2 norm gives a noisier energy landscape not amenable to local optimization.
	//this->correlationPeak = sum( pow2( log(1+sqrt(this->modSqrFFT1)) - log(1+sqrt(this->modSqrFFT2)) ) * this->wedge2 ); // / normalization; 
	//this->correlationPeak = sum( pow2( log(1+sqrt(this->modSqrFFT1)) - log(1+sqrt(this->modSqrFFT2)) ) * this->wedge1 * this->wedge2 ) / normalization; 

	//this->correlationPeak = sum( abs( sqrt(this->modSqrFFT1) - sqrt(this->modSqrFFT2) ) * this->wedge1 * this->wedge2 / normalization ); 
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: execute( vtkTransform * transform )
{
	if ( transform != NULL ){
		this->wedgedInput2->getImage( this->input2, transform );
		this->input2Changed = true;
	}

	if ( this->useMissingWedgeCompensation && this->input1Changed && this->input1ChangedLocally ){
		// build wedge filter
		this->wedge1.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
		this->wedgedInput1->getWedgeImage( this->wedge1 );
	}
	if ( this->useMissingWedgeCompensation &&  this->input2Changed ){
		this->wedge2.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput2->getDimensions() );
		this->wedgedInput2->getWedgeImage( this->wedge2, transform );

		this->wedge2 *= this->wedge1;

		// if overlap is almost entire volume ignore wedge
		this->fourierFilter->wedgeOn( this->wedge2 );
	}

	nbfCorrelationImageMetric< Pixel, 3 > :: execute();
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: executeFourier( vtkTransform * transform )
{
	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");

	if ( transform != NULL ){
		this->wedgedInput2->getImage( this->input2, transform );
		this->input2Changed = true;
	}

	if ( this->input1Changed && this->input1ChangedLocally ){
		// build wedge filter
		this->wedge1.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
		this->wedgedInput1->getWedgeImage( this->wedge1 );
	}
	if ( this->input2Changed ){
		this->wedge2.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput2->getDimensions() );
		this->wedgedInput2->getWedgeImage( this->wedge2, transform );
	
		this->wedge2 *= this->wedge1;

		// if overlap is almost entire volume ignore wedge
		this->fourierFilter->wedgeOn( this->wedge2 );
	}

	Array< double, Dim > A;

	if ( this->input1Changed == true ){
		// window data in real space
		this->imageFilter->execute( this->input1 );

		nbfVTKInterface::vtkToBlitzReference( this->input1, A );

		this->fourierFilter->initializeFFTW( A.shape() );
		real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
		imag( this->fourierFilter->blitzFFT1 ) = 0;

		fftw_execute( this->fourierFilter->fftplan );

		// store FFT1 to avoid recomputation
		this->FFT1saved.resize( this->fourierFilter->blitzFFT2.shape() );
		this->FFT1saved = this->fourierFilter->blitzFFT2;
		this->input1Changed = false;
	} else {
		this->fourierFilter->blitzFFT2 = this->FFT1saved;
	}

	// this->normalizeFourier( this->fourierFilter->blitzFFT2 );
	this->fourierFilter->blitzFFT2 = this->fourierFilter->blitzFFT2 * cast<double>(this->wedge2) * cast<double>(this->fourierFilter->filter);

	// norm squared of input 1
	this->modSqrFFT1.resize( this->wedge1.shape() );
	this->modSqrFFT1 = real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) );
	//if ( this->input1ChangedLocally == true ){
	//	// norm squared of input 1
	//	this->modSqrFFT1.resize( this->FFT1saved.shape() );
	//	this->modSqrFFT1 = log( 1 + sqrt( real( this->FFT1saved * conj( this->FFT1saved ) ) ) );
	//	//this->modSqrFFT1 = log( 1 + sqrt( real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ) ) );
	//	this->input1ChangedLocally = false;
	//}

	//w.write(this->modSqrFFT1);

	if ( this->input2Changed == true ){
	
		// window data in real space
		this->imageFilter->execute( this->input2 );

		nbfVTKInterface::vtkToBlitzReference( this->input2, A );
		this->fourierFilter->initializeFFTW( TinyVector< int, 3 >( A.shape() ) );

		real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
		imag( this->fourierFilter->blitzFFT1 ) = 0;

		fftw_execute( this->fourierFilter->fftplan );

		this->input2Changed = false;
	}

	//this->fourierFilter->updateFilter();
	
	// this->normalizeFourier( this->fourierFilter->blitzFFT2 );
	this->fourierFilter->blitzFFT2 = this->fourierFilter->blitzFFT2 * cast<double>(this->wedge2) * cast<double>(this->fourierFilter->filter);

	this->correlationPeak = sum( abs( log( 1 + sqrt(this->modSqrFFT1) ) - log( 1 + sqrt( real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ) ) ) ) ) / sum( this->wedge2 * this->fourierFilter->filter ); 

	// The L1 norm gives a much smoother energy function for minimization
	// Without normalization the peak is better localized.
	// this->correlationPeak = sum( abs( this->modSqrFFT1 - log( 1 + sqrt( real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ) ) ) ) * this->wedge2 * this->fourierFilter->filter ) / sum( this->wedge2 * this->fourierFilter->filter ); 
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: executeFourierNew( vtkTransform * transform )
{
	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	if ( transform != NULL ){
		this->wedgedInput2->getImage( this->input2, this->imageFilter, transform );
		this->input2Changed = true;
	}

	if ( this->input1Changed && this->input1ChangedLocally ){
		// build wedge filter
		this->wedge1.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
		this->wedgedInput1->getWedgeImage( this->wedge1 );
		//w.write(this->wedge1);
	}

	if ( this->input2Changed ){
		this->wedge2.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput2->getDimensions() );
		this->wedgedInput2->getWedgeImage( this->wedge2, transform );
		this->wedge2 *= this->wedge1;
		//w.write(this->wedge2);

		// if overlap is almost entire volume ignore wedge
		this->fourierFilter->wedgeOn( this->wedge2 );
	}

	Array< double, Dim > A;

	if ( this->input1Changed == true ){
		// window data in real space
		this->imageFilter->execute( this->input1 );

		nbfVTKInterface::vtkToBlitzReference( this->input1, A );

		//w.write(A);

		this->fourierFilter->initializeFFTW( A.shape() );
		real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
		imag( this->fourierFilter->blitzFFT1 ) = 0;

		fftw_execute( this->fourierFilter->fftplan );

		// store FFT1 to avoid recomputation
		this->FFT1saved.resize( this->fourierFilter->blitzFFT2.shape() );
		this->FFT1saved = this->fourierFilter->blitzFFT2;
		this->input1Changed = false;
	}
	//w.write(real(this->FFT1saved));

	this->fourierFilter->updateFilter();

	this->FFT1.resize( this->FFT1saved.shape() );
	this->FFT1 = this->FFT1saved;
	Pixel k1 = 1; // this->normalizeFourier( this->FFT1 );

	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");
	//this->fourierFilter->initializeFFTW( this->FFT1.shape() );	
	//this->fourierFilter->blitzFFT2 = this->FFT1;
	//fftw_execute( this->fourierFilter->ifftplan );
	//real( this->fourierFilter->blitzFFT1 ) = real( this->fourierFilter->blitzFFT1 ) * cast<double>(this->fourierFilter->shift);
	//w.write( real( this->FFT1 ) );

	if ( this->input2Changed == true ){
		// window data in real space
		this->imageFilter->execute( this->input2 );
		this->input2Changed = false;
	}

	// ALWAYS UPDATE INPUT 2
	nbfVTKInterface::vtkToBlitzReference( this->input2, A );
	//w.write(A);
	this->fourierFilter->initializeFFTW( TinyVector< int, 3 >( A.shape() ) );

	real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
	imag( this->fourierFilter->blitzFFT1 ) = 0;

	fftw_execute( this->fourierFilter->fftplan );

	//w.write(real(this->fourierFilter->blitzFFT2));

	Pixel k2 = 1;//this->normalizeFourier( this->fourierFilter->blitzFFT2 );

	//w.write(real(this->fourierFilter->blitzFFT2));

	this->correlationScale = k1 / k2;

	//fftw_execute( this->fourierFilter->ifftplan );
	//real( this->fourierFilter->blitzFFT1 ) = real( this->fourierFilter->blitzFFT1 ) * cast<double>(this->fourierFilter->shift);

	this->fourierFilter->blitzFFT2 = this->FFT1 * conj( this->fourierFilter->blitzFFT2 ) * cast<double>( this->fourierFilter->shift );

	// MCF - mutual correlation function
	if ( this->useMutualCorrelation == true ){
		real( this->fourierFilter->blitzFFT1 ) = pow( real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ), .25 );
		this->fourierFilter->blitzFFT2 = where( real( this->fourierFilter->blitzFFT1 ) > 0, this->fourierFilter->blitzFFT2 / real( this->fourierFilter->blitzFFT1 ), this->fourierFilter->blitzFFT2 );
	}

	// real( fftshift( ifft( imf1 .* conj(imf2) ) ) )
	fftw_execute( this->fourierFilter->ifftplan );

	Array< double, Dim > T( real( this->fourierFilter->blitzFFT1 ) );

	T *= this->fourierFilter->shift;

	// restrict rotation search to admisible range
	TinyVector< Pixel, Dim > center = T.shape() / 2.0;
	Range I( blitz :: extrema :: max( 0, center[0] - this->restrictTranslationSearch ), blitz :: extrema :: min( center[0] + this->restrictTranslationSearch, T.ubound(firstDim) ) );
	Array< double, Dim > vT( T(I,I,I) );

	if ( this->overlapNormalizedDistances == true ){
		this->correlationPeak = max( 0.0, 2 * ( 1 - max( vT ) / ( 1.0 * T.size() ) ) / sum( this->wedge2 * this->fourierFilter->filter ) * sum( this->fourierFilter->filter ) );
	} else {
		this->correlationPeak = max( 0.0, 2 * ( 1 - max( vT ) / ( 1.0 * T.size() ) ) );
	}

	// compute transformation shifts

	TinyVector< double, Dim > cccPeakPosition = this->locateMaxima( vT );
	cccPeakPosition = - ( floor(  vT.shape() / 2.0 ) - cccPeakPosition );

	// apply shifts to input transform

	//cout << -cccPeakPosition[0] << ", " << cccPeakPosition[1] << ", " << -cccPeakPosition[2] << endl;

	if ( Dim < 3 ){
		transform->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], 0 );
	} else{
		transform->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], - cccPeakPosition[thirdDim] );
	}
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: reinitialize( TinyVector< int, Dim > & s )
{
	this->imageFilter->buildMask( s );
	this->fourierFilter->initializeFFTW( s );

	// get mask and its FFT
	this->mask.resize( s );
	this->mask = this->imageFilter->mask * this->fourierFilter->shift;
	//real( this->fourierFilter->blitzFFT1 ) = cast< double >(this->mask);
	//imag( this->fourierFilter->blitzFFT1 ) = 0;
	//fftw_execute( this->fourierFilter->fftplan );
	//this->FFTmask.resize( this->fourierFilter->blitzFFT2.shape() );
	//this->FFTmask = this->fourierFilter->blitzFFT2;

	//this->fourierFilter->blitzFFTreal = cast< double >(this->mask);
	//fftw_execute( this->fourierFilter->fftplanreal );

	//TinyVector< int, 3 > shape( s[0], s[1], s[2]/2+1 );
	//Array< complex< double >, 3 > B( reinterpret_cast<complex<double>*>(this->fourierFilter->blitzFFT.data()), shape, neverDeleteData );

	//this->FFTmaskreal.resize( B.shape() );
	//this->FFTmaskreal = B;

}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: reinitializeHalf( TinyVector< int, Dim > & s )
{
	this->imageFilter->buildMask( s );
	this->fourierFilter->initializeFFTWhalf( s );

	// get mask and its FFT
	this->mask.resize( s );
	this->mask = this->imageFilter->mask;

	this->maskfft.resize( this->mask.shape() );
	int size = this->imageFilter->squareMaskSize;
	if ( ( size > 0 ) && ( size < floor( this->maskfft.shape()[0] / 2.0 ) ) ){
		this->maskfft = 0;
		if ( Dim == 2 ){
			if ( ( size < this->maskfft.rows() ) && ( size < this->maskfft.cols() ) ){
				this->maskfft( Range(size,this->maskfft.ubound(firstDim)-size), Range(size,this->maskfft.ubound(secondDim)-size) ) = 1;
			} else {
				cout << "ERROR - mask edge size does not fit current image size. In " << __FILE__ << "," << __LINE__ << endl;
			}
		} else {
			if ( ( size < this->maskfft.rows() ) && ( size < this->maskfft.cols() ) && ( size < this->maskfft.depth() ) ){
				this->maskfft( Range(size,this->maskfft.ubound(firstDim)-size), Range(size,this->maskfft.ubound(secondDim)-size), Range(size,this->maskfft.ubound(thirdDim)-size) ) = 1;
			} else {
				cout << "ERROR - mask edge size does not fit current image size. In " << __FILE__ << "," << __LINE__ << endl;
			}
		}
	} else {
		cerr << "WARNING - ignoring square image mask due to invalid mask dimensions." << endl;
		this->maskfft = 1;
	}
	this->imageFilter->softenMask( this->maskfft, 2.5 );
	// this->maskfft *= this->fourierFilter->shift;

	this->imageFilter->fileMaskOn();
	this->imageFilter->buildMaskFromFile(s);
	if ( this->imageFilter->mask.size() > 0 ){
		this->referenceMask.resize( this->imageFilter->mask.shape() );
		this->referenceMask = this->imageFilter->mask;
	} else {
		this->referenceMask.resize( this->mask.shape() );
		this->referenceMask = this->mask;
	}
	this->imageFilter->fileMaskOff();
	this->imageFilter->buildMask( s );

	//this->fourierFilter->blitzFFTreal = cast< double >(this->mask);
	//fftw_execute( this->fourierFilter->fftplanreal );

	//TinyVector< int, 3 > shape( s[0], s[1], s[2]/2+1 );
	//Array< complex< double >, 3 > B( reinterpret_cast<complex<double>*>(this->fourierFilter->blitzFFT.data()), shape, neverDeleteData );
	//Array< Pixel, 3 > shiftr( this->fourierFilter->shift( Range :: all(), Range :: all(), Range(fromStart,shape[2]-1) ) );

	//this->FFTmaskreal.resize( B.shape() );
	//this->FFTmaskreal = B;

	//B = this->FFTmaskreal * conj( this->FFTmaskreal ) * cast<double>(shiftr);
	//fftw_execute( this->fourierFilter->ifftplanreal );

	//this->maskDen.resize( this->fourierFilter->blitzFFTreal.shape() );
	//this->maskDen = this->fourierFilter->blitzFFTreal / this->fourierFilter->blitzFFTreal.size();
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: executeFourierNewWindows( vtkTransform * transform )
{
	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	if ( transform != NULL ){
		this->wedgedInput2->getImage( this->input2, this->imageFilter, transform );
		this->input2Changed = true;
	}

	if ( this->input1Changed && this->input1ChangedLocally ){
		// build wedge filter
		this->wedge1.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
		this->wedgedInput1->getWedgeImage( this->wedge1 );
	}

	if ( this->input2Changed ){
		this->wedge2.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput2->getDimensions() );
		this->wedgedInput2->getWedgeImage( this->wedge2, transform );
		this->wedge2 *= this->wedge1;

		// if overlap is almost entire volume ignore wedge
		this->fourierFilter->wedgeOn( this->wedge2 );
	}

	// initialization
	if ( this->mask.size() == 0 ){
		TinyVector< int, Dim > s = this->wedge2.shape();
		this->reinitialize(s);
	}
	this->fourierFilter->updateFilter();
	imag( this->fourierFilter->blitzFFT1 ) = 0;

	Array< double, Dim > A;

	if ( this->input1Changed == true ){
		// apply image filter by-passing the mask
		this->imageFilter->execute( this->input1, true );
		nbfVTKInterface::vtkToBlitzReference( this->input1, A );
		if ( this->useMutualCorrelation == true ){
			this->fourierFilter->execute( A, 1 );
		} else {
			this->fourierFilter->execute( A );
		}
		real( this->fourierFilter->blitzFFT1 ) = A * this->mask;
		fftw_execute( this->fourierFilter->fftplan );
		this->FFT1saved.resize( this->fourierFilter->blitzFFT2.shape() );
		this->FFT1saved = this->fourierFilter->blitzFFT2;
		
		real( this->fourierFilter->blitzFFT1 ) *= A;
		fftw_execute( this->fourierFilter->fftplan );
		this->FFT1sqr.resize( this->fourierFilter->blitzFFT2.shape() );
		this->FFT1sqr = this->fourierFilter->blitzFFT2;

		this->input1Changed = false;
	}

	if ( this->input2Changed == true ){
		// window data in real space
		this->imageFilter->execute( this->input2, true );
		nbfVTKInterface::vtkToBlitzReference( this->input2, A );
		if ( this->useMutualCorrelation == true ){
			this->fourierFilter->execute( A, 1 );
		} else {
			this->fourierFilter->execute( A );
		}
		this->input2Changed = false;
	}

	// ALWAYS UPDATE INPUT 2
	nbfVTKInterface::vtkToBlitzReference( this->input2, A );
	real( this->fourierFilter->blitzFFT1 ) = A * this->mask;
	fftw_execute( this->fourierFilter->fftplan );
	this->FFT2.resize( this->fourierFilter->blitzFFT2.shape() );
	this->FFT2 = this->fourierFilter->blitzFFT2;

	real( this->fourierFilter->blitzFFT1 ) *= A;
	fftw_execute( this->fourierFilter->fftplan );

	// windowed cross correlation
	this->fourierFilter->blitzFFT2 = ( this->FFT1sqr * conj( this->FFTmask )
		- 2.0 * this->FFT1saved * conj( this->FFT2 )
		+ this->FFTmask * conj( this->fourierFilter->blitzFFT2 ) )
		* cast<double>( this->fourierFilter->shift * this->wedge2 );
	//this->fourierFilter->blitzFFT2 = - this->FFT1saved * conj( this->FFT2 ) * cast<double>( this->fourierFilter->shift );
	//this->fourierFilter->blitzFFT2 = where( real( this->fourierFilter->blitzFFT2 * conj(this->fourierFilter->blitzFFT2) ) > 0, this->fourierFilter->blitzFFT2 / sqrt( real( this->fourierFilter->blitzFFT2 * conj(this->fourierFilter->blitzFFT2) ) ), this->fourierFilter->blitzFFT2 );
	//this->fourierFilter->blitzFFT2 = - this->FFT1saved * conj( this->FFT2 ) * cast<double>( this->fourierFilter->shift * this->wedge2 );
	//this->fourierFilter->execute( this->fourierFilter->blitzFFT2 );

	// real( fftshift( ifft( imf1 .* conj(imf2) ) ) )
	fftw_execute( this->fourierFilter->ifftplan );

	// restrict rotation search to admisible range
	TinyVector< Pixel, Dim > center = this->fourierFilter->blitzFFT2.shape() / 2.0;
	Range I( blitz :: extrema :: max( 0, center[0] - this->restrictTranslationSearch ), blitz :: extrema :: min( center[0] + this->restrictTranslationSearch, this->fourierFilter->blitzFFT2.ubound(firstDim) ) );

	Array< double, Dim > T( real( this->fourierFilter->blitzFFT1(I,I,I) ) );
	T *= this->fourierFilter->shift(I,I,I);

	this->correlationScale = 1.0;

	if ( this->overlapNormalizedDistances == true ){
		this->correlationPeak = max( 0.0, min( T ) / pow2( 1.0 * T.size() ) ) / sum( this->wedge2 * this->fourierFilter->filter ) * sum( this->fourierFilter->filter );
	} else {
		this->correlationPeak = max( 0.0, min( T ) / pow2( 1.0 * T.size() ) );
	}

	// compute transformation shifts
	T = - T;
	TinyVector< double, Dim > cccPeakPosition = this->locateMaxima( T );
	cccPeakPosition = - ( floor(  T.shape() / 2.0 ) - cccPeakPosition );

	// apply shifts to input transform

	//cout << "d=" << this->correlationPeak << ", t=[" << -cccPeakPosition[0] << ", " << cccPeakPosition[1] << ", " << -cccPeakPosition[2] << "]"<< endl;

	if ( Dim < 3 ){
		transform->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], 0 );
	} else{
		transform->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], - cccPeakPosition[thirdDim] );
	}
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: executeFourierNewWindowsReal( vtkTransform * transform )
{
	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	if ( transform != NULL ){
		this->wedgedInput2->getImage( this->input2, this->imageFilter, transform );
		this->input2Changed = true;
	}

	if ( this->input1Changed && this->input1ChangedLocally ){
		// build wedge filter
		this->wedge1.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
		this->wedgedInput1->getWedgeImageHalf( this->wedge1 );
	}

	if ( this->input2Changed ){
		this->wedge2.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput2->getDimensions() );
		this->wedgedInput2->getWedgeImageHalf( this->wedge2, transform );
		this->wedge2 *= this->wedge1;

		// if overlap is almost entire volume ignore wedge
		this->fourierFilter->wedgeOn( this->wedge2 );
	}

	// initialization
	if ( this->mask.size() == 0 ){
		TinyVector< int, Dim > s = this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions();
		this->reinitializeHalf(s);
	}
	this->fourierFilter->updateFilterHalf();

	Array< double, Dim > A;

	if ( this->input1Changed == true ){
		// apply image filter by-passing the mask
		this->imageFilter->execute( this->input1, true );
		this->input1Changed = false;
	}
	
	nbfVTKInterface :: vtkToBlitz( this->input1, A );
	if ( this->useMutualCorrelation == true ){
		this->fourierFilter->executeHalf( A, 1 );
	} else {
		this->fourierFilter->executeHalf( A );
	}
	
	this->fourierFilter->blitzFFTreal = A * this->mask;
	fftw_execute( this->fourierFilter->fftplanreal );
	this->FFT1savedreal.resize( this->fourierFilter->blitzFFT.shape() );
	this->FFT1savedreal = this->fourierFilter->blitzFFT;

	this->fourierFilter->blitzFFTreal = A * A * this->mask;
	fftw_execute( this->fourierFilter->fftplanreal );
	this->FFT1sqrreal.resize( this->fourierFilter->blitzFFT.shape() );
	this->FFT1sqrreal = this->fourierFilter->blitzFFT;

	if ( this->input2Changed == true ){
		// window data in real space
		this->imageFilter->execute( this->input2, true );
		nbfVTKInterface::vtkToBlitzReference( this->input2, A );
		if ( this->useMutualCorrelation == true ){
			this->fourierFilter->executeHalf( A, 1 );
		} else {
			this->fourierFilter->executeHalf( A );
		}
		this->input2Changed = false;
	}

	// ALWAYS UPDATE INPUT 2
	nbfVTKInterface :: vtkToBlitzReference( this->input2, A );
	this->fourierFilter->blitzFFTreal = A * this->mask;
	fftw_execute( this->fourierFilter->fftplanreal );
	this->FFT2real.resize( this->fourierFilter->blitzFFT.shape() );
	this->FFT2real = this->fourierFilter->blitzFFT;

	TinyVector< int, 3 > shape( A.rows(), A.cols(), A.depth() / 2 + 1 );
	Array< Pixel, 3 > shiftr( this->fourierFilter->shift( Range :: all(), Range :: all(), Range(fromStart,shape[2]-1) ) );
	Array< complex< double >, 3 > FFT2sqrre( reinterpret_cast<complex<double>*>(this->fourierFilter->blitzFFT.data()), shape, neverDeleteData );
	Array< complex< double >, 3 > FFT1re( reinterpret_cast<complex<double>*>(this->FFT1savedreal.data()), shape, neverDeleteData );
	Array< complex< double >, 3 > FFT1sqrre( reinterpret_cast<complex<double>*>(this->FFT1sqrreal.data()), shape, neverDeleteData );
	Array< complex< double >, 3 > FFT2re( reinterpret_cast<complex<double>*>(this->FFT2real.data()), shape, neverDeleteData );

	this->fourierFilter->blitzFFTreal = A * A * this->mask;
	fftw_execute( this->fourierFilter->fftplanreal );

	// normalized windowed cross correlation
	FFT2sqrre = ( FFT1sqrre * conj( this->FFTmaskreal )
				- 2.0 * FFT1re * conj( FFT2re )
				+ this->FFTmaskreal * conj( FFT2sqrre ) )
				* cast<double>( shiftr );
	//FFT2sqrre = - 2.0 * FFT1re * conj( FFT2re ) * cast<double>( shiftr );
	
	// real( fftshift( ifft( imf1 .* conj(imf2) ) ) )
	fftw_execute( this->fourierFilter->ifftplanreal );

	// restrict rotation search to admisible range
	TinyVector< Pixel, Dim > center = this->fourierFilter->blitzFFTreal.shape() / 2.0;
	Range I( blitz :: extrema :: max( 0, center[0] - this->restrictTranslationSearch ), blitz :: extrema :: min( center[0] + this->restrictTranslationSearch, A.ubound(firstDim) ) );

	// normalized windowed L2 distance
	Array< double, Dim > T( this->fourierFilter->blitzFFTreal(I,I,I) * this->fourierFilter->shift(I,I,I) );
	//Array< double, Dim > T( this->fourierFilter->blitzFFTreal(I,I,I) / this->maskDen(I,I,I) );

	if ( this->overlapNormalizedDistances == true ){
		//this->correlationPeak = max( 0.0, min( T ) / pow2( 1.0 * T.size() ) );
		this->correlationPeak = max( 0.0, min( T ) / pow2( 1.0 * T.size() ) ) / sum( this->wedge2 * this->fourierFilter->filter ) * sum( this->fourierFilter->filter );
	} else {
		this->correlationPeak = max( 0.0, min( T ) / pow2( 1.0 * T.size() ) );
	}
	this->correlationScale = 1.0;

	// compute transformation shifts
	T = - T;
	TinyVector< double, Dim > cccPeakPosition = this->locateMaxima( T );
	cccPeakPosition = - ( floor(  T.shape() / 2.0 ) - cccPeakPosition );

	// apply shifts to input transform

	//cout << "d=" << this->correlationPeak << ", t=[" << -cccPeakPosition[0] << ", " << cccPeakPosition[1] << ", " << -cccPeakPosition[2] << "]"<< endl;

	if ( Dim < 3 ){
		transform->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], 0 );
	} else{
		transform->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], - cccPeakPosition[thirdDim] );
	}
}

//template< class Pixel, int const Dim >
//void nbfFourierImageMetric< Pixel, Dim > :: executeFourierNewHalf( vtkTransform * transform )
//{
//	nbfMatlabWriter w;
//	w.setFileName("p.matlab");
//
//	// apply transformation to second image
//	this->wedgedInput2->getImage( this->input2, transform );
//
//	// get first missing wedge
//	if ( this->wedgedInput1->isWedgeEffective() && this->input1Changed ){
//		this->wedge1.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
//		this->wedgedInput1->getWedgeImageHalf( this->wedge1 );
//	}
//	// get second missing wedge
//	if ( this->wedgedInput2->isWedgeEffective() ){
//		this->wedge2.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput2->getDimensions() );
//		this->wedgedInput2->getWedgeImageHalf( this->wedge2, transform );
//	}
//	// get effective missing wedge
//	if ( this->wedgedInput1->isWedgeEffective() && this->wedgedInput2->isWedgeEffective() ){
//		this->wedge2 *= this->wedge1;
//		this->fourierFilter->wedgeOn( this->wedge2 );
//	} else if ( this->wedgedInput1->isWedgeEffective() ) {
//		this->fourierFilter->wedgeOn( this->wedge1 );
//	} else if ( this->wedgedInput2->isWedgeEffective() ) {
//		this->fourierFilter->wedgeOn( this->wedge2 );
//	} else {
//		this->fourierFilter->wedgeOff();
//	}
//
//	this->input2Changed = false;
//	this->input1Changed = false;
//
//	// initialize image mask
//	if ( this->mask.size() == 0 ){
//		TinyVector< int, Dim > s = this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions();
//		this->reinitializeHalf(s);
//	}
//
//	Array< double, Dim > A;
//
//	vtkImageData * cdata = vtkImageData :: New();
//	cdata->DeepCopy( this->input1 );
//
//	// apply filters in reciprocal space (wedge + bandpass + mutual)
//	nbfVTKInterface :: vtkToBlitzReference( cdata, A );
//	if ( this->useMutualCorrelation ){
//		this->fourierFilter->executeHalf(A,1);
//	} else {
//		this->fourierFilter->executeHalf(A,0);
//	}
//
//	// apply image filter by-passing the mask
//	this->imageFilter->execute( cdata, true );
//
//	// compute FFT of first image
//	this->fourierFilter->blitzFFTreal = A * this->mask;
//	fftw_execute( this->fourierFilter->fftplanreal );
//	this->FFT1savedreal.resize( this->fourierFilter->blitzFFT.shape() );
//	this->FFT1savedreal = this->fourierFilter->blitzFFT;
//
//	// create handles to half-size FFTs
//	TinyVector< int, 3 > shape( A.rows(), A.cols(), A.depth() / 2 + 1 );
//	Array< Pixel, 3 > shiftr( this->fourierFilter->shift( Range :: all(), Range :: all(), Range(fromStart,shape[2]-1) ) );
//	Array< complex< double >, 3 > FFT1re( reinterpret_cast<complex<double>*>(this->FFT1savedreal.data()), shape, neverDeleteData );
//	Array< complex< double >, 3 > FFT2re( reinterpret_cast<complex<double>*>(this->fourierFilter->blitzFFT.data()), shape, neverDeleteData );
//
//	// apply Fourier filter and normalize first image
//	Pixel k1 = this->normalizeFourierHalf( FFT1re );
//	//Pixel k1 = this->normalizeFourierHalf( FFT1re, this->useMutualCorrelation );
//	
//	// compute FFT of second image
//	cdata->DeepCopy( this->input2 );
//
//	nbfVTKInterface :: vtkToBlitzReference( cdata, A );
//	if ( this->useMutualCorrelation ){
//		this->fourierFilter->executeHalf(A,1);
//	} else {
//		this->fourierFilter->executeHalf(A,0);
//	}
//
//	// apply image filter by-passing the mask
//	this->imageFilter->execute( cdata, true );
//
//	this->fourierFilter->blitzFFTreal = A * this->mask;
//	fftw_execute( this->fourierFilter->fftplanreal );
//
//	// apply Fourier filter and normalize
//	Pixel k2 = this->normalizeFourierHalf( FFT2re );
//	//Pixel k2 = this->normalizeFourierHalf( FFT2re, this->useMutualCorrelation );
//
//	// normalized cross correlation
//	FFT2re = FFT1re * conj( FFT2re ) * cast<double>( shiftr );
//	fftw_execute( this->fourierFilter->ifftplanreal );
//
//	// restrict search to admisible range
//	TinyVector< Pixel, Dim > center = this->fourierFilter->blitzFFTreal.shape() / 2.0;
//	Range I( blitz :: extrema :: max( 0, center[0] - this->restrictTranslationSearch ), blitz :: extrema :: min( center[0] + this->restrictTranslationSearch, A.ubound(firstDim) ) );
//	Array< double, Dim > T( 1.0 + this->fourierFilter->blitzFFTreal(I,I,I) * this->fourierFilter->shift(I,I,I) / A.size() );
//	if ( ( this->overlapNormalizedDistances == true ) 
//		&& ( this->wedgedInput1->isWedgeEffective() || this->wedgedInput2->isWedgeEffective() ) ){
//		//this->correlationPeak = max( 0.0, min( T ) / pow2( 1.0 * T.size() ) );
//		this->correlationPeak = max( 0.0, 2.0 * ( 2.0 - max( T ) ) ) / sum( this->fourierFilter->wedge * this->fourierFilter->filter ) * sum( this->fourierFilter->filter );
//	} else {
//		this->correlationPeak = max( 0.0, 2.0 * ( 2.0 - max( T ) ) );
//	}
//	this->correlationScale = 1.0;
//
//	// compute transformation shifts
//	TinyVector< double, Dim > cccPeakPosition = this->locateMaxima( T );
//	cccPeakPosition = - ( floor(  T.shape() / 2.0 ) - cccPeakPosition );
//
//	cdata->Delete();
//
//#ifdef WIN32
//	// cout << "d=" << this->correlationPeak << ", t=[" << -cccPeakPosition[0] << ", " << cccPeakPosition[1] << ", " << -cccPeakPosition[2] << "]"<< endl;
//#endif
//
//	// apply shifts to input transform
//
//	if ( Dim < 3 ){
//		transform->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], 0 );
//	} else{
//		transform->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], - cccPeakPosition[thirdDim] );
//	}
//}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: executeFourierNewHalf( vtkTransform * transform )
{

	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	// initialize image mask
	TinyVector< int, Dim > s = this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions();
	if ( this->mask.size() == 0 ){
        // cout << "Initializing image mask with size = " << s << endl;
        this->reinitializeHalf(s);
	}

	// create handles to half-size FFTs
	TinyVector< int, Dim > shape( this->mask.rows(), this->mask.cols(), this->mask.depth() / 2 + 1 );
	Array< Pixel, Dim > shiftr( this->fourierFilter->shift( Range :: all(), Range :: all(), Range(fromStart,shape[2]-1) ) );
	Array< complex< double >, Dim > FFTre( reinterpret_cast<complex<double>*>(this->fourierFilter->blitzFFT.data()), shape, neverDeleteData );

    // create handles to half-size big FFTs
	TinyVector< int, Dim > shapebig( 2 * this->mask.rows(), 2 * this->mask.cols(), 2 * this->mask.depth() / 2 + 1 );
	Array< complex< double >, Dim > FFTrebig( reinterpret_cast<complex<double>*>(this->fourierFilter->blitzFFTbig.data()), shapebig, neverDeleteData );

	Array< double, Dim > A;
	if ( this->useMissingWedgeCompensation == true ){
		if ( this->input1Changed == true ){
			this->wedge1.resize( s );
			this->wedgedInput1->getWedgeImageHalf( this->wedge1 );
		}

		this->wedge2.resize( s );
		this->wedgedInput2->getWedgeImageHalf( this->wedge2, transform );
		//w.write(this->wedge2);
	
		this->wedge2 *= this->wedge1;
		this->fourierFilter->wedgeOn( this->wedge2 );
	} else {
        // cout << "Not using wedge compensation" << endl;
		this->wedge1.free();
		this->wedge2.free();
		this->fourierFilter->wedgeOff();
	}

	// compute \int ( I_1 - I_2 )^2 \times m_1 in reciprocal space

	// compute FFT of I_1 alone
	nbfVTKInterface :: vtkToBlitzReference( this->input1, A );

    this->fourierFilter->blitzFFTreal = A * this->maskfft * this->fourierFilter->shift;
	fftw_execute( this->fourierFilter->fftplanreal );
	this->input1Changed = false;

	//Array< double, 3 > H( FFTre.shape() );
	//H = real( FFTre );
	//w.write(H);

	// apply band pass + wedge filter + mutual correlation scaling
	this->fourierFilter->executeHalf( FFTre, this->useMutualCorrelation );
	
	// apply normalization in reciprocal space
	Pixel k1 = this->normalizeFourierHalf( FFTre, false );
	
    // invert FFT
	fftw_execute( this->fourierFilter->ifftplanreal );

    // this->fourierFilter->blitzFFTreal *= this->fourierFilter->shift;

	Pixel k1_m = sum( this->fourierFilter->blitzFFTreal * this->fourierFilter->blitzFFTreal * this->referenceMask );

	// apply mask in real space
	this->fourierFilter->blitzFFTreal *= this->referenceMask;
	// this->fourierFilter->blitzFFTreal = this->fourierFilter->blitzFFTreal * this->referenceMask * this->fourierFilter->shift;
	
	// compute final FFT of I_1 \times m_1
	fftw_execute( this->fourierFilter->fftplanreal );
	this->FFT1re.resize( FFTre.shape() );
	this->FFT1re = FFTre;

	// compute FFT of I_2 alone
	this->wedgedInput2->getImage( this->input2, transform );
	this->input2Changed = false;

	// NOT TESTED - POTENTIAL ISSUES WITH SYMMETRY AXIS
	// // Application of symmetry on raw volume before evaluation of CCC
	// if ( this->imageFilter->symmetryFactor > 1 ){
		// this->imageFilter->symmetrize( this->input2, this->imageFilter->symmetryFactor );
	// }	
	
	nbfVTKInterface :: vtkToBlitzReference( this->input2, A );
    // w.write(A);
	this->fourierFilter->blitzFFTreal = A * this->maskfft * this->fourierFilter->shift;
	fftw_execute( this->fourierFilter->fftplanreal );

	//H = real( FFTre );
	//w.write(H);

	// apply band pass + wedge filter + mutual correlation scaling
	this->fourierFilter->executeHalf( FFTre, this->useMutualCorrelation );

	// apply normalization in reciprocal space
	Pixel k2 = this->normalizeFourierHalf( FFTre, false );

	// store normalized FFT of I_2 (fix fftw scaling)
	this->FFT2re.resize( FFTre.shape() );
	double normf = this->fourierFilter->shift.size();
	this->FFT2re = 2.0 * this->FFT1re * conj( FFTre * normf );

	// invert FFT
	fftw_execute( this->fourierFilter->ifftplanreal );

	//H = this->fourierFilter->blitzFFTreal * this->fourierFilter->shift;
	//w.write(H);
    // this->fourierFilter->blitzFFTreal *= this->fourierFilter->shift;

	// FFT of cuadratic term I_2^2
	this->fourierFilter->blitzFFTreal = this->fourierFilter->blitzFFTreal * this->fourierFilter->blitzFFTreal * this->fourierFilter->shift;
	fftw_execute( this->fourierFilter->fftplanreal );
	this->FFT1re.resize( FFTre.shape() );
	this->FFT1re = FFTre;

	// compute FFT of m_1
	this->fourierFilter->blitzFFTreal = this->referenceMask * this->fourierFilter->shift;
	fftw_execute( this->fourierFilter->fftplanreal );
	
	// evaluate quadratic expression in reciprocal space
	FFTre = ( - FFTre * conj( this->FFT1re ) + this->FFT2re ) * cast<double>( shiftr );

	double factor = ( 1.0 * FFTrebig.rows() * FFTrebig.cols() * 2 * ( FFTrebig.depth() - 1 ) ) / ( 1.0 * FFTre.rows() * FFTre.cols() * 2 * ( FFTre.depth() - 1 ) );

	// normalized cross correlation
	//FFTre = this->FFT1re * conj( FFTre ) * cast<double>( shiftr );
	
	Pixel upsampling_factor = 2.0;
	
	// put into up-sampled fft
	FFTrebig = 0;
	FFTrebig( Range( A.rows() - floor(A.rows()/2.0), A.rows() + floor((A.rows()-1)/2.0) ), Range(A.cols()-floor(A.cols()/2.0), A.cols()+floor((A.cols()-1)/2.0)), Range( FFTre.depth() - 1, FFTrebig.ubound(thirdDim) ) ) = FFTre * factor;
	fftw_execute( this->fourierFilter->ifftplanrealbig );
	
	// locate maxima in upsampled correlation volume

	//Array< double, Dim > Tf( 1.0 + this->fourierFilter->blitzFFTrealbig * this->fourierFilter->shiftbig / this->fourierFilter->shiftbig.size() );
	//w.write(Tf);

	// restrict search to admisible range
	TinyVector< Pixel, Dim > center = this->fourierFilter->blitzFFTrealbig.shape() / 2.0;
	Pixel restriction = center[0];
	if ( this->restrictTranslationSearch < numeric_limits< Pixel > :: max() ){
		restriction = blitz :: extrema :: max( this->restrictTranslationSearch * upsampling_factor, 1.0 );
	}
	Range I( blitz :: extrema :: max( 0, center[0] - floor(restriction) ), blitz :: extrema :: min( center[0] + floor(restriction), this->fourierFilter->blitzFFTrealbig.ubound(firstDim) ) );
	Range J( blitz :: extrema :: max( 0, center[1] - floor(restriction) ), blitz :: extrema :: min( center[1] + floor(restriction), this->fourierFilter->blitzFFTrealbig.ubound(secondDim) ) );
	Range K( blitz :: extrema :: max( 0, center[2] - floor(restriction) ), blitz :: extrema :: min( center[2] + floor(restriction), this->fourierFilter->blitzFFTrealbig.ubound(thirdDim) ) );
	////Array< double, Dim > T( 1.0 + this->fourierFilter->blitzFFTrealbig(I,I,I) * this->fourierFilter->shiftbig(I,I,I) / this->fourierFilter->shiftbig.size() );
	Array< double, Dim > T( this->fourierFilter->blitzFFTrealbig(I,J,K) * this->fourierFilter->shiftbig(I,J,K) / this->fourierFilter->shiftbig.size() );
	
	if ( ( A.rows() != A.cols() ) || ( A.cols() != A.depth() ) ){
		cerr << "WARNING: Trying to align volumes with different sizes in each dimension. THIS HAS NOT BEEN TESTED PROPERLY." << endl;
		if ( k1_m < max( T ) ){
			T = -T;
		}
	}
    // w.write(T);

	//Array< double, Dim > T( ( - k1 - k2 + 2.0 * this->fourierFilter->blitzFFTreal(I,I,I) * this->fourierFilter->shift(I,I,I) ) / A.size() );
	Array< Pixel, Dim > radius( T.shape() );
	firstIndex i; secondIndex j; thirdIndex k;
	if ( Dim == 2 ) radius = ( restriction - i ) * ( restriction - i ) + ( restriction - j ) * ( restriction - j );
	//if ( Dim == 3 ) radius = ( restriction - i ) * ( restriction - i ) + ( restriction - j ) * ( restriction - j ) + ( restriction - k ) * ( restriction - k );
	if ( Dim == 3 ) radius = ( T.ubound(0) / 2 - i ) * ( T.ubound(0) / 2 - i ) + ( T.ubound(1) / 2 - j ) * ( T.ubound(1) / 2 - j ) + ( T.ubound(2) / 2 - k ) * ( T.ubound(2) / 2 - k );
	T = where( radius < restriction * restriction, T, - numeric_limits< double > :: max() );
	//w.write(T);
	if ( ( this->overlapNormalizedDistances == true ) && ( this->fourierFilter->wedge.size() > 0 ) ){
		this->correlationPeak = max( 0.0, ( k1_m - max( T ) ) ) / sum( this->fourierFilter->wedge * this->fourierFilter->filter ) * sum( this->fourierFilter->filter );
		//this->correlationPeak = max( 0.0, 2.0 * ( 2.0 - max( T ) ) ) / sum( this->fourierFilter->wedge * this->fourierFilter->filter ) * sum( this->fourierFilter->filter );
		//this->correlationPeak = max( 0.0, - max( T ) ) / sum( this->fourierFilter->wedge * this->fourierFilter->filter ) * sum( this->fourierFilter->filter );
	} else {
		this->correlationPeak = max( 0.0, ( k1_m - max( T ) ) );
		//this->correlationPeak = m- max( T );
		//this->correlationPeak = max( 0.0, 2.0 * ( 2.0 - max( T ) ) );
		//this->correlationPeak = max( 0.0, - max( T ) );
	}
	this->correlationPeak = this->correlationPeak / this->fourierFilter->shiftbig.size() / this->fourierFilter->shiftbig.size();
	this->correlationScale = 1.0;

	TinyVector< double, Dim > cccPeakPosition = this->locateMaxima( T );
	cccPeakPosition = cccPeakPosition - ( T.shape() - 1 ) / 2.0;
	cccPeakPosition /= upsampling_factor;

	//// sub-pixel
	//TinyVector< int, 3 > shift;
	//for ( int i = 0; i < Dim; i++ ){
	//	if ( cccPeakPosition[0] > floor( A.rows() / 2.0 ) ){
	//		shift[i] = cccPeakPosition[i] - A.extent(i) - 1;
	//	} else {
	//		shift[i] = cccPeakPosition[i] - 1;
	//	}
	//}
	//shift /= 2.0;

	//int upsampling = 10;
	//// shift = round( shift * upsampling ) / upsampling; 
	//int dftshift = floor( ceil( upsampling * 1.5 ) / 2.0 ); // Center of output array at dftshift+1

	//fourthIndex l;

	//Array< complex< double >, 2 > E( FFTre.depth(), upsampling );
	//real(E) = cos( - 2 * vtkMath::Pi() * ( i - FFTre.rows() / 2 ) * j / FFTre.rows() );
	//imag(E) = sin( - 2 * vtkMath::Pi() * ( i - FFTre.rows() / 2 ) * j / FFTre.rows() );

	//Array< complex< double >, 3 > Fw( FFTre.rows(), FFTre.cols(), upsampling );
	//Fw = sum( FFTre(i,j,l) * E(l,k), l );

	//Array< complex< double >, 3 > Fv( FFTre.rows(), upsampling, upsampling );
	//Fv = sum( Fw(i,l,j) * E(l,k), l );

	//Array< complex< double >, 3 > F( upsampling, upsampling, upsampling );
	//F = sum( Fv(l,i,j) * E(l,k), l );
	//F = conj(F);

	//TinyVector< double, Dim > cccPeakPositionSuper = this->locateMaxima( real(F) );
	//cccPeakPositionSuper -= ( dftshift + 1 );
	//shift += ( cccPeakPositionSuper / ( 1.0 * upsampling ) );

	//cdata->Delete();

//#ifdef NBF_VERBOSE
	// cout << " d=" << this->correlationPeak << ", t=[" << -cccPeakPosition[0] << ", " << cccPeakPosition[1] << ", " << -cccPeakPosition[2] << "]"<< endl;
//#endif

	// apply shifts to input transform

	vtkTransform * result = vtkTransform :: New();
	if ( Dim < 3 ){
		result->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], 0 );
	} else {
		result->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], - cccPeakPosition[thirdDim] );
	}

	vtkMatrix4x4 * concatenation = vtkMatrix4x4 :: New();
	vtkMatrix4x4 :: Multiply4x4( transform->GetMatrix(), result->GetMatrix(), concatenation );

	transform->SetInput( (vtkTransform*)NULL );
	transform->SetMatrix( concatenation );
	result->Delete();
	concatenation->Delete();
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: executeFourierNewHalf2D( vtkTransform * transform )
{

	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	// initialize image mask
	int dims[3];
	this->input1->GetDimensions(dims);
	TinyVector< int, Dim > s = this->imageFilter->getPaddingFactor() * TinyVector< int, Dim >( dims[0], dims[1] );
	if ( this->mask.size() == 0 ){
		this->reinitializeHalf(s);
	}

	// create handles to half-size FFTs
	TinyVector< int, 2 > shape( this->mask.rows(), this->mask.cols() / 2 + 1 );
	Array< Pixel, 2 > shiftr( this->fourierFilter->shift( Range :: all(), Range(fromStart,shape[1]-1) ) );
	Array< complex< double >, 2 > FFTre( reinterpret_cast<complex<double>*>(this->fourierFilter->blitzFFT.data()), shape, neverDeleteData );

	// create handles to half-size big FFTs
	TinyVector< int, 2 > shapebig( 2 * this->mask.rows(), 2 * this->mask.cols() / 2 + 1 );
	Array< complex< double >, 2 > FFTrebig( reinterpret_cast<complex<double>*>(this->fourierFilter->blitzFFTbig.data()), shapebig, neverDeleteData );

	Array< double, Dim > A;
	if ( this->useMissingWedgeCompensation == true ){
		if ( this->input1Changed == true ){
			this->wedge1.resize( s );
			this->wedgedInput1->getWedgeImageHalf( this->wedge1 );
			//w.write(this->wedge1);
		}

		this->wedge2.resize( s );
		this->wedgedInput2->getWedgeImageHalf( this->wedge2, transform );
		//w.write(this->wedge2);
		this->wedge2 *= this->wedge1;
		this->fourierFilter->wedgeOn( this->wedge2 );
	} else {
		this->wedge1.free();
		this->wedge2.free();
		this->fourierFilter->wedgeOff();
	}

	// compute \int ( I_1 - I_2 )^2 \times m_1 in reciprocal space

	// compute FFT of I_1 alone
	nbfVTKInterface :: vtkToBlitzReference( this->input1, A );
	// w.write(A);
	this->fourierFilter->blitzFFTreal = A * this->maskfft * this->fourierFilter->shift;
	fftw_execute( this->fourierFilter->fftplanreal );
	this->input1Changed = false;

	//Array< double, 3 > H( FFTre.shape() );
	//H = real( FFTre );
	//w.write(H);

	// apply band pass + wedge filter + mutual correlation scaling
	this->fourierFilter->executeHalf( FFTre, this->useMutualCorrelation );
	
	// apply normalization in reciprocal space
	Pixel k1 = this->normalizeFourierHalf( FFTre, false );
	
	// invert FFT
	fftw_execute( this->fourierFilter->ifftplanreal );

	Pixel k1_m = sum( this->fourierFilter->blitzFFTreal * this->fourierFilter->blitzFFTreal * this->referenceMask );

	// apply mask in real space
	this->fourierFilter->blitzFFTreal *= this->referenceMask;

	// compute final FFT of I_1 \times m_1
	fftw_execute( this->fourierFilter->fftplanreal );
	this->FFT1re.resize( FFTre.shape() );
	this->FFT1re = FFTre;

	// compute FFT of I_2 alone
	//this->wedgedInput2->getImage( this->input2, transform );
	this->input2Changed = false;
	nbfVTKInterface :: vtkToBlitzReference( this->input2, A );
	// w.write(A);
	this->fourierFilter->blitzFFTreal = A * this->maskfft * this->fourierFilter->shift;
	fftw_execute( this->fourierFilter->fftplanreal );

	//H = real( FFTre );
	//w.write(H);

	// apply band pass + wedge filter + mutual correlation scaling
	this->fourierFilter->executeHalf( FFTre, this->useMutualCorrelation );

	// apply normalization in reciprocal space
	Pixel k2 = this->normalizeFourierHalf( FFTre, false );

	// store normalized FFT of I_2 (fix fftw scaling)
	this->FFT2re.resize( FFTre.shape() );
	double normf = this->fourierFilter->shift.size();
	this->FFT2re = 2.0 * this->FFT1re * conj( FFTre * normf );

	// invert FFT
	fftw_execute( this->fourierFilter->ifftplanreal );

	//H = this->fourierFilter->blitzFFTreal * this->fourierFilter->shift;
	//w.write(H);

	// FFT of cuadratic term I_2^2
	this->fourierFilter->blitzFFTreal = this->fourierFilter->blitzFFTreal * this->fourierFilter->blitzFFTreal * this->fourierFilter->shift;
	fftw_execute( this->fourierFilter->fftplanreal );
	this->FFT1re.resize( FFTre.shape() );
	this->FFT1re = FFTre;

	// compute FFT of m_1
	this->fourierFilter->blitzFFTreal = this->referenceMask * this->fourierFilter->shift;
	fftw_execute( this->fourierFilter->fftplanreal );
	
	// evaluate quadratic expression in reciprocal space
	FFTre = ( - FFTre * conj( this->FFT1re ) + this->FFT2re ) * cast<double>( shiftr );

	double factor = ( 1.0 * FFTrebig.rows() * 2 * ( FFTrebig.cols() - 1 ) ) / ( 1.0 * FFTre.rows() * 2 * ( FFTre.cols() - 1 ) );

	// normalized cross correlation
	//FFTre = this->FFT1re * conj( FFTre ) * cast<double>( shiftr );
	
	Pixel upsampling_factor = 2.0;
	
	// put into up-sampled fft
	FFTrebig = 0;
	FFTrebig( Range( A.rows() - floor(A.rows()/2.0), A.rows() + floor((A.rows()-1)/2.0) ), Range( FFTre.cols() - 1, FFTrebig.ubound(secondDim) ) ) = FFTre * factor;
	fftw_execute( this->fourierFilter->ifftplanrealbig );
	
	// locate maxima in upsampled correlation volume

	//Array< double, Dim > Tf( 1.0 + this->fourierFilter->blitzFFTrealbig * this->fourierFilter->shiftbig / this->fourierFilter->shiftbig.size() );
	//w.write(Tf);

	// restrict search to admisible range
	TinyVector< Pixel, Dim > center = this->fourierFilter->blitzFFTrealbig.shape() / 2.0;
	Pixel restriction = center[0];
	if ( this->restrictTranslationSearch < numeric_limits< Pixel > :: max() ){
		restriction = blitz :: extrema :: max( this->restrictTranslationSearch * upsampling_factor, 1.0 );
	}
	Range Ix( blitz :: extrema :: max( 0, center[0] - floor(restriction) ), blitz :: extrema :: min( center[0] + floor(restriction), this->fourierFilter->blitzFFTrealbig.ubound(firstDim) ) );
	Range Iy( blitz :: extrema :: max( 0, center[1] - floor(restriction) ), blitz :: extrema :: min( center[1] + floor(restriction), this->fourierFilter->blitzFFTrealbig.ubound(secondDim) ) );
	////Array< double, Dim > T( 1.0 + this->fourierFilter->blitzFFTrealbig(I,I,I) * this->fourierFilter->shiftbig(I,I,I) / this->fourierFilter->shiftbig.size() );
	Array< double, Dim > T( this->fourierFilter->blitzFFTrealbig(Ix,Iy) * this->fourierFilter->shiftbig(Ix,Iy) / this->fourierFilter->shiftbig.size() );
	//w.write(T);

	//Array< double, Dim > T( ( - k1 - k2 + 2.0 * this->fourierFilter->blitzFFTreal(I,I,I) * this->fourierFilter->shift(I,I,I) ) / A.size() );
	Array< Pixel, Dim > radius( T.shape() );
	firstIndex i; secondIndex j; thirdIndex k;
	if ( Dim == 2 ) radius = ( restriction - i ) * ( restriction - i ) + ( restriction - j ) * ( restriction - j );
	if ( Dim == 3 ) radius = ( restriction - i ) * ( restriction - i ) + ( restriction - j ) * ( restriction - j ) + ( restriction - k ) * ( restriction - k );
	if ( k1_m < max( T ) ){
		T = -T;
	}
	T = where( radius < restriction * restriction, T, - numeric_limits< double > :: max() );
	//w.write(T);
	if ( ( this->overlapNormalizedDistances == true ) && ( this->fourierFilter->wedge.size() > 0 ) ){
		this->correlationPeak = max( 0.0, ( k1_m - max( T ) ) ) / sum( this->fourierFilter->wedge * this->fourierFilter->filter ) * sum( this->fourierFilter->filter );
		//this->correlationPeak = max( 0.0, 2.0 * ( 2.0 - max( T ) ) ) / sum( this->fourierFilter->wedge * this->fourierFilter->filter ) * sum( this->fourierFilter->filter );
		//this->correlationPeak = max( 0.0, - max( T ) ) / sum( this->fourierFilter->wedge * this->fourierFilter->filter ) * sum( this->fourierFilter->filter );
	} else {
		this->correlationPeak = max( 0.0, ( k1_m - max( T ) ) );
		//this->correlationPeak = max( 0.0, 2.0 * ( 2.0 - max( T ) ) );
		//this->correlationPeak = max( 0.0, - max( T ) );
	}
	this->correlationPeak = this->correlationPeak / this->fourierFilter->shiftbig.size() / this->fourierFilter->shiftbig.size();
	this->correlationScale = 1.0;

	TinyVector< double, Dim > cccPeakPosition = this->locateMaxima( T );
	cccPeakPosition = cccPeakPosition - ( T.shape() - 1 ) / 2.0;
	cccPeakPosition /= upsampling_factor;

	// apply shifts to input transform

	vtkTransform * result = vtkTransform :: New();
	if ( Dim < 3 ){
		result->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], 0 );
	} else {
		result->Translate( - cccPeakPosition[firstDim], - cccPeakPosition[secondDim], - cccPeakPosition[thirdDim] );
	}

	vtkMatrix4x4 * concatenation = vtkMatrix4x4 :: New();
	vtkMatrix4x4 :: Multiply4x4( transform->GetMatrix(), result->GetMatrix(), concatenation );
	transform->SetInput( (vtkTransform*)NULL );
	transform->SetMatrix( concatenation );
	result->Delete();
	concatenation->Delete();
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: executeFourierComplex( vtkTransform * transform )
{
	if ( transform != NULL ){
		this->wedgedInput2->getImage( this->input2, transform );
		this->input2Changed = true;
	}

	if ( this->input1Changed && this->input1ChangedLocally ){
		// build wedge filter
		this->wedge1.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
		this->wedgedInput1->getWedgeImage( this->wedge1 );
	}
	if ( this->input2Changed ){
		this->wedge2.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput2->getDimensions() );
		this->wedgedInput2->getWedgeImage( this->wedge2, transform );

		this->wedge2 *= this->wedge1;

		// if overlap is almost entire volume ignore wedge
		this->fourierFilter->wedgeOn( this->wedge2 );
	}

	Array< double, Dim > A;

	if ( this->input1Changed == true ){
		// window data in real space
		this->imageFilter->execute( this->input1 );

		nbfVTKInterface::vtkToBlitzReference( this->input1, A );

		this->fourierFilter->initializeFFTW( A.shape() );
		real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
		imag( this->fourierFilter->blitzFFT1 ) = 0;

		fftw_execute( this->fourierFilter->fftplan );

		// store FFT1 to avoid recomputation
		this->FFT1saved.resize( this->fourierFilter->blitzFFT2.shape() );
		this->FFT1saved = this->fourierFilter->blitzFFT2;
		this->input1Changed = false;
	}

	if ( this->input2Changed == true ){
	
		// window data in real space
		this->imageFilter->execute( this->input2 );

		nbfVTKInterface::vtkToBlitzReference( this->input2, A );
		this->fourierFilter->initializeFFTW( TinyVector< int, 3 >( A.shape() ) );

		real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
		imag( this->fourierFilter->blitzFFT1 ) = 0;

		fftw_execute( this->fourierFilter->fftplan );

		this->input2Changed = false;
	}

	this->fourierFilter->updateFilter();
	
	//w.write(this->wedge1);
	//w.write(this->wedge2);
	//
	//this->wedge1 = sqrt( real( this->FFT1saved * conj( this->FFT1saved ) ) );
	//w.write( this->wedge1 );
	//
	//this->wedge1 = sqrt( real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ) );
	//w.write( this->wedge1 );

	this->fourierFilter->blitzFFT2 = this->fourierFilter->blitzFFT2 - this->FFT1saved;
	this->correlationPeak = sum( sqrt( real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ) ) * this->wedge2 * this->fourierFilter->filter ) / sum( this->wedge2 * this->fourierFilter->filter );
}

/*
template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: execute( vtkTransform * transform )
{
	this->wedgedInput2->getImage( this->input2, transform );
	this->input2Changed = true;

	Array< Pixel, Dim > A;

	//Array< Pixel, Dim > I1, I2;
	//nbfVTKInterface::vtkToBlitz( this->input1, I1 );
	//nbfVTKInterface::vtkToBlitz( this->input2, I2 );
	//this->correlationPeak = sum( abs( I1 - I2 ) );

	//return;

	if ( this->input1Changed == true ){
		
		// window data in real space
		this->imageFilter->execute( this->input1 );

		nbfVTKInterface::vtkToBlitzReference( this->input1, A );

		this->fourierFilter->initializeFFTW( A.shape() );
		real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
		imag( this->fourierFilter->blitzFFT1 ) = 0;

		fftw_execute( this->fourierFilter->fftplan );

		this->modSqrFFT1.resize( this->fourierFilter->blitzFFT2.shape() );
		this->modSqrFFT1 = real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) );

		this->wedge1.resize( this->wedgedInput1->getDimensions() );
		this->wedgedInput1->getWedgeImage( this->wedge1 );

		this->input1Changed = false;
	}

	if ( this->input2Changed == true ){	
		// window data in real space
		this->imageFilter->execute( this->input2 );

		nbfVTKInterface::vtkToBlitzReference( this->input2, A );
		this->fourierFilter->initializeFFTW( A.shape() );

		real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
		imag( this->fourierFilter->blitzFFT1 ) = 0;

		fftw_execute( this->fourierFilter->fftplan );

		this->modSqrFFT2.resize( this->fourierFilter->blitzFFT2.shape() );
		this->modSqrFFT2 = real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) );

		this->wedge2.resize( this->wedgedInput2->getDimensions() );
		this->wedgedInput2->getWedgeImage( this->wedge2, transform );

		this->input2Changed = false;
	}

	//Pixel normalization = sum( this->wedge2 ) / ( 1.0 * this->wedge1.size() );
	//Pixel normalization = sum( where( abs( this->wedge2 ) > 0, 1, 0 ) ) / ( 1.0 * this->wedge1.size() );
	TinyVector< int, 3 > center = floor( this->wedge1.shape() / 2.0 );
	this->modSqrFFT1( center ) = 0;
	this->modSqrFFT2( center ) = 0;

	// compute x
	Array< Pixel, 3 > X( this->modSqrFFT1.shape() );
	//X = abs( log(1+sqrt(this->modSqrFFT1)) - log(1+sqrt(this->modSqrFFT2)) ) * this->wedge2 * this->wedge1;
	//X = pow2( log(1+sqrt(this->modSqrFFT1)) - log(1+sqrt(this->modSqrFFT2)) ) * this->wedge2 * this->wedge1;
	X = abs( sqrt(this->modSqrFFT1) - sqrt(this->modSqrFFT2) ) * this->wedge1 * this->wedge2;
	
	// median(x)
	vector< Pixel > vectorized;
	Array< Pixel, 3 > :: iterator iter = X.begin();
	while ( iter != X.end() ){
		if ( (*iter) > 0 ){
			vectorized.push_back( *iter );
		}
		++iter;
	}
	sort( vectorized.begin(), vectorized.end() );

	// compute mad
	Array< Pixel, 3 > Y( X.shape() );
	Y = where( X != 0, abs( X - vectorized[ floor( vectorized.size() / 2.0 ) ] ), 0 );

	iter = Y.begin();
	vectorized.clear();
	while ( iter != Y.end() ){
		if ( (*iter) != 0 ){
			vectorized.push_back( *iter );
		}
		++iter;
	}
	sort( vectorized.begin(), vectorized.end() );

	Pixel sigma = 1.4826 * vectorized[ floor( vectorized.size() / 2.0 ) ];
	
	this->fourierFilter->updateFilter();

	Pixel normalization = sum( where( this->wedge1 * this->wedge2 > 0, 1, 0 ) ) / ( 1.0 * this->wedge1.size() );

	this->wedge2 = this->wedge1 * this->wedge2 * this->fourierFilter->filter;

	// tukey - errors bigger than sigma are ignored, so if many pixels have big error (because they are discarded)
	// the value of the distance metric will be small!
	// X = where( abs(X) > sigma, 0, X * pow2( 1.0 - X * X / sigma / sigma ) ) * this->wedge2;
	
	// huber
	X = where( abs(X) > sigma, 1, X / sigma ) * this->wedge2;
	this->correlationPeak = sum( abs( X ) ) / normalization;
}
*/

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: refine( vtkTransform * t, Pixel tolerance )
{
	// set seed transform
	this->setSeedTransform( t );

	vnl_vector< double > seed(3);
	seed[0] = seed[1] = seed[2] = 0;

	//vtkTransform * t0 = vtkTransform::New();
	//this->executeFourier(t0);
	//t0->Delete();
	//Pixel maxPeak = this->getCorrelationPeak();
	//int maxDim = 0;
	//int maxOffset = 0;

	//for ( int dim = 0; dim < Dim; dim++ ){
	//	for ( int i = -1; i < 2; i+=2 ){
	//		vtkTransform * t = vtkTransform::New();
	//		switch ( dim ){
	//			case firstDim:
	//              t->RotateX( i );
	//				break;
	//			case secondDim:
	//				t->RotateY( i );
	//				break;
	//			case thirdDim:
	//				t->RotateZ( i );
	//				break;
	//		}
	//		this->executeFourier(t);
	//		t->Delete();
	//		if ( this->getCorrelationPeak() < maxPeak ){
	//			maxPeak = this->getCorrelationPeak();
	//			maxDim = dim; maxOffset = i;
	//		}
	//	}
	//}
	//seed[maxDim] = maxOffset;

	//this->lastX = this->lastY = this->lastZ = this->lastDistance = - numeric_limits< Pixel > :: max();
	//this->oneBeforeLastX = this->oneBeforeLastY = this->oneBeforeLastZ = this->oneBeforeLastDistance = - numeric_limits< Pixel > :: max();

	// refine maximum estimate with powell algorithm
	vnl_powell optimizer( this );
	optimizer.set_max_function_evals(1);
#if _DEBUG
	optimizer.set_verbose(true);
	optimizer.set_trace(true);
#endif
	optimizer.set_linmin_xtol( 1e-2 );
	optimizer.set_initial_step( 5.0 ); // changed from .5
	optimizer.set_x_tolerance( 1e-1 );
	optimizer.set_f_tolerance( 1e-1 );
	
	// refine only if required
	if ( tolerance > 0 ){
		optimizer.minimize( seed );
	}

	//// refine maximum estimate with simplex algorithm
	//vnl_amoeba optimizer( *this );
	////optimizer.set_relative_diameter(1e-1);
	////optimizer.set_x_tolerance(.1);
	//optimizer.minimize( seed );

	if ( optimizer.get_end_error() != numeric_limits< Pixel > :: max() ){

		vtkTransform * vt = vtkTransform :: New();
		vt->SetInput( this->seedTransform );
		vt->RotateX( seed[0] );
		vt->RotateY( seed[1] );
		vt->RotateZ( seed[2] );

		// get final alignment
		//this->execute( vt );
		this->executeFourierNewHalf( vt );
		t->DeepCopy( vt );
		vt->Delete();
	
#if 0
		cerr << "R[" << seed[0] << "," << seed[1] << "," << seed[2] << "],D=" << this->correlationPeak << ", ";
		Pixel translation[3];
		t->GetPosition(translation);
		cerr << "T=[" << translation[0] << "," << translation[1] << "," << translation[2] << "], D=" << this->correlationPeak << endl;
#endif
	} else {
		// get final alignment (compute score directly, no shift, no rotation)
		//this->execute( t );
		this->executeFourierNewHalf( t );
		cerr << "WARNING - Could not find minimizer.\n Aligning ";
		if ( this->wedgedInput1->getTypeId() == NBF_WEDGED_SUB_IMAGE_3D ){
			cerr << reinterpret_cast< nbfWedgedSubImage3D< Pixel > * >( this->wedgedInput1 )->getFileName() << " to ";
		} else {
			cerr << "average to ";
		}
		if ( this->wedgedInput2->getTypeId() == NBF_WEDGED_SUB_IMAGE_3D ){
			cerr << reinterpret_cast< nbfWedgedSubImage3D< Pixel > * >( this->wedgedInput2 )->getFileName();
		} else {
			cerr << " average";
		}
		cerr << ", In " << __FILE__ << ":" << __LINE__ << endl;

		cerr << "R[" << seed[0] << "," << seed[1] << "," << seed[2] << "],D=" << this->correlationPeak << ", ";
		Pixel translation[3];
		t->GetPosition(translation);
		cerr << "T=[" << translation[0] << "," << translation[1] << "," << translation[2] << "], D=" << this->correlationPeak << endl;
	}
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: execute()
{
	// type check
	if ( this->wedgedInput2->getTypeId() != NBF_WEDGED_SUB_IMAGE_3D ){
		cerr << "ERROR - Second input is expected to be of type nbfWedgedSubImage3D. In " << __FILE__ << ", " << __LINE__ << endl;
		return;
	}

	vtkTransform * t = vtkTransform::New();

	reinterpret_cast< nbfWedgedSubImage3D< Pixel > * >(this->wedgedInput2)->getTransform(t);

	vtkTransform * store = vtkTransform::New();
	store->DeepCopy( t );
	
	reinterpret_cast< nbfWedgedSubImage3D< Pixel > * >(this->wedgedInput2)->setTransform( (vtkTransform*)NULL );

	this->refine(t);
	//this->input1Changed = false; this->input2Changed = false;

	// restore original transform
	reinterpret_cast< nbfWedgedSubImage3D< Pixel > * >(this->wedgedInput2)->setTransform( store );

	this->transform->DeepCopy(t);

	//// compute the actual incremental change
	//vtkMatrix4x4 * incremental = vtkMatrix4x4 :: New();
	//store->Inverse();
	//vtkMatrix4x4::Multiply4x4( store->GetMatrix(), t->GetMatrix(), incremental );
	//
	//store->Delete();

	//// set transform attribute
	//vtkTransform * tmp = vtkTransform::New();
	//tmp->Concatenate( incremental );
	//this->transform->DeepCopy( tmp );
	//tmp->Delete();

	// store results in attributes

	this->candidateTransforms.clear();
	this->candidateCorrelationPeaks.clear();
	this->candidateCorrelationScales.clear();
	this->candidateWedgeOverlaps.clear();

	// copy transformation matrix into blitz array for storage
	Array< double, 1 > currentMatrix(16);

	double matrix[16];
	vtkMatrix4x4 :: DeepCopy( matrix, this->transform->GetMatrix() );
	//vtkMatrix4x4 :: DeepCopy( matrix, incremental );
	//incremental->Delete();
	
	for ( int i = 0; i < 16; i++ ){
		currentMatrix(i) = matrix[i];
	}

	this->candidateTransforms.push_back( currentMatrix );

	this->candidateCorrelationPeaks.push_back( this->getCorrelationPeak() );
	this->candidateCorrelationScales.push_back( this->getCorrelationScale() );
	if ( this->computeOverlaps ){
		this->candidateWedgeOverlaps.push_back( this->getWedgeOverlap( t ) );
	} else {
		this->candidateWedgeOverlaps.push_back( 1.0 );
	}
	t->Delete();
}


template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: getLowDimensionRepresentation( Array< complex< Pixel >, 1 > & C, Array< Pixel, 1 > & W )
{
	Array< double, Dim > A;
	if ( this->input1Changed == true ){
		
		// window data in real space
		this->imageFilter->execute( this->input1 );

		nbfVTKInterface::vtkToBlitzReference( this->input1, A );

		// compute fourier transform
		this->fourierFilter->initializeFFTWhalf( A.shape() );
		this->fourierFilter->blitzFFTreal = A * this->fourierFilter->shift;
		
		fftw_execute( this->fourierFilter->fftplanreal );

		// store FFT1 to avoid recomputation
		this->FFT1savedreal.resize( this->fourierFilter->blitzFFT.shape() );
		this->FFT1savedreal = this->fourierFilter->blitzFFT / sqrt( 1.0 * A.size() );
		this->input1Changed = false;

		// build wedge filter
		this->wedge1.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
		this->wedgedInput1->getWedgeImageHalf( this->wedge1 );

	}
	//else {
	//	this->fourierFilter->blitzFFT2 = this->FFT1saved;
	//}

	this->fourierFilter->wedgeOff();
	this->fourierFilter->updateFilterHalf();

	vector< complex< double > > newImageRepresentation;
	vector< Pixel > newWedgeRepresentation;
	
	TinyVector< int, 3 > shape( A.rows(), A.cols(), A.depth() / 2 + 1 );
	Array< complex< double >, 3 > FFTre( reinterpret_cast<complex<double>*>(this->FFT1savedreal.data()), shape, neverDeleteData );
	typename Array< complex< double >, 3 > :: iterator iterImage = FFTre.begin();
	typename Array< Pixel, 3 > :: iterator iterFilter = this->fourierFilter->filter.begin();
	typename Array< Pixel, 3 > :: iterator iterWedge = this->wedge1.begin();
	
	// kill DC component explicitly
	TinyVector< int, Dim > center( A.rows() / 2, A.cols() / 2, A.depth() - 1 );
	this->fourierFilter->filter( center ) = 0;

	FFTre *= this->fourierFilter->filter;

	//nbfMatlabWriter w;
	//w.setFileName("wedge.mat");
	//w.write(this->wedge1);
	//w.setFileName("filter.mat");
	//w.write(this->fourierFilter->filter);
	//w.setFileName("image.mat");
	//this->wedge1=cast<float>(sqrt(real(this->FFT1saved*conj(this->FFT1saved))));
	//w.write(this->wedge1);

	// use the symmetry property of the FFT and only store half size vectors
	while ( iterFilter != this->fourierFilter->filter.end() ){
		if ( fabs(*iterFilter) > 0.5 ){
			newImageRepresentation.push_back( *iterImage );
			newWedgeRepresentation.push_back( *iterWedge );
		}
		++iterImage; ++iterFilter; ++iterWedge;
	}

	C.resize( newImageRepresentation.size() );
	W.resize( newWedgeRepresentation.size() );
	for ( int i = 0; i < newImageRepresentation.size(); i++ ){
		C(i) = newImageRepresentation[i];
		if ( this->useMutualCorrelation == true ){
			Pixel n = sqrt( sqrt( real( C(i) * conj(C(i)) ) ) );
			if ( n > 0 ){
				C(i) = C(i) / n;
			}
		}
		W(i) = newWedgeRepresentation[i] > .5;
	}
	//w.write(real(C));
	//w.write(imag(C));
	//w.write(W);
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: getLowDimensionRepresentationHalf( Array< complex< Pixel >, 1 > & C, Array< Pixel, 1 > & W, int fold )
{
	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");

	Array< double, Dim > A;

	TinyVector< int, Dim > s = this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions();

	// initialize image mask
	if ( this->mask.size() == 0 ){
		this->reinitializeHalf(s);
	}

	A.resize( s );

	// nbfVTKInterface :: vtkToBlitzReference( this->input1, A );

	TinyVector< int, 3 > shape( A.rows(), A.cols(), A.depth() / 2 + 1 );
	Array< complex< double >, 3 > FFTre( reinterpret_cast<complex<double>*>(this->fourierFilter->blitzFFT.data()), shape, neverDeleteData );
	Array< Pixel, 3 > shiftr( this->fourierFilter->shift( Range :: all(), Range :: all(), Range(fromStart,shape[2]-1) ) );

	// apply image filter by-passing the mask
	this->imageFilter->execute( this->input1, true );

	// apply symmetry
	this->imageFilter->symmetrize( this->input1, fold );

	// build wedge filter as norm of FFT
	//this->fourierFilter->blitzFFTreal = A * this->maskfft;
	//fftw_execute( this->fourierFilter->fftplanreal );
	//this->wedge1.resize( FFTre.shape() );
	//this->wedge1 = sqrt( sqrt( real( FFTre * conj(FFTre) ) ) );
	//this->wedge1 = this->wedge1 / max(this->wedge1);

	// get wedge filter from geometry
	this->wedge1.resize( s );
	if ( fold > 1 ){
		this->wedge1 = 1;
	} else {
		this->wedgedInput1->getWedgeImageHalf( this->wedge1 );
	}

	// apply filters in reciprocal space (bandpass + mutual)
	//this->fourierFilter->wedgeOff();
	//this->fourierFilter->executeHalf( FFTre, this->useMutualCorrelation );
	//this->normalizeFourierHalf( FFTre );

	//// compute fourier transform
	//fftw_execute( this->fourierFilter->ifftplanreal );
	//this->fourierFilter->blitzFFTreal *= this->mask;
	//fftw_execute( this->fourierFilter->fftplanreal );

	nbfVTKInterface :: vtkToBlitzReference( this->input1, A );

	//this->fourierFilter->wedgeOff();
	//this->fourierFilter->executeHalf( A );

	//this->fourierFilter->blitzFFTreal = A * this->maskfft;
	this->fourierFilter->blitzFFTreal = A * this->fourierFilter->shift;
	fftw_execute( this->fourierFilter->fftplanreal );

	//w.write(real(FFTre));

	this->fourierFilter->wedgeOff();
	this->fourierFilter->executeHalf( FFTre, this->useMutualCorrelation );
	fftw_execute( this->fourierFilter->ifftplanreal );
	this->fourierFilter->blitzFFTreal *= this->mask;
	fftw_execute( this->fourierFilter->fftplanreal );

	//this->fourierFilter->blitzFFTreal = A * this->mask * this->fourierFilter->shift;
	//fftw_execute( this->fourierFilter->fftplanreal );

	// FFTre = FFTre * cast<double>(shiftr);

	//this->fourierFilter->blitzFFTreal = A * this->mask * this->fourierFilter->shift;
	//fftw_execute( this->fourierFilter->fftplanreal );
	//this->fourierFilter->wedgeOff();
	//this->fourierFilter->executeHalf( FFTre, this->useMutualCorrelation );
	//FFTre = FFTre * cast<double>(shiftr);

	//w.write(real(FFTre));
	//w.write(this->wedge1);
	//w.write(this->fourierFilter->filter);

	vector< complex< double > > newImageRepresentation;
	vector< Pixel > newWedgeRepresentation;
	
	//// kill DC component explicitly
	//TinyVector< int, Dim > center( this->fourierFilter->filter.rows() / 2, this->fourierFilter->filter.cols() / 2, this->fourierFilter->filter.ubound(thirdDim) );
	//this->fourierFilter->filter( center ) = 0;

	typename Array< complex< double >, 3 > :: iterator iterImage = FFTre.begin();
	typename Array< Pixel, 3 > :: iterator iterFilter = this->fourierFilter->filter.begin();
	typename Array< Pixel, 3 > :: iterator iterWedge = this->wedge1.begin();
	
	//w.write(this->wedge1);
	//w.setFileName("filter.mat");
	//w.write(this->fourierFilter->filter);
	//w.setFileName("image.mat");
	//this->wedge1=cast<float>(sqrt(real(this->FFT1saved*conj(this->FFT1saved))));

	// store half size FFT vectors
	while ( iterFilter != this->fourierFilter->filter.end() ){
		if ( fabs(*iterFilter) > 0.5 ){
			newImageRepresentation.push_back( *iterImage );
			newWedgeRepresentation.push_back( *iterWedge );
		}
		++iterImage; ++iterFilter; ++iterWedge;
	}

	C.resize( newImageRepresentation.size() );
	W.resize( newWedgeRepresentation.size() );
	for ( int i = 0; i < newImageRepresentation.size(); i++ ){
		C(i) = newImageRepresentation[i];
		// W(i) = newWedgeRepresentation[i];
		W(i) = newWedgeRepresentation[i] > .5;
	}
	//w.write(real(C));
	//w.write(imag(C));
	//w.write(W);
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: get2DimensionRepresentationHalf( Array< complex< Pixel >, 1 > & C, Array< Pixel, 2 > & M, int imageIndex, bool shifted )
{
	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	// Setup 2D Fourier filter
	nbfFourierFilter< Pixel, 2 > fourierFilter2D;
	fourierFilter2D.bandPassOn( this->fourierFilter->getBandLowCut(), this->fourierFilter->getBandHighCut(), this->fourierFilter->getBandLowVariance(), this->fourierFilter->getBandHighVariance() );
	if ( this->fourierFilter->getBfactorExponent() != 0 ){
		fourierFilter2D.bfactorOn( this->fourierFilter->getBfactorExponent() );
	}
	
	char currentslice[400];
	sprintf( currentslice, "%s_%04d.mrc", reinterpret_cast<nbfWedgedSubImage3D<float>*>(this->wedgedInput1)->getFileName().c_str(), imageIndex );
	
	char command[400];
	sprintf( command, "newstack %s %s -secs %d >> /dev/null", reinterpret_cast<nbfWedgedSubImage3D<float>*>(this->wedgedInput1)->getFileName().c_str(), currentslice, imageIndex );
	system( command );
	
	nbfMrcReader reader;
	reader.setFileName( currentslice );
	vtkImageData * data = vtkImageData::New();
	reader.read( data );
	Array< Pixel, 3 > A;
	nbfVTKInterface :: vtkToBlitz( data, A );
	
	// delete temporal slice
	sprintf( command, "rm %s >> /dev/null", currentslice );
	system( command );
	
	// Setup 2D image filter
	nbfImageFilter< Pixel, 2 > imageFilter2D;
	imageFilter2D.paddingOn( this->imageFilter->getPaddingFactor() );
	int win = A.rows() / 2.0 - 10;
	imageFilter2D.setMaskSize( win, win, win, 5, false );

	TinyVector< int, Dim > s = this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions();
	s = this->imageFilter->getPaddingFactor() * A.shape();

	// // initialize image mask
	// if ( this->mask.size() == 0 ){
		// this->reinitializeHalf(s);
	// }

	TinyVector< int, 2 > s2Dhalf( A.rows(), A.cols() );
	imageFilter2D.buildMask( s2Dhalf );

	TinyVector< int, 2 > s2D( s[0], s[1] );
	fourierFilter2D.initializeFFTWhalf( s2D );
	// imageFilter2D.buildMask( s2D );
	
	// A.resize( s );

	TinyVector< int, 2 > shape( s[0], s[1] / 2 + 1 );
	Array< complex< double >, 2 > FFTre( reinterpret_cast<complex<double>*>(fourierFilter2D.blitzFFT.data()), shape, neverDeleteData );
	Array< Pixel, 2 > shiftr( fourierFilter2D.shift( Range :: all(), Range(fromStart,shape[1]-1) ) );

	Array< Pixel, 2 > Ip( A( Range :: all(), Range :: all(), 0 ) );
	Array< double, 2 > I( Ip.shape() );
	I = cast<double>(Ip);
	data->Delete();

	// apply image filter by-passing the mask
	imageFilter2D.execute( I, false );

	// Array< Pixel, 3 > test( I.rows(), I.cols(), 1 );
	// test( Range :: all(), Range :: all(), 0 ) = cast<Pixel>(I);
	// w.setFileName("imagefiltered.bin");
	// w.write(test);

	// // apply symmetry
	// imageFilter2D.symmetrize( I, fold );

	// build wedge filter as norm of FFT
	//this->fourierFilter->blitzFFTreal = A * this->maskfft;
	//fftw_execute( this->fourierFilter->fftplanreal );
	//this->wedge1.resize( FFTre.shape() );
	//this->wedge1 = sqrt( sqrt( real( FFTre * conj(FFTre) ) ) );
	//this->wedge1 = this->wedge1 / max(this->wedge1);

	// // get wedge filter from geometry
	// this->wedge1.resize( s );
	// if ( fold > 1 ){
		// this->wedge1 = 1;
	// } else {
		// this->wedgedInput1->getWedgeImageHalf( this->wedge1 );
	// }

	// apply filters in reciprocal space (bandpass + mutual)
	//this->fourierFilter->wedgeOff();
	//this->fourierFilter->executeHalf( FFTre, this->useMutualCorrelation );
	//this->normalizeFourierHalf( FFTre );

	//// compute fourier transform
	//fftw_execute( this->fourierFilter->ifftplanreal );
	//this->fourierFilter->blitzFFTreal *= this->mask;
	//fftw_execute( this->fourierFilter->fftplanreal );

	// nbfVTKInterface :: vtkToBlitzReference( this->input1, A );

	//this->fourierFilter->wedgeOff();
	//this->fourierFilter->executeHalf( A );

	
	//this->fourierFilter->blitzFFTreal = A * this->maskfft;
	if ( shifted == true ){
		fourierFilter2D.blitzFFTreal = I * fourierFilter2D.shift;
	} else {
		fourierFilter2D.blitzFFTreal = I;
	}
	fftw_execute( fourierFilter2D.fftplanreal );

	//w.write(real(FFTre));

	fourierFilter2D.wedgeOff();
	fourierFilter2D.executeHalf( FFTre, this->useMutualCorrelation );
	// fftw_execute( fourierFilter2D.ifftplanreal );
	// fourierFilter2D.blitzFFTreal *= imageFilter2D.mask;
	// fftw_execute( fourierFilter2D.fftplanreal );

	M.resize( FFTre.shape() );
	M = real( FFTre * conj( FFTre ) );
	// M = cast< Pixel >( I );
	
	nbfPolarDomain< Pixel, 2 > polar;
	TinyVector< Pixel , 2 > center( (M.rows()+1)/2, M.cols()-2 );
	polar.setCenter( center );
	polar.setMinRho( M.cols() * fourierFilter2D.getBandLowCut() );
	polar.setMaxRho( M.cols() * fourierFilter2D.getBandHighCut()  );
	polar.setResRho( 2 * M.cols() );
	polar.setResTheta( 360 );
	Array< Pixel, 2 > P;
	Array< bool, 2 > B;
	polar.cartesian2polar( M, P, B );
	Array< Pixel, 1 > Ps( P.rows() );
	secondIndex j;
	// Ps = sum( P, j );
	Ps = sum( P( Range :: all(), Range( Ps.cols()/2, toEnd ) ), j );

	//this->fourierFilter->blitzFFTreal = A * this->mask * this->fourierFilter->shift;
	//fftw_execute( this->fourierFilter->fftplanreal );

	// FFTre = FFTre * cast<double>(shiftr);

	//this->fourierFilter->blitzFFTreal = A * this->mask * this->fourierFilter->shift;
	//fftw_execute( this->fourierFilter->fftplanreal );
	//this->fourierFilter->wedgeOff();
	//this->fourierFilter->executeHalf( FFTre, this->useMutualCorrelation );
	//FFTre = FFTre * cast<double>(shiftr);

	//w.write(real(FFTre));
	//w.write(this->wedge1);
	//w.write(this->fourierFilter->filter);

	vector< complex< double > > newImageRepresentation;
	
	//// kill DC component explicitly
	//TinyVector< int, Dim > center( this->fourierFilter->filter.rows() / 2, this->fourierFilter->filter.cols() / 2, this->fourierFilter->filter.ubound(thirdDim) );
	//this->fourierFilter->filter( center ) = 0;

	// cout << "FFTre = " << FFTre.shape() << endl;
	// cout << "ffilter = " << fourierFilter2D.filter.shape() << endl;
	
	typename Array< complex< double >, 2 > :: iterator iterImage = FFTre.begin();
	typename Array< Pixel, 2 > :: iterator iterFilter = fourierFilter2D.filter.begin();
	
	//w.write(this->wedge1);
	// w.setFileName("filter.mat");
	// w.write(this->fourierFilter->filter);
	// w.setFileName("image.mat");
	// w.write(I);
	//this->wedge1=cast<float>(sqrt(real(this->FFT1saved*conj(this->FFT1saved))));

	// nbfMrcWriter m;
	// m.setFileName("filter.mrc");
	// Array< Pixel, 3 > F( fourierFilter2D.filter.rows(), fourierFilter2D.filter.cols(), 1 );
	// F( Range :: all(), Range :: all(), 0 ) = cast< Pixel >( fourierFilter2D.filter );
	// m.write(F);
	
	// store half size FFT vectors
	while ( iterFilter != fourierFilter2D.filter.end() ){
		if ( fabs(*iterFilter) > 0.5 ){
			newImageRepresentation.push_back( *iterImage );
		}
		++iterImage; ++iterFilter;
	}

	C.resize( newImageRepresentation.size() );
	for ( int i = 0; i < newImageRepresentation.size(); i++ ){
		C(i) = newImageRepresentation[i];
	}
	//w.write(real(C));
	//w.write(imag(C));
	
	C.resize( Ps.shape() );
	real(C) = Ps;
	imag(C) = 0;
}


template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: putLowDimensionRepresentationReal( Array< Pixel, 1 > & C, Pixel bin_factor )
{
	/*
	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");

	Array< double, 3 > A;
	nbfVTKInterface :: vtkToBlitzReference( this->input1, A );

	//TinyVector< int, Dim > s = A.shape();
	//this->reinitializeHalf( s );
	
	//this->fourierFilter->wedgeOff();
	//this->fourierFilter->executeHalf( A );
	//w.write(A);
	
	// switch temporarily to file mask
	this->imageFilter->fileMaskOn();
	this->imageFilter->execute( this->input1, true );

	vtkImageData * mask_bin_vtk = vtkImageData :: New();
	if ( this->imageFilter->mask.size() == 0 ){
		this->imageFilter->buildMask( s );
	}

	Array< Pixel, 3 > mask_bin;
	
	if ( bin_factor != 1 ){
		nbfVTKInterface :: blitzToVtk( this->imageFilter->mask, mask_bin_vtk );
		vtkImageResample * resample = vtkImageResample :: New();
		resample->SetInput( mask_bin_vtk );
		resample->SetAxisMagnificationFactor( 0, 1.0 / bin_factor );
		resample->SetAxisMagnificationFactor( 1, 1.0 / bin_factor );
		resample->SetAxisMagnificationFactor( 2, 1.0 / bin_factor );
		resample->SetInterpolationModeToCubic();
		resample->Update();

		nbfVTKInterface :: vtkToBlitz( resample->GetOutput(), mask_bin );
		mask_bin_vtk->Delete();

		resample->SetInput( this->input1 );
		resample->Update();
		nbfVTKInterface :: vtkToBlitz( resample->GetOutput(), A );

		resample->Delete();
	} else {
		nbfVTKInterface :: vtkToBlitz( this->input1, A );
		mask_bin.resize( this->imageFilter->mask.shape() );
		mask_bin = this->imageFilter->mask;
	}

	typename Array< double, 3 > :: iterator iter = A.begin();
	typename Array< Pixel, 3 > :: iterator miter = mask_bin.begin();
	typename Array< Pixel, 1 > :: iterator citer = C.begin();

	while ( iter != A.end() ){
		if ( *miter > .5 ){
			*iter = *citer;
			citer++;
		}
		++iter;
		++miter;
	}
	
	// unbinning
	if ( bin_factor != 1 ){
		vtkImageResample * resample = vtkImageResample :: New();
		resample->SetInput( A );
		resample->SetAxisMagnificationFactor( 0, bin_factor );
		resample->SetAxisMagnificationFactor( 1, bin_factor );
		resample->SetAxisMagnificationFactor( 2, bin_factor );
		resample->SetInterpolationModeToCubic();
		resample->Update();
		nbfVTKInterface :: vtkToBlitz( resample->GetOutput(), this->input1 );
		resample->Delete();
	}

	// revert sigma to previous value
	this->imageFilter->fileMaskOff();
	*/
}

template< class Pixel, int const Dim >
void nbfFourierImageMetric< Pixel, Dim > :: getLowDimensionRepresentationReal( Array< Pixel, 1 > & C, Pixel bin_factor, int fold )
{
	// wegde-aware symmetrization
	if ( fold > 1 ){

		// build wedged average with single volume
		vector< nbfWedgedSubImage3D< Pixel > > v;
		nbfWedgedSubImage3D< Pixel > * volume = reinterpret_cast< nbfWedgedSubImage3D< Pixel > * >( this->wedgedInput1 );
		v.push_back( *volume );
		nbfWedgedAverageImage3D< Pixel > average( v );
		
		// assign alignments
		Array< Pixel, 3 > referenceAlignments( v.size(), 17, 1 );
		referenceAlignments = 0;
		
		// set unit weight
		referenceAlignments( Range::all(), 0, 0 ) = 1;

		vtkTransform * t = vtkTransform::New();
		volume->getTransform(t);
		double matrix[16];
		vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
		t->Delete();
		for ( int e = 0; e < 16; e++ ){
			referenceAlignments( 0, 1 + e, 0 ) = matrix[e];
		}
		average.setAlignments( referenceAlignments );

		// f-fold-symmetrization
		Array< Pixel, 3 > newAlignments( average.weights.rows(), 17, fold );
		for ( int f = 0; f < fold; f++ ){

			newAlignments( Range::all(), 0, f ) = average.weights( Range::all(), 0 ) / fold;

			vtkTransform * t = vtkTransform :: New();
			t->RotateZ( f * 360.0 / fold );

			for ( int i = 0; i < average.weights.rows(); i++ ){

				// retrieve current transformation
				vtkMatrix4x4 * mat1 = vtkMatrix4x4 :: New();

				double mat[16];
				for ( int m = 0; m < 16; m++ ){
					mat[m] = average.multipleAlignments( i, m, 0 );
				}
				mat1->DeepCopy( mat );

				vtkMatrix4x4 * mat3 = vtkMatrix4x4 :: New();
				vtkMatrix4x4 :: Multiply4x4( mat1, t->GetMatrix(), mat3 );

				vtkMatrix4x4::DeepCopy( mat, mat3 );
				for ( int m = 0; m < 16; m++ ){
					newAlignments( i, 1 + m, f ) = mat[m];
				}
				mat3->Delete();
				mat1->Delete();
			}

			t->Delete();
		}
		average.setAlignments( newAlignments );
		average.getImage( this->input1 );
	}
	////////////////////////////////////////			
	
	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");

	Array< double, 3 > A;
	nbfVTKInterface :: vtkToBlitzReference( this->input1, A );

	TinyVector< int, Dim > s = A.shape();
	this->reinitializeHalf( s );
	//A = ( A - mean(A) ) * this->maskfft * this->fourierFilter->shift;
	//w.write(A);
	//A = A * this->maskfft * this->fourierFilter->shift;
	
	this->fourierFilter->wedgeOff();
	this->fourierFilter->executeHalf( A );
	//w.write(A);
	
	//// store current sigma value
	//Pixel sigma = this->imageFilter->getSigma();
	//this->imageFilter->setSigma(0.0);

	// switch temporarily to file mask
	this->imageFilter->fileMaskOn();
	this->imageFilter->execute( this->input1, true );

	this->imageFilter->symmetrize( this->input1, fold );
	//nbfVTKInterface :: vtkToBlitzReference( this->input1, A );
	//w.write(A);

	//nbfVTKInterface :: vtkToBlitz( this->input1, A );
	//Array< double, 3 > Abin( A.shape() / 2 );
	//for ( int i = 0; i < Abin.rows(); i++ ){
	//	for ( int j = 0; j < Abin.cols(); j++ ){
	//		for ( int k = 0; k < Abin.depth(); k++ ){
	//			Abin(i,j,k) = A(2*i,2*j,2*k) + A(2*i+1,2*j,2*k) + A(2*i,2*j+1,2*k) + A(2*i,2*j,2*k+1) + 
	//				A(2*i+1,2*j+1,2*k) + A(2*i+1,2*j,2*k+1) + A(2*i,2*j+1,2*k+1) + A(2*i+1,2*j+1,2*k+1);
	//		}
	//	}
	//}
	//nbfVTKInterface :: vtkToBlitz( this->input1, A );
	//w.write(A);
	
	vtkImageData * mask_bin_vtk = vtkImageData :: New();
	if ( this->imageFilter->mask.size() == 0 ){
		this->imageFilter->buildMask( s );
	}

	Array< Pixel, 3 > mask_bin;
	
	if ( bin_factor != 1 ){
		nbfVTKInterface :: blitzToVtk( this->imageFilter->mask, mask_bin_vtk );
		vtkImageResample * resample = vtkImageResample :: New();
		resample->SetInput( mask_bin_vtk );
		resample->SetAxisMagnificationFactor( 0, 1.0 / bin_factor );
		resample->SetAxisMagnificationFactor( 1, 1.0 / bin_factor );
		resample->SetAxisMagnificationFactor( 2, 1.0 / bin_factor );
		resample->SetInterpolationModeToCubic();
		resample->Update();

		nbfVTKInterface :: vtkToBlitz( resample->GetOutput(), mask_bin );
		mask_bin_vtk->Delete();

		resample->SetInput( this->input1 );
		resample->Update();
		nbfVTKInterface :: vtkToBlitz( resample->GetOutput(), A );
		resample->Delete();
	} else {
		nbfVTKInterface :: vtkToBlitz( this->input1, A );
		mask_bin.resize( this->imageFilter->mask.shape() );
		mask_bin = this->imageFilter->mask;
	}

	vector< double > newImageRepresentation;

	typename Array< double, 3 > :: iterator iter = A.begin();
	typename Array< Pixel, 3 > :: iterator miter = mask_bin.begin();

	while ( iter != A.end() ){
		if ( *miter > .5 ){
			newImageRepresentation.push_back( *iter );
		}
		++iter;
		++miter;
	}

	C.resize( newImageRepresentation.size() );
	for ( int i = 0; i < newImageRepresentation.size(); i++ ){
		C(i) = newImageRepresentation[i];
	}
	
	//// revert sigma to previous value
	//this->imageFilter->setSigma( sigma );
	this->imageFilter->fileMaskOff();
}
