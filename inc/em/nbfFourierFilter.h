#pragma once

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkTransform.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkImageFourierCenter.h>

#include <io/nbfVTKInterface.h>

#include <fftw3.h>

using namespace blitz;

/** Fourier Filter 3D. For efficiency, filtering is conceived as an all-in-one filter. 
	Normally we apply several different filters in cascade, which corresponds to multiplication
	in reciprocal space. This class allows activation and deactivation of several filters to be
	applied simultaneously. Supported filters are low and high pass, and for dealing with the wedge.
	Bandpass also supported.
	Two different wedges can be applied simultaneously wedge1 and wedge2.
*/
template< class Pixel, const int Dim >
class nbfFourierFilter
{
public:

	nbfFourierFilter();

	~nbfFourierFilter();

	/// Activate use of low pass filter with X, Y and Z cuttof frequencies and order.
	void lowPassOn( TinyVector< Pixel, Dim > &, int = 2 );

	/// Deactive use of low pass filter.
	void lowPassOff(){ this->useLowPass = false; this->filterUpdated = false; this->filterUpdatedHalf = false; }

	/// Activate use of high pass filter with X, Y and Z cuttof frequencies and order.
	void highPassOn( TinyVector< Pixel, Dim > &, int = 2 );

	/// Deactivate use of high pass filter.
	void highPassOff(){ this->useHighPass = false; this->filterUpdated = false; this->filterUpdatedHalf = false; }

	/// Activate use of band pass filter with lowCut, highCut, lowVariance and highVariance.
	void bandPassOn( Pixel, Pixel, Pixel, Pixel );

	/// Deactivate use of high pass filter.
	void bandPassOff(){ this->useBandPass = false; this->filterUpdated = false; this->filterUpdatedHalf = false; }	

	/// Activate use of gaussian filter with center frequency and variance.
	void gaussianOn( Pixel, Pixel );

	/// Deactivate use of gaussian filter.
	void gaussianOff(){ this->useGaussian = false; this->filterUpdated = false; this->filterUpdatedHalf = false; }

	/// Activate use of exponential bfactor.
	void bfactorOn( Pixel p ){ this->bfactorExponent = p; this->filterUpdated = false; this->filterUpdatedHalf = false; }

	/// Deactivate use of exponential bfactor.
	void bfactorOff(){ this->bfactorExponent = 0; this->filterUpdated = false; this->filterUpdatedHalf = false; }

	void wedgeOn( Array< Pixel, Dim > & W ){ this->wedge.reference(W); this->dummyFilter = false; }
	void wedgeOff(){ this->wedge.free(); }

	/// Apply filter to real space input data. Arguments: (input,output). 
	void execute( vtkImageData * );
	void execute( Array< double, Dim > &, int = 0 );

	/// Apply filter to reciprocal space input data. Arguments: (input,output)
	void execute( Array< Pixel, Dim > &, Array< Pixel, Dim > & );
	void execute( Array< complex< double >, Dim > & );

	void executeHalf( Array< complex< double >, Dim > &, bool );
	void executeHalf( Array< double, Dim > &, int = 0 );
	void executeHalf( Array< Pixel, Dim > &, Array< Pixel, Dim > & );


	void  fft( Array< complex< double >, Dim > &, Array< complex< double >, Dim > & );
	void  fft( Array< Pixel, Dim > &, Array< complex< double >, Dim > & );
	void ifft( Array< complex< double >, Dim > &, Array< complex< double >, Dim > & );

	/// Make available for direct FFT manipulation
	Array< complex< double >, Dim > blitzFFT1;
	Array< complex< double >, Dim > blitzFFT2;
	fftw_plan  fftplan;
    fftw_plan ifftplan;

	Array< double, Dim > blitzFFT;						// internal FFTW in-place array
	Array< double, Dim > blitzFFTreal;					// real view
	Array< double, Dim > blitzFFTbig;						// internal FFTW in-place array
	Array< double, Dim > blitzFFTrealbig;					// real view
	fftw_plan  fftplanreal;
    fftw_plan ifftplanreal;
    fftw_plan ifftplanrealbig;

	void initializeFFTW( const TinyVector< int, Dim > & );
	void initializeFFTWhalf( const TinyVector< int, Dim > &, bool = false );

	void centerFFT( vtkImageData * );

	template< class otherPixel >
	void centerFFT( Array< otherPixel, Dim > & );

	Array< Pixel, Dim > shift;
	Array< Pixel, Dim > shiftbig;
	Array< Pixel, Dim > shiftreal;

	// store LPF * HPF * BPF
	Array< Pixel, Dim > filter;
	Array< Pixel, Dim > wedge;

	void updateFilter();
	void updateFilterHalf();

	Pixel getBandLowCut(){ if ( this->useBandPass ){ return this->bandLowCut; } else { return 0; } }
	Pixel getBandHighCut(){ if ( this->useBandPass ){ return this->bandHighCut; } else { return 1;} }
	Pixel getBandLowVariance(){ return this->bandLowVariance; }
	Pixel getBandHighVariance(){ return this->bandHighVariance; }
	Pixel getBfactorExponent(){ return this->bfactorExponent; }
	
protected:

	void buildLowPass();
	void buildHighPass();
	void buildBandPass();
	void buildBFactor();
	void buildGaussian();

	void buildLowPassHalf(){};
	void buildHighPassHalf(){};
	void buildBandPassHalf();
	void buildBFactorHalf();
	void buildGaussianHalf(){};

	Array< Pixel, Dim > lowPass;
	Array< Pixel, Dim > highPass;
	Array< Pixel, Dim > bandPass;
	Array< Pixel, Dim > bfactor;
	Array< Pixel, Dim > gaussian;

	/// Used to shift filter image to match that of the FFT geometry.
	TinyVector< int, Dim > center;

	bool useLowPass, useHighPass, useBandPass, useGaussian, filterUpdated, filterUpdatedHalf;

	bool dummyFilter;

	Pixel lowOrder;
	TinyVector< Pixel, Dim > lowCut;

	Pixel highOrder;
	TinyVector< Pixel, Dim > highCut;
	
	double bandLowCut, bandHighCut, bandLowVariance, bandHighVariance;
	double gaussianMean, gaussianVariance;

	double bfactorExponent;

	// store filter geometry
	TinyVector< int, Dim > size;

	// FFTW structures
	bool fftwNotInitialized, fftwNotInitializedHalf;

	fftw_complex * inFFT;
	fftw_complex * outFFT;
	double * FFTreal;
	double * FFTrealbig;

	vtkImageFourierCenter * centerfft;
};


template< class Pixel, const int Dim >
nbfFourierFilter< Pixel, Dim > :: nbfFourierFilter()
: useLowPass(false), useHighPass(false), useBandPass(false), useGaussian(false), filterUpdated(false), filterUpdatedHalf(false), fftwNotInitialized(true),fftwNotInitializedHalf(true), dummyFilter(true)
{ 
	this->size = 0;
	this->lowOrder = 2;
	this->lowCut = .5;
	this->highOrder = 2;
	this->highCut = .5;
	this->bandLowCut = 0;
	this->bandHighCut = 1;
	this->bandLowVariance = .1;
	this->bandHighVariance = .1;
	this->gaussianMean = .5;
	this->gaussianVariance = .1;
	this->bfactorExponent = 0;
	// calling New() statically causes a crash before main()
	//this->centerfft = vtkImageFourierCenter::New();
	this->centerfft = NULL;
}

template< class Pixel, const int Dim >
nbfFourierFilter< Pixel, Dim > :: ~nbfFourierFilter()
{ 
	if ( this->fftwNotInitialized == false ){
		fftw_destroy_plan( this->fftplan );
		fftw_destroy_plan( this->ifftplan );
	}
	if ( this->fftwNotInitializedHalf == false ){
		fftw_destroy_plan( this->fftplanreal );
		fftw_destroy_plan( this->ifftplanreal );
		fftw_destroy_plan( this->ifftplanrealbig );
	}
	if (this->centerfft != NULL) {
		this->centerfft->Delete();
	}
}

template< class Pixel, int const Dim  >
void nbfFourierFilter< Pixel, Dim > :: initializeFFTW( const TinyVector< int, Dim > & A )
{
	if ( ( sum( fabs( this->size - A ) ) != 0 ) || this->fftwNotInitialized ){

		this->size = A;
		
		// check if all dimensions are even (the fftshift as implemented requires this)
		for ( int i = 0; i < Dim; i++ ){
			if ( this->size(i) % 2 != 0 ){
				cerr << "ERROR" << __FILE__ << " : " << __LINE__ << endl;
				cerr << "Expecting even image dimensions." << endl;
				assert(0);
			}
		}

		this->blitzFFT1.resize( this->size );
		this->blitzFFT2.resize( this->size );
		this->shift.resize( this->size );
		this->shiftbig.resize( 2 * this->size );

		typename Array< Pixel, Dim > :: iterator iter = this->shift.begin();
		while ( iter != this->shift.end() ){
			(*iter) = pow( -1.0, (int)sum( iter.position() ) );
			++iter;
		}

		iter = this->shiftbig.begin();
		while ( iter != this->shiftbig.end() ){
			(*iter) = pow( -1.0, (int)sum( iter.position() ) );
			++iter;
		}

		this->inFFT  = reinterpret_cast<fftw_complex*>( this->blitzFFT1.data() );
		this->outFFT  = reinterpret_cast<fftw_complex*>( this->blitzFFT2.data() );
		
		if ( this->fftwNotInitialized == false ){
			fftw_destroy_plan( this->fftplan );
			fftw_destroy_plan( this->ifftplan );
		}

		if ( Dim == 2 ){
			this->fftplan = fftw_plan_dft_2d( size[0], size[1], this->inFFT, this->outFFT, FFTW_FORWARD, FFTW_MEASURE );
			this->ifftplan = fftw_plan_dft_2d( size[0], size[1], this->outFFT, this->inFFT, FFTW_BACKWARD, FFTW_MEASURE );
		}
		if ( Dim == 3 ){
			this->fftplan = fftw_plan_dft_3d( size[0], size[1], size[2], this->inFFT, this->outFFT, FFTW_FORWARD, FFTW_MEASURE );
			this->ifftplan = fftw_plan_dft_3d( size[0], size[1], size[2], this->outFFT, this->inFFT, FFTW_BACKWARD, FFTW_MEASURE );
		}

		this->fftwNotInitialized = false;
		this->fftwNotInitializedHalf = true;

	}
}

template< class Pixel, int const Dim  >
void nbfFourierFilter< Pixel, Dim > :: initializeFFTWhalf( const TinyVector< int, Dim > & A, bool ignoreBig )
{
	if ( sum( fabs( this->size - A ) ) != 0 ){

		this->size = A;
		
		// check if all dimensions are even (the fftshift as implemented requires this)
		for ( int i = 0; i < Dim; i++ ){
			if ( this->size(i) % 2 != 0 ){
				cerr << "ERROR" << __FILE__ << " : " << __LINE__ << endl;
				cerr << "Expecting even image dimensions." << endl;
				assert(0);
			}
		}

		if ( Dim == 2 ){
			this->blitzFFT.resize( this->size[0], 2 * ( this->size[1] / 2 + 1 ) );
			this->blitzFFTreal.reference( this->blitzFFT( Range :: all(), Range(fromStart, this->size[1] - 1 ) ) );
			if ( ignoreBig == false ){
				this->blitzFFTbig.resize( 2 * this->size[0], 2 * ( 2 * this->size[1] / 2 + 1 ) );
				this->blitzFFTrealbig.reference( this->blitzFFTbig( Range :: all(), Range(fromStart, 2 * this->size[1] - 1 ) ) );
			}
		}
		if ( Dim == 3 ){
			this->blitzFFT.resize( this->size[0], this->size[1], 2 * ( this->size[2] / 2 + 1 ) );
			this->blitzFFTreal.reference( this->blitzFFT( Range :: all(), Range :: all(), Range(fromStart, this->size[2] - 1 ) ) );
			if ( ignoreBig == false ){
				this->blitzFFTbig.resize( 2 * this->size[0], 2 * this->size[1], 2 * ( 2 * this->size[2] / 2 + 1 ) );
				this->blitzFFTrealbig.reference( this->blitzFFTbig( Range :: all(), Range :: all(), Range(fromStart, 2 * this->size[2] - 1 ) ) );
			}
		}

		this->shift.resize( this->size );
		if ( ignoreBig == false ){
			this->shiftbig.resize( 2 * this->size );
		}

		typename Array< Pixel, Dim > :: iterator iter = this->shift.begin();
		while ( iter != this->shift.end() ){
			(*iter) = pow( -1.0, (int)sum( iter.position() ) );
			++iter;
		}

		if ( ignoreBig == false ){
			iter = this->shiftbig.begin();
			while ( iter != this->shiftbig.end() ){
				(*iter) = pow( -1.0, (int)sum( iter.position() ) );
				++iter;
			}
		}

		this->FFTreal = this->blitzFFT.data();

		if ( ignoreBig == false ){
			this->FFTrealbig  = this->blitzFFTbig.data();
		}
		
		if ( this->fftwNotInitializedHalf == false ){
			fftw_destroy_plan( this->fftplanreal );
			fftw_destroy_plan( this->ifftplanreal );
			if ( this->ifftplanrealbig != NULL ){
				fftw_destroy_plan( this->ifftplanrealbig );
			}
		}

		if ( Dim == 2 ){
			this->fftplanreal = fftw_plan_dft_r2c_2d( size[0], size[1], this->FFTreal, reinterpret_cast<fftw_complex*>(this->FFTreal), FFTW_MEASURE );
			this->ifftplanreal = fftw_plan_dft_c2r_2d( size[0], size[1], reinterpret_cast<fftw_complex*>(this->FFTreal), this->FFTreal, FFTW_MEASURE );
			if ( ignoreBig == false ){
				this->ifftplanrealbig = fftw_plan_dft_c2r_2d( 2 * size[0], 2 * size[1], reinterpret_cast<fftw_complex*>(this->FFTrealbig), this->FFTrealbig, FFTW_MEASURE );
			} else {
				this->ifftplanrealbig = NULL;
			}
		}
		if ( Dim == 3 ){
			this->fftplanreal = fftw_plan_dft_r2c_3d( size[0], size[1], size[2], this->FFTreal, reinterpret_cast<fftw_complex*>(this->FFTreal), FFTW_MEASURE );
			this->ifftplanreal = fftw_plan_dft_c2r_3d( size[0], size[1], size[2], reinterpret_cast<fftw_complex*>(this->FFTreal), this->FFTreal, FFTW_MEASURE );
			if ( ignoreBig == false ){
				this->ifftplanrealbig = fftw_plan_dft_c2r_3d( 2 * size[0], 2 * size[1], 2 * size[2], reinterpret_cast<fftw_complex*>(this->FFTrealbig), this->FFTrealbig, FFTW_MEASURE );
			} else {
				this->ifftplanrealbig = NULL;
			}
		}

		this->fftwNotInitialized = true;
		this->fftwNotInitializedHalf = false;

	}
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: centerFFT( vtkImageData * data )
{
	//vtkImageFourierCenter * centerfft = vtkImageFourierCenter::New();
	//centerfft->SetInput( data );
	//centerfft->Update();
	//data->DeepCopy( centerfft->GetOutput() );
	//centerfft->Delete();
	//return;


	if (this->centerfft == NULL) {
		this->centerfft = vtkImageFourierCenter::New();
	}

	this->centerfft->SetInput( data );
	this->centerfft->Update();
	data->DeepCopy( this->centerfft->GetOutput() );
}

#if 0
template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: centerFFT( Array< double, Dim > & A )
{
	//Array< Pixel, 3 > B( A.shape() );
	//B = A;
	vtkImageData * data = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( A, data );
	this->centerFFT( data );
	nbfVTKInterface::vtkToBlitz( data, A );
	//A = B;
	data->Delete();
}
#else
template< class Pixel, const int Dim >
template< class otherPixel >
void nbfFourierFilter< Pixel, Dim > :: centerFFT( Array< otherPixel, Dim > & A )
{
	// assume size(A) is even
	TinyVector< int, Dim > middle = floor( A.shape() / 2 );
	Array< otherPixel, Dim > S( middle );

	Range lowerI( 0, S.rows() - 1 ), upperI( S.rows(), A.ubound(firstDim) );
	Range lowerJ( 0, S.cols() - 1 ), upperJ( S.cols(), A.ubound(secondDim) );

	if ( Dim == 2 ){
		S = A( lowerI, lowerJ );
		A( lowerI, lowerJ ) = A( upperI, upperJ );
		A( upperI, upperJ ) = S;

		S = A( lowerI, upperJ );
		A( lowerI, upperJ ) = A( upperI, lowerJ );
		A( upperI, lowerJ ) = S;
	}
	if ( Dim == 3 ){
		Range lowerK( 0, S.depth() - 1 ), upperK( S.depth(), A.ubound(thirdDim) );
		
		S = A( lowerI, lowerJ, lowerK );
		A( lowerI, lowerJ, lowerK ) = A( upperI, upperJ, upperK );
		A( upperI, upperJ, upperK ) = S;

		S = A( lowerI, lowerJ, upperK );
		A( lowerI, lowerJ, upperK ) = A( upperI, upperJ, lowerK );
		A( upperI, upperJ, lowerK ) = S;

		S = A( lowerI, upperJ, lowerK );
		A( lowerI, upperJ, lowerK ) = A( upperI, lowerJ, upperK );
		A( upperI, lowerJ, upperK ) = S;

		S = A( lowerI, upperJ, upperK );
		A( lowerI, upperJ, upperK ) = A( upperI, lowerJ, lowerK );
		A( upperI, lowerJ, lowerK ) = S;
	}
}
#endif

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: updateFilter()
{
	if ( this->filterUpdated == false ){
		
		/// initialization
		this->filter.resize( this->size );
		this->center = floor( this->size / 2.0 );
		//this->center = floor( this->size / 2.0 ) - 1;

		this->filter = 1;
		//if ( this->wedge.size() > 0 ){
		//	this->filter = wedge;
		//	this->dummyFilter = false;
		//}
		//else{
		//	this->filter = 1;
		//}

		if ( this->useLowPass == true ){
			this->buildLowPass();
			this->filter *= this->lowPass;
			this->dummyFilter = false;
		}
		if ( this->useHighPass == true ){
			this->buildHighPass();
			this->filter *= this->highPass;
			this->dummyFilter = false;
		}
		if ( this->useBandPass == true ){
			this->buildBandPass();
			this->filter *= this->bandPass;
			this->dummyFilter = false;
		}
		if ( this->useGaussian == true ){
			this->buildGaussian();
			this->filter *= this->gaussian;
			this->dummyFilter = false;
		}
		if ( this->bfactorExponent != 0 ){
			this->buildBFactor();
			this->filter *= this->bfactor;
			this->dummyFilter = false;
		}

		//nbfMatlabWriter w;
		//w.setFileName("p.matlab");
		//w.write(this->filter);
		//w.write(this->wedge);

		this->filterUpdated = true;
		this->filterUpdatedHalf = false;
	}
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: updateFilterHalf()
{
	if ( this->filterUpdatedHalf == false ){
		
		/// initialization
		if ( Dim == 2 ){
			this->filter.resize( this->size[0], this->size[1] / 2 + 1 );
		}
		if ( Dim == 3 ){
			this->filter.resize( this->size[0], this->size[1], this->size[2] / 2 + 1 );
		}
		this->center = floor( this->size / 2.0 );

		this->filter = 1;
		if ( Dim == 2 ){
			this->filter( this->filter.rows() / 2, this->filter.cols() - 1 ) = 0;
		}
		if ( Dim == 3 ){
			this->filter( this->filter.rows() / 2, this->filter.cols() / 2, this->filter.depth() - 1 ) = 0;
		}

		if ( this->useLowPass == true ){
			this->buildLowPassHalf();
			this->filter *= this->lowPass;
			this->dummyFilter = false;
		}
		if ( this->useHighPass == true ){
			this->buildHighPassHalf();
			this->filter *= this->highPass;
			this->dummyFilter = false;
		}
		if ( this->useBandPass == true ){
			this->buildBandPassHalf();
			this->filter *= this->bandPass;
			this->dummyFilter = false;
		}
		if ( this->useGaussian == true ){
			this->buildGaussianHalf();
			this->filter *= this->gaussian;
			this->dummyFilter = false;
		}
		if ( this->bfactorExponent != 0 ){
			this->buildBFactorHalf();
			this->filter *= this->bfactor;
			this->dummyFilter = false;
		}

		//nbfMatlabWriter w;
		//w.setFileName("p.matlab");
		//w.write(this->filter);
		//w.write(this->wedge);

		this->filterUpdated = false;
		this->filterUpdatedHalf = true;
	}
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: execute( Array< Pixel, Dim > & real, Array< Pixel, Dim > & imag )
{
	this->updateFilter();

	// apply filter
	if ( this->dummyFilter == false ){
		if ( this->wedge.size() > 0 ){
			real *= ( this->filter * this->wedge );
			imag *= ( this->filter * this->wedge );
		} else {
			real *= this->filter;
			imag *= this->filter;
		}
	}
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: execute( Array< complex< double >, Dim > & complex )
{
	this->updateFilter();

	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");
	//w.write( this->filter );
	//w.write( real(complex) );
	//w.write( imag(complex) );
	// apply filter
	if ( this->dummyFilter == false ){
		if ( this->wedge.size() > 0 ){
			complex *= ( this->filter * this->wedge );
		} else {
			complex *= this->filter;
		}
	}
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: executeHalf( Array< complex< double >, Dim > & complex, bool mutual )
{
	this->updateFilterHalf();

	if ( mutual == true ){
		// take sqrt of modulo and apply filter
		Array< double, Dim > C( complex.shape() );
		C = pow( real( complex * conj( complex ) ), .25 );
		complex = where( C > 0, complex / C, complex );
	}

	// apply filter
	if ( this->dummyFilter == false ){
		if ( this->wedge.size() > 0 ){
			complex *= ( this->filter * this->wedge );
		} else {
			complex *= this->filter;
		}
	}
}


template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: execute( vtkImageData * re )
{
	Array< double, Dim > A;
	nbfVTKInterface :: vtkToBlitzReference( re, A );
	this->execute(A);
	//nbfVTKInterface :: blitzToVtk( A, re );
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: execute( Array< double, Dim > & A, int code )
{
	this->initializeFFTW( A.shape() );

	real( this->blitzFFT1 ) = A * this->shift;
	imag( this->blitzFFT1 ) = 0;

	fftw_execute(this->fftplan);

	switch ( code ){
	// apply filter as usual
	case 0:
		// apply filter and take sqrt of modulo
		this->execute(this->blitzFFT2);
		break;
	case 1:
		this->execute(this->blitzFFT2);
		real( this->blitzFFT1 ) = real( this->blitzFFT2 * conj( this->blitzFFT2 ) );
		this->blitzFFT2 = where( real( this->blitzFFT1 ) > 0, this->blitzFFT2 / pow( real( this->blitzFFT1 ), .25 ), this->blitzFFT2 );
		break;
		// do not apply filter and take sqrt of modulo
	case 2:
		real( this->blitzFFT1 ) = real( this->blitzFFT2 * conj( this->blitzFFT2 ) );
		this->blitzFFT2 = where( real( this->blitzFFT1 ) > 0, this->blitzFFT2 / pow( real( this->blitzFFT1 ), .25 ), this->blitzFFT2 );
		break;
	default:
		cerr << "ERROR: Fourier filter is undefined. In " __FILE__ << "," << __LINE__ << endl;
	}
	
	fftw_execute(this->ifftplan);
	
	A = real(this->blitzFFT1) * this->shift / A.size();
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: executeHalf( Array< double, Dim > & A, int code )
{
	this->initializeFFTWhalf( A.shape() );

	this->blitzFFTreal = A * this->shift;
	
	fftw_execute(this->fftplanreal);

	TinyVector< int, Dim > shape;
	if ( Dim == 2 ){
		shape = TinyVector< int, Dim >( A.rows(), A.cols() / 2 + 1 );
	}
	if ( Dim == 3 ){
		shape = TinyVector< int, Dim >( A.rows(), A.cols(), A.depth() / 2 + 1 );
	}
	Array< complex< double >, Dim > B( reinterpret_cast< complex<double>* >( this->blitzFFT.data() ), shape, neverDeleteData );

	this->executeHalf( B, code );
	//Array< double, Dim > C( B.shape() );
	//switch ( code ){
	//case 0:
	//	// apply filter as usual
	//	this->executeHalf(B);
	//	break;
	//case 1:
	//	// take sqrt of modulo and apply filter
	//	C = pow( real( B * conj( B ) ), .25 );
	//	B = where( C > 0, B / C, B );
	//	this->executeHalf(B);
	//	break;
	//case 2:
	//	// do not apply filter and take sqrt of modulo
	//	C = pow( real( B * conj( B ) ), .25 );
	//	B = where( C > 0, B / C, B );
	//	break;
	//default:
	//	cerr << "ERROR: Fourier filter is undefined. In " __FILE__ << "," << __LINE__ << endl;
	//}
	
	fftw_execute(this->ifftplanreal);
	
	A = this->blitzFFTreal * this->shift / A.size();
}


template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: bandPassOn( Pixel lowCut, Pixel highCut, Pixel lowVariance, Pixel highVariance )
{
	if ( ( lowCut > 0.0 ) || ( highCut < 1.0 ) ){
		this->bandLowCut = lowCut;
		this->bandHighCut = highCut;
		this->bandLowVariance = lowVariance;
		this->bandHighVariance = highVariance;

		this->useBandPass = true;
		this->filterUpdated = false;
		this->filterUpdatedHalf = false;
	}
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: buildBandPass(){

	this->bandPass.resize( this->size );

	Array< Pixel, Dim > radius( this->size );

	//firstIndex i; secondIndex j; thirdIndex k;
	//radius = sqrt( pow2( ( i - this->center[0] ) / (double)this->center[0] ) + 
	//	           pow2( ( j - this->center[1] ) / (double)this->center[1] ) +
	//			   pow2( ( k - this->center[2] ) / (double)this->center[2] ) + 0.0 );

	typename Array< Pixel, Dim > :: iterator iter = radius.begin();
	while ( iter != radius.end() ){
		TinyVector< Pixel, Dim > tmp( ( 1.0 * iter.position() - this->center ) / this->center );
		(*iter) = sqrt( sum( tmp * tmp ) );
		++iter;
	}

	this->bandPass = where( radius < this->bandLowCut, exp( - pow2(this->bandLowCut-radius) / 2.0 / this->bandLowVariance / this->bandLowVariance ), 1 );
	this->bandPass = where( radius > this->bandHighCut, exp( - pow2(radius-this->bandHighCut) / 2.0 / this->bandHighVariance / this->bandHighVariance ), this->bandPass );
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: buildBandPassHalf()
{
	if ( Dim == 2 ){
        this->bandPass.resize( this->size[0], this->size[1] / 2 + 1 );
	}
	if ( Dim == 3 ){
        this->bandPass.resize( this->size[0], this->size[1], this->size[2] / 2 + 1 );
	}

	Array< Pixel, Dim > radius( this->bandPass.shape() );

	typename Array< Pixel, Dim > :: iterator iter = radius.begin();
	while ( iter != radius.end() ){
		TinyVector< Pixel, Dim > tmp( ( 1.0 * iter.position() - this->center ) / this->center );
		(*iter) = sqrt( sum( tmp * tmp ) );
		++iter;
	}

	this->bandPass = where( radius < this->bandLowCut, exp( - pow2(this->bandLowCut-radius) / 2.0 / this->bandLowVariance / this->bandLowVariance ), 1 );
	this->bandPass = where( radius > this->bandHighCut, exp( - pow2(radius-this->bandHighCut) / 2.0 / this->bandHighVariance / this->bandHighVariance ), this->bandPass );
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: buildBFactor(){

	this->bfactor.resize( this->size );

	Array< Pixel, Dim > radius( this->size );

	typename Array< Pixel, Dim > :: iterator iter = radius.begin();
	while ( iter != radius.end() ){
		TinyVector< Pixel, Dim > tmp( ( 1.0 * iter.position() - this->center ) / this->center );
		(*iter) = sqrt( sum( tmp * tmp ) );
		++iter;
	}

	this->bfactor = exp( - this->bfactorExponent * pow2( radius ) );
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: buildBFactorHalf()
{
	if ( Dim == 2 ){
        this->bfactor.resize( this->size[0], this->size[1] / 2 + 1 );
	}
	if ( Dim == 3 ){
        this->bfactor.resize( this->size[0], this->size[1], this->size[2] / 2 + 1 );
	}

	Array< Pixel, Dim > radius( this->bfactor.shape() );

	typename Array< Pixel, Dim > :: iterator iter = radius.begin();
	while ( iter != radius.end() ){
		TinyVector< Pixel, Dim > tmp( ( 1.0 * iter.position() - this->center ) / this->center );
		(*iter) = sqrt( sum( tmp * tmp ) );
		++iter;
	}

	this->bfactor = exp( - this->bfactorExponent * pow2( radius ) );
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: gaussianOn( Pixel gaussianMean, Pixel gaussianVariance )
{
	this->gaussianMean = gaussianMean;
	this->gaussianVariance = gaussianVariance;

	this->useGaussian = true;
	this->filterUpdated = false;
	this->filterUpdatedHalf = false;
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: buildGaussian(){
	this->gaussian.resize( this->size );

	Array< Pixel, Dim > radius( this->size );

	typename Array< Pixel, Dim > :: iterator iter = radius.begin();
	while ( iter != radius.end() ){
		TinyVector< Pixel, Dim > tmp( ( 1.0 * iter.position() - this->center ) / this->center );
		(*iter) = sqrt( sum( tmp * tmp ) );
		++iter;
	}
	this->gaussian = exp( - pow2( radius - this->gaussianMean ) / 2.0 / this->gaussianVariance / this->gaussianVariance );
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: lowPassOn( TinyVector< Pixel, Dim > & cut, int order )
{
	this->lowOrder = order;
	this->lowCut = cut;

	this->useLowPass = true;
	this->filterUpdated = false;
	this->filterUpdatedHalf = false;
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: buildLowPass(){

	this->lowPass.resize( this->size );

	//firstIndex i; secondIndex j; thirdIndex k;

	//this->lowPass = 1.0 / ( 1.0 + pow( pow2( ( i - this->center[0] ) / 2.0 / (double)(this->center[0]) / this->lowCutX ) + 
	//	                               pow2( ( j - this->center[1] ) / 2.0 / (double)(this->center[1]) / this->lowCutY ) + 
	//					               pow2( ( k - this->center[2] ) / 2.0 / (double)(this->center[2]) / this->lowCutZ ),
	//						           this->lowOrder ) );

	typename Array< Pixel, Dim > :: iterator iter = this->lowPass.begin();
	while ( iter != this->lowPass.end() ){
		TinyVector< Pixel, Dim > tmp( ( 1.0 * iter.position() - this->center ) / 2.0 / this->center / this->lowCut );
		(*iter) = 1.0 / ( 1.0 + pow( (Pixel)sum( tmp * tmp ), this->lowOrder ) );
		++iter;
	}
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: highPassOn( TinyVector< Pixel, Dim > & cut, int order )
{
	this->highOrder = order;
	this->highCut = cut;

	this->useHighPass = true;
	this->filterUpdated = false;
	this->filterUpdatedHalf = false;
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: buildHighPass(){

	this->highPass.resize( this->size );

	//firstIndex i; secondIndex j; thirdIndex k;

	//this->highPass = 1.0 / ( 1.0 + 1.0 / pow( pow2( ( i - this->center[0] ) / 2.0 / (double)(this->center[0]) / this->highCutX ) + 
	//	                                 pow2( ( j - this->center[1] ) / 2.0 / (double)(this->center[1]) / this->highCutY ) + 
	//					                 pow2( ( k - this->center[2] ) / 2.0 / (double)(this->center[2]) / this->highCutZ ),
	//						             this->highOrder ) );

	typename Array< Pixel, Dim > :: iterator iter = this->highPass.begin();
	while ( iter != this->highPass.end() ){
		TinyVector< Pixel, Dim > tmp( ( 1.0 * iter.position() - this->center ) / 2.0 / this->center / this->highCut );
		(*iter) = 1.0 / ( 1.0 + 1.0 / pow( (Pixel)sum( tmp * tmp ), this->lowOrder ) );
		++iter;
	}
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: fft( Array< complex< double >, Dim > & in, Array< complex< double >, Dim > & out )
{
	this->initializeFFTW( in.shape() );
	this->blitzFFT1 = in;
	fftw_execute(fftplan);
	out = this->blitzFFT2;
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: fft( Array< Pixel, Dim > & in, Array< complex< double >, Dim > & out )
{
	this->initializeFFTW( in.shape() );
	real(this->blitzFFT1) = in;
	imag(this->blitzFFT1) = 0;
	fftw_execute(fftplan);
	out = this->blitzFFT2;
}

template< class Pixel, const int Dim >
void nbfFourierFilter< Pixel, Dim > :: ifft( Array< complex< double >, Dim > & in, Array< complex< double >, Dim > & out )
{
	this->initializeFFTW( in.shape() );
	this->blitzFFT1 = in;
	fftw_execute(ifftplan);
	out = this->blitzFFT2;
}
