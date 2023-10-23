#pragma once

#include <io/nbfVTKInterface.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfFourierImageMetric.h>
#include <em/nbfFourierFilter.h>
#include <fm/nbfFastMarching.h>
#include <bs/nbfBordStrategyConst.h>
#include <bs/nbfBordStrategyMirror.h>
#include <nbfPolarDomain.h>
#include <nbfCylindricalDomain3.h>
#include <nbfMorphologyFilters.h>

#include <vtkSphericalTransform.h>
#include <vtkCylindricalTransform.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageContinuousDilate3D.h>
#include <vtkImageConstantPad.h>
#include <vtkImageDilateErode3D.h>
#include <vtkStructuredPointsWriter.h>

#include <nbfTimer.h>
// #define REAL double_debug

using namespace blitz;

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <fftw3.h>
#include "complex.h" 
#include "csecond.h"

extern "C"
{
#include "so3_correlate_fftw.h"
#include "soft_fftw.h"
#include "FST_semi_memo.h"
#include "cospmls.h"
#include "primitive_FST.h"
#include "seminaive.h"
#include "utils_so3.h"
#include "rotate_so3.h"
}

//BZ_DECLARE_STENCIL2(localMaxima3D3x3,A,B)
//	if ( ( A(0,0,0) >=   A(0,0,1) ) && ( A(0,0,0) >=   A(0,0,-1) ) &&
//		 ( A(0,0,0) >=   A(0,1,1) ) && ( A(0,0,0) >=   A(0,1,-1) ) && ( A(0,0,0) >=  A(0,1,0) ) &&
//		 ( A(0,0,0) >=  A(0,-1,1) ) && ( A(0,0,0) >=  A(0,-1,-1) ) && ( A(0,0,0) >= A(0,-1,0) ) &&
//		 ( A(0,0,0) >=   A(1,0,1) ) && ( A(0,0,0) >=   A(1,0,-1) ) && ( A(0,0,0) >=  A(1,0,0) ) &&
//		 ( A(0,0,0) >=   A(1,1,1) ) && ( A(0,0,0) >=   A(1,1,-1) ) && ( A(0,0,0) >=  A(1,1,0) ) &&
//		 ( A(0,0,0) >=  A(1,-1,1) ) && ( A(0,0,0) >=  A(1,-1,-1) ) && ( A(0,0,0) >= A(1,-1,0) ) &&
//		 ( A(0,0,0) >=  A(-1,0,1) ) && ( A(0,0,0) >=  A(-1,0,-1) ) && ( A(0,0,0) >= A(-1,0,0) ) &&
//		 ( A(0,0,0) >=  A(-1,1,1) ) && ( A(0,0,0) >=  A(-1,1,-1) ) && ( A(0,0,0) >= A(-1,1,0) ) &&
//		 ( A(0,0,0) >= A(-1,-1,1) ) && ( A(0,0,0) >= A(-1,-1,-1) ) && ( A(0,0,0) >= A(-1,-1,0) ) ){
//		B(0,0,0) = 1;
//	//if ( ( A(0,0,0) >   A(0,0,1) ) && ( A(0,0,0) >   A(0,0,-1) ) &&
//	//	 ( A(0,0,0) >   A(0,1,0) ) && ( A(0,0,0) >   A(0,-1,0) ) &&
//	//	 ( A(0,0,0) >  A(1,0,0) ) && ( A(0,0,0) >  A(-1,0,0) ) ){
//	//	B = A(0,0,0);
//	} else {
//		B(0,0,0) = 0;
//	}
//BZ_STENCIL_END

BZ_DECLARE_STENCIL_OPERATOR1(localMaxima3D3x3, A)
    return ( ( A(0,0,0) >= A(0,0,1)   ) && ( A(0,0,0) >= A(0,0,-1)   ) &&
		     ( A(0,0,0) >= A(0,1,1)   ) && ( A(0,0,0) >= A(0,1,-1)   ) && ( A(0,0,0) >= A(0,1,0)   ) &&
		     ( A(0,0,0) >= A(0,-1,1)  ) && ( A(0,0,0) >= A(0,-1,-1)  ) && ( A(0,0,0) >= A(0,-1,0)  ) &&
		     ( A(0,0,0) >= A(1,0,1)   ) && ( A(0,0,0) >= A(1,0,-1)   ) && ( A(0,0,0) >= A(1,0,0)   ) &&
		     ( A(0,0,0) >= A(1,1,1)   ) && ( A(0,0,0) >= A(1,1,-1)   ) && ( A(0,0,0) >= A(1,1,0)   ) &&
		     ( A(0,0,0) >= A(1,-1,1)  ) && ( A(0,0,0) >= A(1,-1,-1)  ) && ( A(0,0,0) >= A(1,-1,0)  ) &&
		     ( A(0,0,0) >= A(-1,0,1)  ) && ( A(0,0,0) >= A(-1,0,-1)  ) && ( A(0,0,0) >= A(-1,0,0)  ) &&
		     ( A(0,0,0) >= A(-1,1,1)  ) && ( A(0,0,0) >= A(-1,1,-1)  ) && ( A(0,0,0) >= A(-1,1,0)  ) &&
		     ( A(0,0,0) >= A(-1,-1,1) ) && ( A(0,0,0) >= A(-1,-1,-1) ) && ( A(0,0,0) >= A(-1,-1,0) ) ) * A(0,0,0);
BZ_END_STENCIL_OPERATOR

BZ_ET_STENCIL(localMaxima3D3x3,P_numtype)

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
class nbfProjectionRotationMetric3D : public nbfCorrelationImageMetric< Pixel, 3 >
{
public:

	/// Constructor.
	nbfProjectionRotationMetric3D( nbfImageFilter< Pixel, 3 > * = NULL, nbfFourierFilter< Pixel, 3 > * = NULL );

	virtual ~nbfProjectionRotationMetric3D();

	virtual int getId(){ return NBF_IMAGE_METRIC_PROJECTION; }

	/// Redefine from nbfImageMetric
	void execute();

	void soft( Array< REAL, 2 > &, Array< REAL, 2 > &, Array< REAL, 3 > &, bool = false );
	// same as soft but with proper scaling
	void softScaled( Array< REAL, 2 > &, Array< REAL, 2 > &, Array< REAL, 3 > &, bool = false );

	void initializeHarmonics( int,int );
	void harmonics( Array< REAL, 2 > &, REAL, REAL, REAL, Array< REAL, 2 > & );
	void harmonics( Array< Pixel, 2 > & A, vtkTransform * t, Array< Pixel, 2 > & R );

	// Given a correlation volume in SO3, produce rotated images on the sphere corresponding to the maxima and minima of correlation
	void rotatedSphericalImageMinMax( Array< REAL, 3 > &, Array< REAL, 2 > &, Array< REAL, 2 > &, Array< REAL, 2 > &, REAL &, REAL & );

	void cartesianToSpherical( vtkImageData *, Array< Pixel, 2 > & );

	// Instead of selecting the single best peak of spherical CCC volume, 
	// set to try given number of best peaks and keep the best one. Default = 10;
	void setNumberOfCandidatePeaksToSearch( int a ){ this->candidatePeaks = a; }
	int getNumberOfCandidatePeaksToSearch(){ return this->candidatePeaks; }
	
	// Set maximum translation allowance. If translation is biger that value, distance is set to 1.
	// Default = 20.
	void setTranslationAllowance( Pixel a ){ this->translationAllowance = a; }

	///// Redefine from nbfImageMetric
	//Pixel getWedgeOverlap();

	/// Distance value between inputs (function of the correlation value). Redefine from parent.
	virtual Pixel getDistance(){ Pixel p = nbfImageMetric< Pixel, 3 > :: getDistance(); return this->correlationPeak; }

	// spherical bandwidth
	float B;

	void conditionSphericalWindow( Array< Pixel, 2 > & );

protected:

	void getTrueLocalMaximaPositions( Array< REAL, 3 > &,
		PriorityQueue< TinyVector< Pixel, 4 >, greaterPD< Pixel, 3 > , vector< TinyVector< Pixel, 4 > > > &,
		Pixel = 0.25 ); // this error is the norm between the transformation matrices

	void eliminateEquivalentAngles( vector< TinyVector< Pixel, 4 > > &,
									PriorityQueue< TinyVector< Pixel, 4 >, greaterPD< Pixel, 3 > , vector< TinyVector< Pixel, 4 > > > &,
									Pixel = 0.25 );

	int candidatePeaks;
	Pixel translationAllowance;

	//nbfPolarDomain< Pixel, 2 > polar;

	//Array< Pixel, 2 > proj1, proj2;

	//Array< Pixel, 2 > polar1, polar2;
	//Array< bool, 2 > B;

	//nbfFourierFilter< Pixel, 2 > filter2D;

	//nbfCorrelationImageMetric< Pixel, 2 > ccc;

	Array< Pixel, 2 > Pp1saved;

	vtkImageConstantPad * pad;
	vtkImageChangeInformation * change;
	vtkSphericalTransform * spherical;
	vtkImageReslice * reslice;
	vtkImageData * magnitude;

private:

	// spherical harmonics
	bool harmonicsInitialized;
	int bwIn, bwOut, degLim ;
	fftw_complex *workspace1, *workspace2  ;
	REAL *workspace3 ;
	REAL *sigR, *sigI ;
	REAL *sigCoefR, *sigCoefI ;
	REAL *patCoefR, *patCoefI ;
	fftw_complex *so3Sig, *so3Coef ;
	fftw_plan p1 ;
	REAL *seminaive_naive_tablespace  ;
	REAL **seminaive_naive_table ;
	int na[2], inembed[2], onembed[2] ;
	int rank, howmany, istride, idist, ostride, odist ;

	// rotation on the unit sphere
	bool rotation_spherical_harmonics_initialized;
	double *rot_seminaive_naive_tablespace, *rot_trans_seminaive_naive_tablespace2;
	double *rot_seminaive_naive_tablespace2 ;
	double **rot_seminaive_naive_table2,**rot_seminaive_naive_table ;
	double **rot_trans_seminaive_naive_table2;
	REAL *rot_sigInR, *rot_sigInI, *rot_sigOutR, *rot_sigOutI ;
	REAL *rot_scratch ;

};

template< class Pixel >
nbfProjectionRotationMetric3D< Pixel > :: nbfProjectionRotationMetric3D( nbfImageFilter< Pixel, 3 > * imf, nbfFourierFilter< Pixel, 3 > * fff )
: nbfCorrelationImageMetric< Pixel, 3 >( imf, fff )
{
	// this->filter.bandPassOn(.02,.5,.01);
	// this->ccc.windowPolarOn();
	//this->ccc.windowOff();
	this->harmonicsInitialized = false;
	this->rotation_spherical_harmonics_initialized = false;
	
	this->B = 64;
	this->candidatePeaks = 5;
	this->translationAllowance = 20;

	this->pad = vtkImageConstantPad::New();
	this->change = vtkImageChangeInformation::New();
	this->change->CenterImageOn();
	this->spherical = vtkSphericalTransform::New();
	this->reslice = vtkImageReslice::New();
	this->reslice->SetResliceTransform( this->spherical );
	this->reslice->SetInterpolationModeToCubic();
	this->reslice->WrapOn();
	this->magnitude = vtkImageData::New();
}

template< class Pixel >
nbfProjectionRotationMetric3D< Pixel > :: ~nbfProjectionRotationMetric3D()
{
	if ( this->harmonicsInitialized == true ){
		free( this->seminaive_naive_table ) ;
		free( this->seminaive_naive_tablespace ) ;
		fftw_free( this->so3Coef ) ;
		free( this->patCoefI );
		free( this->patCoefR );
		free( this->sigCoefI );
		free( this->sigCoefR );
		free( this->workspace3 );
		fftw_free( this->workspace2 );
		fftw_free( this->workspace1 );
		fftw_free( this->so3Sig ) ;
		free( this->sigI );
	}

	this->pad->Delete();
	this->change->Delete();
	this->spherical->Delete();
	this->reslice->Delete();
	this->magnitude->Delete();

	if ( this->rotation_spherical_harmonics_initialized == true ){
		free(this->rot_trans_seminaive_naive_table2);
		free(this->rot_seminaive_naive_table2);
		free(this->rot_seminaive_naive_table);

		free(this->rot_seminaive_naive_tablespace2);
		free(this->rot_trans_seminaive_naive_tablespace2);
		free(this->rot_seminaive_naive_tablespace);

		free(this->rot_scratch);
		free(this->rot_sigOutI);
		free(this->rot_sigOutR);
		free(this->rot_sigInI);
		free(this->rot_sigInR);
	}
}

template< class Pixel >
void nbfProjectionRotationMetric3D< Pixel > :: cartesianToSpherical( vtkImageData * input, Array< Pixel, 2 > & Pp )
{
	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	// window data in real space
	this->imageFilter->execute( input );

	// PAD INPUT IMAGES IF NECCESARY
	int padc = 2;
	if ( this->imageFilter->getPaddingFactor() > 1 ){
		padc = 1;
	}

	int extent[6];
	input->GetExtent(extent);
	if ( padc == 1 ){
		this->pad->SetOutputWholeExtent( extent[0], extent[1]+0, extent[2], extent[3]+0, extent[4], extent[5]+0 );
	} else if ( padc == 2 ){
		this->pad->SetOutputWholeExtent( extent[0], padc*extent[1]+1, extent[2], padc*extent[3]+1, extent[4], padc*extent[5]+1 );
	}

	// set pad constant to mean of un-windowed image (the windowed image has this value in its boundary)
	this->pad->SetConstant( input->GetScalarComponentAsDouble(0,0,0,0) );
	this->pad->SetInput( input );
	this->pad->Update();

	Array< double, 3 > A;
	if ( input->GetScalarType() == VTK_DOUBLE ){
		nbfVTKInterface::vtkToBlitzReference( this->pad->GetOutput(), A );
	} else {
		nbfVTKInterface::vtkToBlitz( this->pad->GetOutput(), A );
	}

	//w.write(A);

#if 0
	this->fourierFilter->initializeFFTW( TinyVector< int, 3 >( A.shape() ) );
	real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
	imag( this->fourierFilter->blitzFFT1 ) = 0;
	fftw_execute( this->fourierFilter->fftplan );
	// this->fourierFilter->blitzFFT2 /= this->fourierFilter->blitzFFT2.size();
	// this->fourierFilter->updateFilter();
	// Use mutual correlation function
	if ( this->useMutualCorrelation == true ){
		A = sqrt( sqrt( real( ( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ) ) ) );
	} else {
		A = sqrt( real( ( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ) ) );
	}
	//A = log( 1 + real( ( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ) ) );
#else
	this->fourierFilter->initializeFFTWhalf( TinyVector< int, 3 >( A.shape() ), true );
	this->fourierFilter->blitzFFTreal = A * this->fourierFilter->shift;
	fftw_execute( this->fourierFilter->fftplanreal );
	TinyVector< int, 3 > shape( A.rows(), A.cols(), A.depth() / 2 + 1 );
	Array< complex< double >, 3 > FFTre( reinterpret_cast<complex<double>*>(this->fourierFilter->blitzFFT.data()), shape, neverDeleteData );
	if ( this->useMutualCorrelation == true ){
		real( FFTre ) = sqrt( sqrt( real( ( FFTre * conj( FFTre ) ) ) ) );
	} else {
		real( FFTre ) = sqrt( real( ( FFTre * conj( FFTre ) ) ) );
	}
	// compose full size FFT
	A( Range :: all(), Range :: all(), Range( fromStart, A.depth() / 2 ) ) = real( FFTre );
	A( 0, Range :: all(), Range :: all() ) = 0;
	int r2 = A.rows() / 2;
	int c2 = A.cols() / 2;
	int d2 = A.depth() / 2;
	for ( int i = - r2 + 1; i < r2; i++ ){
		for ( int j = - c2 + 1; j < c2; j++ ){
			for ( int k = 1; k < d2; k++ ){
				A( r2 + i, c2 + j, d2 + k ) = A( r2 - i, c2 - j, d2 - k );
			}
		}
	}
#endif

	//w.write(A);

#if _DEBUG
	//Array< float, 3 > Bol(A.shape());
	//Bol = cast<float>(A);
	//w.write(A);
	//nbfMrcWriter wr;
	//wr.setFileName("f.mrc");
	//vtkImageData * d = vtkImageData::New();
	//nbfVTKInterface::blitzToVtk(Bol,d);
	//wr.write(d);
#endif

	//TinyVector< int, 3 > cen(B,B,B);
	//A( cen ) = 0;
	
	nbfVTKInterface::blitzToVtk( A, this->magnitude );

	this->change->SetInput( this->magnitude );
	this->change->Update();

	REAL origin[3];
	this->change->GetOutput()->GetOrigin(origin);
	origin[0] -= .5;
	origin[1] -= .5;
	origin[2] -= .5;
	this->change->GetOutput()->SetOrigin(origin);

	int dims[3];
	this->magnitude->GetDimensions(dims);

	this->reslice->SetInput( this->change->GetOutput() );

	int bandL = max( 4, (int)ceil( this->fourierFilter->getBandLowCut() * ( dims[0] - 1 ) ) );
	int bandU = ceil( blitz::extrema::min( .5, this->fourierFilter->getBandHighCut() ) * ( dims[0] - 1 ) );
	//reslice->SetOutputExtent( bandL, bandU, 0, 2*B-1, 0, 2*B-1 );
	//reslice->SetOutputExtent( 0, (dims[0]-1.0)/2.0, 0, 2*B-1, 0, 2*B-1 );
	this->reslice->SetOutputExtent( bandL, bandU, 0, 2*B-1, 0, 2*B-1 );
	
	this->reslice->SetOutputOrigin(0, vtkMath::Pi() / 4.0 / this->B, 0 );
	this->reslice->SetOutputSpacing( 1, 2.0 * vtkMath::Pi() / 4.0 / this->B, vtkMath::Pi() / this->B );
	this->reslice->Update();
	
	if ( this->reslice->GetOutput()->GetScalarType() == VTK_DOUBLE ){
		nbfVTKInterface::vtkToBlitzReference( this->reslice->GetOutput(), A );
	} else {
		nbfVTKInterface::vtkToBlitz( this->reslice->GetOutput(), A );
	}

	//A = exp( A );
#if _DEBUG
	//w.write(A);
#endif

	//for ( int i = 0; i < A.rows(); i++ ){
	//	A( i, Range::all(), Range::all() ) =  A( i, Range::all(), Range::all() ) - min( A( i, Range::all(), Range::all() ) );
	//	A( i, Range::all(), Range::all() ) =  A( i, Range::all(), Range::all() ) / max( A( i, Range::all(), Range::all() ) );
	//}
	//w.write(A);

	firstIndex i;
	secondIndex j;
	thirdIndex k;

	Pp.resize( A.cols(), A.depth() );
	Pp = sum( A(k,i,j),k);
	//Pp = cast< Pixel >( A( 9, Range::all(), Range::all() ) );
	//Pp = cast< Pixel >( A( 10, Range::all(), Range::all() ) );

	//// normalization
	//Pp = Pp - min(Pp);
	//Pp = Pp / max(Pp);
}

template< class Pixel >
void nbfProjectionRotationMetric3D< Pixel > :: execute()
{
#if 0
	// Metric to be used for minimization
	nbfFourierImageMetric< Pixel, 3 > cccMetricDummy( this->imageFilter, this->fourierFilter );
	cccMetricDummy.setInput1( this->wedgedInput1 );
	cccMetricDummy.setInput2( this->wedgedInput2 );
	vtkTransform * dummy = vtkTransform::New();
	cccMetricDummy.executeFourierNew( dummy );
	this->correlationPeak = cccMetricDummy.getCorrelationPeak();
	this->correlationScale = cccMetricDummy.getCorrelationScale();
	this->transform->DeepCopy( dummy );
	dummy->Delete();
	return;
#endif

	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	Array< Pixel, 2 > Pp1, Pp2;

	//nbfTimer timer;
	//timer.start();

	// use shape mask for reference volume if available
	this->imageFilter->fileMaskOn();


	// project |FFT| of first image onto sphere
	if ( this->input1Changed == true ){
		this->cartesianToSpherical( this->input1, Pp1 );
		this->Pp1saved.resize( Pp1.shape() );
		this->Pp1saved = Pp1;
		this->input1Changed = false;
	} else {
		Pp1.resize( this->Pp1saved.shape() );
		Pp1 = this->Pp1saved;
	}


	// switch back to isotropic mask
	this->imageFilter->fileMaskOff();

	// project |FFT| of second image onto sphere
	if ( this->input2Changed == true ){
		int dims[3];
		this->input2->GetDimensions(dims);
		this->cartesianToSpherical( this->input2, Pp2 );
		this->input2Changed = false;
	}


	Pixel M = max( max(Pp1), max(Pp2) );
	Pp1 = Pp1 / M;
	Pp2 = Pp2 / M;

	// wedge images on the sphere
	Array< Pixel, 2 > W1( Pp1.shape() ), W2( Pp2.shape() );
	
	bool useWindow1 = false;
	bool useWindow2 = false;


	if ( this->useMissingWedgeCompensation == true ){
		// get binary wedge images
		this->wedgedInput1->getSphericalWedgeImage( W1, (vtkTransform*)NULL, this );
		this->wedgedInput2->getSphericalWedgeImage( W2, (vtkTransform*)NULL, this );

		// smooth wedge images (make them bandlimited)
		this->conditionSphericalWindow(W1);
		this->conditionSphericalWindow(W2);

		// if missing wedge not big enough
		Pixel wedgeThreshold = .95;
		useWindow1 = sum(W1) <= ( wedgeThreshold * W1.size() );
		useWindow2 = sum(W2) <= ( wedgeThreshold * W2.size() );
	} else {
		W1 = Pp1;
		W2 = Pp2;
	}

#ifdef _DEBUG
	w.write(Pp1);
	w.write(Pp2);
	w.write(W1);
	w.write(W2);
#endif

	// images on the sphere for SOFT
	Array< REAL, 2 > F1( Pp1.shape() ), F2( Pp2.shape() );

	//this->wedgedInput1->getImage( this->input1 );
	//Array< Pixel, 3 > K;
	//nbfVTKInterface::vtkToBlitz( this->input1, K );
	//w.setFileName("i1.matlab");
	//w.write(K);

	//// save spherical wedge image + corresponding wedge image
	//this->wedgedInput2->getImage( this->input2 );
	//nbfVTKInterface::vtkToBlitz( this->input2, K );
	//w.setFileName("i2.matlab");
	//w.write(K);
	//w.setFileName("p.matlab");

	// Pixel a1 = 3.0 * vtkMath::Pi();

	Array< REAL, 3 > C1, C2;


	// if both volumes have seizable windows
	if ( useWindow1 && useWindow2 ){

		F1 = Pp1 * Pp1 * W1;
		F2 = cast<REAL>(W2);
		//F1 = W1 * cos( a1 * Pp1 );
		//F2 = W2 * cos( a1 * Pp2 );
		this->softScaled(F1,F2,C1);

		F1 = cast<REAL>(W1);
		F2 = Pp2 * Pp2 * W2;
		//F1 = W1 * sin( a1 * Pp1 );
		//F2 = W2 * sin( a1 * Pp2 );
		this->softScaled(F1,F2,C2);

		C1 = sqrt( C1 ) * sqrt( C2 );

		F1 = Pp1 * W1;
		F2 = Pp2 * W2;
		this->softScaled(F1,F2,C2);

		C1 = 2.0 * ( 1.0 - C2 / C1 );

		// compute overlap
		F1 = cast<REAL>(W1);
		F2 = cast<REAL>(W2);
		this->softScaled(F1,F2,C2);
		
		// invert to find local minima by local maxima algorithm
		C1 = - C1 / C2;
	} 
	else if ( useWindow1 == true ){	// if only first volume has seizable missing wedge
		
		F1 = cast<REAL>(W1); 
		F2 = Pp2 * Pp2;
		this->softScaled(F1,F2,C2);

		F1 = W1 * Pp1;
		F2 = cast<REAL>(Pp2);
		this->softScaled(F1,F2,C1);

		C1 = C1 / sqrt( C2 );
	}
	else if ( useWindow2 == true ){ // if only second volume has seizable misssing wedge

		F1 = Pp1 * Pp1;
		F2 = cast<REAL>(W2); 
		this->softScaled(F1,F2,C2);

		F1 = cast<REAL>(Pp1);
		F2 = W2 * Pp2;
		this->softScaled(F1,F2,C1);

		C1 = C1 / sqrt( C2 );
	} 
	else { // if neither volume has seizable missing wedge
		F1 = cast< REAL >(Pp1);
		F2 = cast< REAL >(Pp2);
		this->soft(F1,F2,C1);
	}

////#ifdef NBF_VERBOSE
//	Array< float, 3 > C1f( C1.shape() );
//	C1f = cast<float>(C1);
//	vtkImageData * c1vtk = vtkImageData :: New();
//	nbfVTKInterface :: blitzToVtk( C1f, c1vtk );
//	nbfMrcWriter mw;
//	mw.setFileName("cost.mrc");
//	mw.write(c1vtk);
//	c1vtk->Delete();
////#endif

	//Array< Pixel, 3 > subC1 = C1( Range::all(), Range(fromStart, blitz::extrema::min( floor( this->restrictAngularSearch * this->bwOut / 90.0 - 1.0 ), C1.ubound(secondDim) ) ), Range::all() );
	//Array< Pixel, 3 > subC1( C1 );
	//Array< Pixel, 3 > subC3( subC1.shape() );

	// restrict search to valid range
	C1( Range::all(), Range( blitz::extrema::max( 1, blitz::extrema::min( floor( this->restrictRotationSearch * this->bwOut / 90.0 - 1.0 ), C1.ubound(secondDim) ) ), toEnd ), Range::all() ) = - numeric_limits< Pixel > :: max();

	// bord strategy to apply local maxima search stencil
	BordStrategyConst< REAL, 3 > bsForC1( C1, 2, - numeric_limits< Pixel > :: max() );
	bsForC1.refresh();

	// compute local maxima in 3x3x3 neighborhood
	C2.resize( C1.shape() );
	C2 = localMaxima3D3x3(C1);

	//// add global maxima in restricted range since this point may not be a local maxima.
	//TinyVector< int, 3 > lm = maxIndex( C1( Range::all(), Range(fromStart, blitz::extrema::min( blitz::extrema::max( 1, floor( this->restrictRotationSearch * this->bwOut / 90.0 - 1.0 ) ), C1.ubound(secondDim) ) ), Range::all() ) ); 
	//C2( lm ) = C1( lm );

	// build priority queue with local maxima points
	PriorityQueue< TinyVector< Pixel, 4 >, greaterPD< Pixel, 3 > , vector< TinyVector< Pixel, 4 > > > pqNB1;

	this->getTrueLocalMaximaPositions( C2, pqNB1 );

	// second priority queue with local maxima points
	PriorityQueue< TinyVector< Pixel, 4 >, greaterPD< Pixel, 3 > , vector< TinyVector< Pixel, 4 > > > pqNB2;


	// if using windows compute peaks without wedge compensation in case they provide better alignment candidates
	if ( useWindow1 || useWindow2 ){
		F1 = cast< REAL >(Pp1);
		F2 = cast< REAL >(Pp2);
		this->soft(F1,F2,C2);

		// restrict search to valid range
		C2( Range::all(), Range( blitz::extrema::max( 1, blitz::extrema::min( floor( this->restrictRotationSearch * this->bwOut / 90.0 - 1.0 ), C2.ubound(secondDim) ) ), toEnd ), Range::all() ) = - numeric_limits< Pixel > :: max();

		C1 = C2;

		bsForC1.refresh();

		// compute local maxima in 3x3x3 neighborhood
		C2 = localMaxima3D3x3(C1);

		//// add global maxima in restricted range since this point may not be a local maxima.
		//TinyVector< int, 3 > lm = maxIndex( C1( Range::all(), Range(fromStart, blitz::extrema::min( blitz::extrema::max( 1, floor( this->restrictRotationSearch * this->bwOut / 90.0 - 1.0 ) ), C1.ubound(secondDim) ) ), Range::all() ) ); 
		//C2( lm ) = C1( lm );

		this->getTrueLocalMaximaPositions( C2, pqNB2 );
	}


	// Metric to be used for minimization
	nbfFourierImageMetric< Pixel, 3 > cccMetric( this->imageFilter, this->fourierFilter );
	cccMetric.setTranslationSearchRestriction( this->restrictTranslationSearch );
	cccMetric.setRotationSearchRestriction( this->restrictRotationSearch );
	cccMetric.setToUseMutualCorrelation( this->getToUseMutualCorrelation() );
	cccMetric.setMissingWedgeCompensation( this->getMissingWedgeCompensation() );
	cccMetric.setToComputeOverlapNormalizedDistances( this->getToComputeOverlapNormalizedDistances() );
	
	cccMetric.setInput1( this->wedgedInput1 );
	cccMetric.setInput2( this->wedgedInput2 );

	PriorityQueue< TinyVector< Pixel, 4 >, greaterPD< Pixel, 3 > , vector< TinyVector< Pixel, 4 > > > pqNBcombined;
	
	PriorityQueue< TinyVector< Pixel, 4 >, greaterPD< Pixel, 3 > , vector< TinyVector< Pixel, 4 > > > * pqNBcurrent;

	// set the number of active queues
	int numberOfActiveQueues = 1;
	if ( useWindow1 || useWindow2 ){
		numberOfActiveQueues = 2;
	}

	// store combination of candidate transforms
	vector< TinyVector< Pixel, 4 > > combinedBestRotations;

	if ( pqNB1.size() + pqNB2.size() == 0 ){
		cerr << "WARNING - No valid rotations were identified in the SH computations." << endl;
	}

	// try three dimensional metric at spherical candidates
	for ( int queues = 0; queues < numberOfActiveQueues; queues++ ){
		if ( queues == 0 ){
			pqNBcurrent = &pqNB1;
		} else if ( queues == 1 ){
			pqNBcurrent = &pqNB2;
		}

		int index = 0;

		// search all peaks only for wedge-aware spherical volume
		// search half the number of peaks for the version without wedge
		int peaks = this->candidatePeaks;
		if ( queues == 1 ){
			peaks = blitz::extrema::max( 1, floor( this->candidatePeaks / 3.0 ) );
		}

		while ( ( index < peaks ) && ( pqNBcurrent->size() > 0 ) ){

			// get current peak from priority queue
			TinyVector< Pixel, 4 > p = pqNBcurrent->top();

			if ( p[3] == numeric_limits< Pixel > :: max() ){
				break;
			}

			// retrieve rotation angles
			Pixel alpha = p[0];
			Pixel beta  = p[1];
			Pixel gamma = p[2];

			vtkTransform * t1 = vtkTransform::New();
			t1->RotateZ( - gamma );
			t1->RotateY( - beta );
			t1->RotateZ( - alpha );

			vtkTransform * tdirect = vtkTransform::New();
			tdirect->SetMatrix( t1->GetMatrix() );
			//this->cartesianToSpherical( this->input2, Ppt );
			//w.write(Ppt);
			//this->wedgedInput2->getSphericalWedgeImage( Ppt, t1 );
			//this->conditionSphericalWindow(Ppt);


			//cccMetric.execute(tdirect);
			cccMetric.executeFourierNewHalf(tdirect);
			Pixel currentD = cccMetric.getCorrelationPeak();
			


	//Array< double, 3 > L;
	//this->wedgedInput2->getImage( this->input2, tdirect );
	//nbfVTKInterface :: vtkToBlitzReference( this->input2, L );

	//Array< Pixel, 3 > Lf( L.shape() );
	//Lf = cast<float>(L);
	//vtkImageData * data=vtkImageData::New();
	//nbfVTKInterface :: blitzToVtk( Lf, data );
	//mw.setFileName("i2.mrc");
	//mw.write( data );
	//data->Delete();


			// now check if translation is within allowable limits
			if ( this->isTransformValid( tdirect ) == true ){

				// try mirror image to discard false positive (this ambiguity comes from using the magnitude of the FFT)
				vtkTransform * mirror = vtkTransform::New();
				mirror->Scale(-1,-1,-1);
				mirror->PreMultiply();
				mirror->Concatenate( t1 );

				cccMetric.executeFourierNewHalf(mirror);
				//cccMetric.execute(mirror);
				mirror->Delete();

				// count relevant minima only
				if ( cccMetric.getCorrelationPeak() >= currentD ){
					index++;
				} else {
//					// skipped
//#ifdef NBF_VERBOSE
//					cout << "False peak, skipping." << endl;
//#endif
				}

#ifdef NBF_VERBOSE
				cout << "R=[" << alpha << "," << beta << "," << gamma << "],\tD=" << currentD << ",\tSH=" << p[3] << endl;		
#endif
				p = TinyVector< Pixel, 4 >( alpha, beta, gamma, currentD );

				// pqNBcombined.push(p);
				
				combinedBestRotations.push_back(p);
			}
			pqNBcurrent->pop();
			t1->Delete();
			tdirect->Delete();
		}
	}


	//// now extract candidates in vector (sorted by fitness)

	//while ( pqNBcombined.size() > 0 ){
	//	combinedBestRotations.push_back( pqNBcombined.top() );
	//	pqNBcombined.pop();
	//}
	
	if ( combinedBestRotations.size() == 0 ){
		cerr << "WARNING - No valid translations resulted from the candidate rotations.\n Aligning ";
		if ( this->wedgedInput1->getTypeId() == NBF_WEDGED_SUB_IMAGE_3D ){
			cerr << reinterpret_cast< nbfWedgedSubImage3D< Pixel > * >( this->wedgedInput1 )->getFileName() << " to ";
		} else {
			cerr << " reference with " << sum( where( reinterpret_cast< nbfWedgedAverageImage3D< Pixel > * >( this->wedgedInput1 )->multipleAlignments, 1, 0 ) ) << " elements to ";
		}
		if ( this->wedgedInput2->getTypeId() == NBF_WEDGED_SUB_IMAGE_3D ){
			cerr << reinterpret_cast< nbfWedgedSubImage3D< Pixel > * >( this->wedgedInput2 )->getFileName() << endl;
		} else {
			cerr << " reference with " << sum( where( reinterpret_cast< nbfWedgedAverageImage3D< Pixel > * >( this->wedgedInput2 )->multipleAlignments, 1, 0 ) ) << " elements.";
		}
		//this->transform->Identity();
		//this->correlationPeak = numeric_limits< Pixel > :: max();
		//return;
		TinyVector< Pixel, 4 > p( 0, 0, 0, 1.0 );
		
		combinedBestRotations.push_back(p);
	}


	// eliminate duplicate angles from combinedBestRotations and store valid positions in pqNBcombined
	this->eliminateEquivalentAngles( combinedBestRotations, pqNBcombined, 2.5 );

	vtkTransform * bestTransform = vtkTransform::New();

	int index = 0;
	
	Pixel bestDistance = numeric_limits< Pixel > :: max();

	this->candidateTransforms.clear();
	this->candidateCorrelationPeaks.clear();
	this->candidateCorrelationScales.clear();
	this->candidateWedgeOverlaps.clear();

	// refine best candidates
	while ( ( pqNBcombined.size() > 0 ) && ( index < this->numberOfCandidates ) ){
		TinyVector< Pixel, 4 > p = pqNBcombined.top();

#ifdef NBF_VERBOSE
		cout << "Refining R=[" << p[0] << "," << p[1] << "," << p[2] << "], D=" << p[3] << endl;
#endif
		vtkTransform * t = vtkTransform::New();
		t->RotateZ( - p[2] ); // - gamma
		t->RotateY( - p[1] ); // - beta
		t->RotateZ( - p[0] ); // - alpha
		
		// descent refinement
		if ( this->refinement == true ){
			cccMetric.refine(t);

			int refinementIterations = 0;

			// new final alignment (computationally cheaper than doing new refine and accurate enough)
			for ( int refine = 0; refine < refinementIterations; refine++ ){
				cccMetric.executeFourierNewHalf( t );
			}

		} else {
			cccMetric.refine(t,0);
		}

#ifdef NBF_VERBOSE
		Pixel translation[3];
		t->GetPosition(translation);
		cout << "Final refined T=[" << translation[0] << "," << translation[1] << "," << translation[2] << "], D=" << cccMetric.getDistance() << endl;
#endif

		if ( cccMetric.getDistance() < bestDistance ){
			bestDistance = cccMetric.getDistance();
			//vtkTransform * tmp = vtkTransform::New();
			bestTransform->DeepCopy( t );
			//tmp->Delete();
			//bestTransform->SetInput( t );

			// assign best distance and scale to attributes
			this->correlationPeak = cccMetric.getCorrelationPeak();
			this->correlationScale = cccMetric.getCorrelationScale();
		}

		// copy transformation matrix into blitz array for storage
		Array< double, 1 > currentMatrix(16);

		double matrix[16];
		vtkMatrix4x4::DeepCopy( matrix, t->GetMatrix() );

		for ( int i = 0; i < 16; i++ ){
			currentMatrix(i) = matrix[i];
		}

		this->candidateTransforms.push_back( currentMatrix );

		this->candidateCorrelationPeaks.push_back( cccMetric.getCorrelationPeak() );
		this->candidateCorrelationScales.push_back( cccMetric.getCorrelationScale() );
		this->candidateWedgeOverlaps.push_back( 1.0 );
		// this->candidateWedgeOverlaps.push_back( cccMetric.getWedgeOverlap( t ) );

	//Array< double, 3 > L;
	//this->wedgedInput2->getImage( this->input2, t );
	//nbfVTKInterface :: vtkToBlitzReference( this->input2, L );
	//Array< float, 3 > Lf( L.shape() );
	//vtkImageData * data = vtkImageData :: New();
	//Lf = cast<float>(L);
	//nbfVTKInterface :: blitzToVtk( Lf, data );
	//stringstream name;
	//name << "2BF1_peak_" << index << ".mrc";
	//mw.setFileName( name.str().c_str() );
	//mw.write( data );
	//data->Delete();


		pqNBcombined.pop();
		t->Delete();

		index++;
	}

	// assign final transformation to attribute
	//vtkTransform * tmp = vtkTransform::New();
	//this->transform->DeepCopy( tmp );
	//tmp->Delete();
	this->transform->DeepCopy( bestTransform );

#undef FITTINGHIV
#ifdef FITTINGHIV
	// BUILD BRUTE FORCE VOLUME TO SEE METRIC LANDSCAPE
	Pixel size = 30.0; // 30.0
	Pixel delta = 3.0; // 3.0

	Array< Pixel, 3 > C( 2*size+1, 2*size+1, 2*size+1 );
	C = 0;

	// get rid of translation
	Pixel translation[3];
	bestTransform->GetPosition(translation);
	translation[0] *= -1;
	translation[1] *= -1;
	translation[2] *= -1;
	bestTransform->PostMultiply();
	bestTransform->Translate( translation );

	cccMetric.setSeedTransform( bestTransform );
	vnl_vector< double > seed(3);
	seed[0] = seed[1] = seed[2] = 0;

	Pixel minimizer = cccMetric.f(seed);
	int minX = 0;
	int minY = 0;
	int minZ = 0;

	cccMetric.setRotationSearchRestriction( 180 );

	Pixel saved[3];

	int indexI, indexJ, indexK;
	indexI = indexJ = indexK = 0;
	Pixel i, j, k;
	for ( i = - size; i <= size; i++ ){
		for ( j = - size; j <= size; j++ ){
			for ( k = - size; k <= size; k++ ){
				vnl_vector< double > current(3);
				current[0] = seed[0] + 2 * delta * i;
				current[1] = seed[1] + delta * j;
				current[2] = seed[2] + 2 * delta * k;
				C(indexI,indexJ,indexK) = cccMetric.f(current);
				if ( ( indexI == 53 ) && ( indexJ == 28 ) && ( indexK == 38 ) ){
					cout << current[0] << ", " << current[1] << ", " << current[2] << endl;
					saved[0] = current[0];
					saved[1] = current[1];
					saved[2] = current[2];
					cout << cccMetric.f(current) << endl;
				}
				if ( C(indexI,indexJ,indexK) < minimizer ){
					minimizer = C(indexI,indexJ,indexK);
					minX = current[0]; minY = current[1]; minZ = current[2];
				}
				// cout << "[" << current[0] << "," << current[1] << "," << current[2] << "] - " << C(indexI,indexJ,indexK) << endl; 
				indexK++;
			}
			indexK = 0;
			indexJ++;
		}
		indexJ = 0;
		indexI++;
	}
	cout << "New minima = [" << minX << "," << minY << "," << minZ << "] - " << minimizer << endl;
	nbfMrcWriter mw;
	mw.setFileName("ccc_3_degrees.mrc");
	vtkImageData * data1 = vtkImageData :: New();
	nbfVTKInterface :: blitzToVtk( C, data1 );
	mw.write( data1 );
	data1->Delete();
	return;
#endif

	//bestTransform->RotateX( saved[0] );
	//bestTransform->RotateY( saved[1] );
	//bestTransform->RotateZ( saved[2] );
	//cccMetric.executeFourierNewHalf( bestTransform );
	//this->transform->DeepCopy( bestTransform );

	bestTransform->Delete();

	//cout << *(this->transform->GetMatrix()) << endl;
#ifdef NBF_VERBOSE
	//nbfMrcWriter mw;
	this->fourierFilter->wedgeOff();

	this->wedgedInput1->getImage( this->input1 );
	//this->imageFilter->fileMaskOn();
	//this->imageFilter->execute( this->input1 );
	Array< double, 3 > L;
	nbfVTKInterface :: vtkToBlitzReference( this->input1, L );

	Array< float, 3 > Lf( L.shape() );
	vtkImageData * data = vtkImageData :: New();
	Lf = cast<float>(L);
	nbfVTKInterface :: blitzToVtk( Lf, data );
	nbfMrcWriter mw;
	mw.setFileName("i1.mrc");
	mw.write( data );
	
	this->imageFilter->execute( this->input1 );
	nbfVTKInterface :: vtkToBlitzReference( this->input1, L );

	this->fourierFilter->executeHalf( L);

	Lf = cast<float>(L);
	nbfVTKInterface :: blitzToVtk( Lf, data );
	mw.setFileName("i1f.mrc");
	mw.write( data );

	// save spherical wedge image + corresponding wedge image
	this->wedgedInput2->getImage( this->input2, this->transform );
	nbfVTKInterface :: vtkToBlitzReference( this->input2, L );

	Lf = cast<float>(L);
	nbfVTKInterface :: blitzToVtk( Lf, data );

	mw.setFileName("i2.mrc");
	mw.write( data );
	//this->imageFilter->execute( this->input2 );
	//nbfVTKInterface :: vtkToBlitzReference( this->input2, L );
	//this->fourierFilter->executeHalf( L);
	//w.setFileName("i2f.matlab");
	//w.write(L);

	this->imageFilter->execute( this->input2 );
	nbfVTKInterface :: vtkToBlitzReference( this->input2, L );
	this->fourierFilter->executeHalf( L);
	Lf = cast<float>(L);
	nbfVTKInterface :: blitzToVtk( Lf, data );
	mw.setFileName("i2f.mrc");
	mw.write( data );

	this->cartesianToSpherical( this->input2, Pp2 );
	w.setFileName("c.matlab");
	w.write(Pp2);

	this->wedgedInput2->getSphericalWedgeImage( Pp2, this->transform, this );
	w.setFileName("wc.matlab");
	w.write(Pp2);
	data->Delete();


#endif
}

template< class Pixel >
void nbfProjectionRotationMetric3D< Pixel > :: initializeHarmonics( int bwIn, int bwOut)
{
	if ( this->harmonicsInitialized == false ){
		this->bwIn = bwIn;
		this->bwOut = bwOut;
		this->degLim = this->bwOut - 1;
		int n = 2 * this->bwIn;

		//sigR = (REAL *) calloc( n * n, sizeof(REAL) );
		this->sigI = (REAL *) calloc( n * n, sizeof(REAL) );

		// IMPROVE
		//Array< Pixel, 2 > Ci( this->sigI, A.shape(), neverDeleteData );
		//Ci = 0;

		this->so3Sig = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * (8*this->bwOut*this->bwOut*this->bwOut) );
		this->workspace1 = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * (8*this->bwOut*this->bwOut*this->bwOut) );
		this->workspace2 = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * ((14*this->bwIn*this->bwIn) + (48 * this->bwIn)));
		this->workspace3 = (REAL *) malloc( sizeof(REAL) * (12*n + n*this->bwIn));
		this->sigCoefR = (REAL *) malloc( sizeof(REAL) * this->bwIn * this->bwIn ) ;
		this->sigCoefI = (REAL *) malloc( sizeof(REAL) * this->bwIn * this->bwIn ) ;
		this->patCoefR = (REAL *) malloc( sizeof(REAL) * this->bwIn * this->bwIn ) ;
		this->patCoefI = (REAL *) malloc( sizeof(REAL) * this->bwIn * this->bwIn ) ;
		this->so3Coef = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * ((4*this->bwOut*this->bwOut*this->bwOut-this->bwOut)/3) ) ;


		this->seminaive_naive_tablespace =
			(double *) malloc(sizeof(double) *
			(Reduced_Naive_TableSize(this->bwIn,this->bwIn) +
			Reduced_SpharmonicTableSize(this->bwIn,this->bwIn)));

		n = 2 * bwIn;

		/****
		At this point, check to see if all the memory has been
		allocated. If it has not, there's no point in going further.
		****/

		if ( (this->seminaive_naive_tablespace == NULL) ||
			(this->sigI == NULL) ||
			(this->so3Coef == NULL) ||
			(this->workspace1 == NULL) || (this->workspace2 == NULL) ||
			(this->workspace3 == NULL) ||
			(this->sigCoefR == NULL) || (this->sigCoefI == NULL) ||
			(this->patCoefR == NULL) || (this->patCoefI == NULL) ||
			(this->so3Sig == NULL) )
		{
			perror("Error in allocating memory");
			exit( 1 ) ;
		}

		/* create plan for inverse SO(3) transform */
		n = 2 * this->bwOut ;
		this->howmany = n*n ;
		this->idist = n ;
		this->odist = n ;
		this->rank = 2 ;
		this->inembed[0] = n ;
		this->inembed[1] = n*n ;
		this->onembed[0] = n ;
		this->onembed[1] = n*n ;
		this->istride = 1 ;
		this->ostride = 1 ;
		this->na[0] = 1 ;
		this->na[1] = n ;

		this->p1 = fftw_plan_many_dft( this->rank, this->na, this->howmany,
			this->workspace1, this->inembed,
			this->istride, this->idist,
			this->so3Sig, this->onembed,
			this->ostride, this->odist,
			FFTW_FORWARD, FFTW_ESTIMATE );

		// fprintf(stdout,"Generating seminaive_naive tables...\n");
		this->seminaive_naive_table = SemiNaive_Naive_Pml_Table(this->bwIn, this->bwIn,
			this->seminaive_naive_tablespace,
			(double *) this->workspace2);
		this->harmonicsInitialized = true;
	}
}

template< class Pixel >
void nbfProjectionRotationMetric3D< Pixel > :: soft( Array< REAL, 2 > & A, Array< REAL, 2 > & B, Array< REAL, 3 > & C, bool inverse )
{
	int bandwidth = A.rows() / 2;

	if ( bandwidth == 256 )
		this->initializeHarmonics( bandwidth, bandwidth / 2 );
	else
		this->initializeHarmonics( bandwidth, bandwidth );

	//int maxloc, ii, jj, kk ;
	
  //bwIn = atoi( argv[3] );
  //bwOut = atoi( argv[4] );
  //degLim = atoi( argv[5] );

  //bwIn = A.rows() / 2;
  //bwOut = bwIn;
  //degLim = bwIn - 1;

  //bwIn = 8;
  //bwOut = bwIn;
  //degLim = bwIn - 1;

  //printf("Reading in signal file\n");
  /* read in SIGNAL samples */
  /* first the signal */
  //FILE * fp = fopen("C:\\home\\code\\soft10\\randomS2sigA_bw8.dat","r");
  //for ( i = 0 ; i < n * n ; i ++ )
  //  {
  //    fscanf(fp,"%lf", sigR + i);
  //    fscanf(fp,"%lf", sigI + i);
  //  }
  //fclose( fp );

  this->sigR = A.dataZero();
  Array< REAL, 2 > Ci( this->sigI, A.shape(), neverDeleteData );
  Ci = 0;

  //Array< Pixel, 2 > C1( sigR, A.shape(), neverDeleteData );

  //nbfMatlabWriter w;
  //w.setFileName("p.matlab");
  //w.write(C1);

  int n = 2 * bwIn;

  //printf("now taking spherical transform of signal\n");
  FST_semi_memo( this->sigR, this->sigI,
		 this->sigCoefR, this->sigCoefI,
		 n, this->seminaive_naive_table,
		 (double *) this->workspace2, 1, this->bwIn ) ;

  //printf("Reading in pattern file\n");
  /* read in SIGNAL samples */
  /* first the signal */
  //fp = fopen("C:\\home\\code\\soft10\\randomS2sigA_bw8.dat","r");
  ////fp = fopen("C:\\home\\code\\soft10\\joder.dat","r");
  //for ( i = 0 ; i < n * n ; i ++ )
  //  {
  //    fscanf(fp,"%lf", sigR + i);
  //    fscanf(fp,"%lf", sigI + i);
  //  }
  //fclose( fp );

  sigR = B.dataZero();
  //Ci = 0;
  //Array< Pixel, 2 > C2( sigR, B.shape(), neverDeleteData );

  //w.write(C2);

  //printf("now taking spherical transform of pattern\n");
  FST_semi_memo( this->sigR, this->sigI,
		 this->patCoefR, this->patCoefI,
		 n, this->seminaive_naive_table,
		 (double *) this->workspace2, 1, this->bwIn ) ;

  //printf("freeing seminaive_naive_table and seminaive_naive_tablespace\n");
  //printf("about to combine coefficients\n");

  /* combine coefficients */
  so3CombineCoef_fftw( this->bwIn, this->bwOut, this->degLim,
		       this->sigCoefR, this->sigCoefI,
		       this->patCoefR, this->patCoefI,
		       this->so3Coef ) ;
 
  //Array< complex< Pixel >, 1 > dummy( reinterpret_cast<complex<Pixel>*>(so3Coef), (4*bwOut*bwOut*bwOut-bwOut)/3 );

  //if ( inverse == true ){
	 // so3Coef = reinterpret_cast< fftw_complex * >( C.dataZero() );
  //}
  //else{
		//C.resize( dummy.shape() );
	 //   C = dummy;
		//// return;
  //}

  //printf("about to inverse so(3) transform\n");

  /* now inverse so(3) */
  Inverse_SO3_Naive_fftw( this->bwOut,
			  this->so3Coef,
			  this->so3Sig,
			  this->workspace1,
			  this->workspace2,
			  this->workspace3,
			  &(this->p1),
			  1 ) ;
  //printf("finished inverse so(3) transform\n");

  n = 2 * this->bwOut;
  //double * distance;
  //distance = new double[n*n*n];
  //for ( int i = 0 ; i < 8*this->bwOut*this->bwOut*this->bwOut ; i ++ )
  //  {
		//distance[i] = so3Sig[i][0];
  //  }

  GeneralArrayStorage<3> storage;
  storage.ordering() = thirdDim, firstDim, secondDim;
  //Array< Pixel, 3 > D( distance, shape(n,n,n), neverDeleteData, storage );

  //C.resize( D.shape() );
  //C = D;

  Array< complex< double >, 3 > H( reinterpret_cast< complex<double> * >(so3Sig), shape(n,n,n), neverDeleteData, storage );

  C.resize( H.shape() );
  //C = cast< Pixel >( real(H) );
  C = real(H);

  //delete [] distance;
  //w.write(C);

  //ii = floor( maxloc / (4.*bwOut*bwOut) );
  //tmp = maxloc - (ii*4.*bwOut*bwOut);
  //jj = floor( tmp / (2.*bwOut) );
  //tmp = maxloc - (ii *4*bwOut*bwOut) - jj*(2*bwOut);
  //kk = tmp ;

  //TinyVector< int, 3 > con = maxIndex(D);
  //cout << con << ", " << D(con) << endl;

  //printf("ii = %d\tjj = %d\tkk = %d\n", ii, jj, kk);

  //printf("alpha = %f\nbeta = %f\ngamma = %f\n",
	 //M_PI*jj/((REAL) bwOut),
	 //M_PI*(2*ii+1)/(4.*bwOut),
	 //M_PI*kk/((REAL) bwOut) );

  //printf("alpha = %f\nbeta = %f\ngamma = %f\n",
	 // M_PI*jj/((REAL) bwOut) * vtkMath::RadiansToDegrees(),
	 //M_PI*(2*ii+1)/(4.*bwOut) * vtkMath::RadiansToDegrees(),
	 //M_PI*kk/((REAL) bwOut) * vtkMath::RadiansToDegrees() );

  //printf("alpha = %f\nbeta = %f\ngamma = %f\n",
	 // M_PI*con[0]/((REAL) bwOut) * vtkMath::RadiansToDegrees(),
	 //M_PI*(2*con[1]+1)/(4.*bwOut) * vtkMath::RadiansToDegrees(),
	 //M_PI*con[2]/((REAL) bwOut) * vtkMath::RadiansToDegrees() );
}

template< class Pixel >
void nbfProjectionRotationMetric3D< Pixel > :: softScaled( Array< REAL, 2 > & F1, Array< REAL, 2 > & F2, Array< REAL, 3 > & C, bool inverse )
{
	this->soft( F1, F2, C, inverse );

	// function to compute integrals on the unit sphere
	Array< REAL, 2 > theta( F1.shape() );
	firstIndex i;
	theta = sin( vtkMath::Pi() * ( 2 * i + 1 ) / 4 / this->B );

	Array< REAL, 2 > Rmin( F1.shape() );
	Array< REAL, 2 > Rmax( F1.shape() );

	REAL minValue, maxValue;

	// adjust volume range
	this->rotatedSphericalImageMinMax(C,F2,Rmin,Rmax,minValue,maxValue);
	
	REAL wMin = sum( F1 * Rmin * theta );
	REAL wMax = sum( F1 * Rmax * theta );

	C = ( C - minValue ) / ( maxValue - minValue ) * ( wMax - wMin ) + wMin;
}

template< class Pixel >
void nbfProjectionRotationMetric3D< Pixel > :: conditionSphericalWindow( Array< Pixel, 2 > & I )
{
	// do dilation to fill-in the small wedge gap at low frequencies
	//Array< bool, 2 > Ainput( I.shape() ), Aoutput( I.shape() );
	//BordStrategyMirrorSimple< bool, 2 > bsForAinput( Ainput, 1 );

	//Ainput = I > 0;
	//bsForAinput.refresh();
	//Aoutput = dilate2D(Ainput);
	//I = where( Aoutput == true, 1, 0 );

	vtkImageData * window = vtkImageData::New();
	vtkImageGaussianSmooth * filter = vtkImageGaussianSmooth::New();

	// filter down the window functions (to improve spherical harmonics representation)
	nbfVTKInterface::blitzToVtk( I, window );
	filter->SetInput( window );
	filter->SetRadiusFactors(1,1,1);
	filter->Update();
	nbfVTKInterface::vtkToBlitz( filter->GetOutput(), I );

	filter->Delete();
	window->Delete();
}

template< class Pixel >
void nbfProjectionRotationMetric3D< Pixel > :: harmonics( Array< Pixel, 2 > & A, vtkTransform * t, Array< Pixel, 2 > & R )
{
	if ( t == NULL ){
		if ( sum( abs( A.shape() - R.shape() ) ) > 0 ){
			cerr << "ERROR - Wedges geometries do not match. In " << __FILE__ "," << __LINE__ << endl;
		} else {
			R = A;
		}
		return;
	}

	vtkTransform * rotation = vtkTransform :: New();
	rotation->SetInput( t );
	double position[3];
	t->GetPosition(position);
	position[0] = - position[0];
	position[1] = - position[1];
	position[2] = - position[2];
	rotation->PostMultiply();
	rotation->Translate(position);

	Pixel phi, theta, psi;

	// see if denominator is zero
	if ( abs( rotation->GetMatrix()->GetElement(2,0) ) < numeric_limits< Pixel > :: min() ){
		// see if numerator is also zero
		if ( abs( rotation->GetMatrix()->GetElement(2,1) ) < numeric_limits< Pixel > :: min() ){
			psi = 0;
		// else, sign of numerator determines angle
		} else if ( rotation->GetMatrix()->GetElement(2,1) < - numeric_limits< Pixel > :: min() ){
			psi = vtkMath :: Pi() / 2.0;
		} else {
			psi = - vtkMath :: Pi() / 2.0;
		}
	} else {
		psi = atan( - rotation->GetMatrix()->GetElement(2,1) / rotation->GetMatrix()->GetElement(2,0) );
	}

	// see if denominator is zero
	if ( abs( rotation->GetMatrix()->GetElement(0,2) ) < numeric_limits< Pixel > :: min() ){
		// see if numerator is also zero
		if ( abs( rotation->GetMatrix()->GetElement(1,2) ) < numeric_limits< Pixel > :: min() ){
			phi = 0;
		// else, sign of denominator determines angle
		} else if ( rotation->GetMatrix()->GetElement(1,2) > numeric_limits< Pixel > :: min() ){
			phi = vtkMath :: Pi() / 2.0;
		} else {
			phi = - vtkMath :: Pi() / 2.0;
		}
	} else {
		phi = atan( rotation->GetMatrix()->GetElement(1,2) / rotation->GetMatrix()->GetElement(0,2) );
	}

	theta = acos( rotation->GetMatrix()->GetElement(2,2) );

	Pixel psi_1, psi_2;
	psi_1 = psi;
	if ( psi > 0 ){
		psi_2 = psi - vtkMath::Pi();
	} else {
		psi_2 = psi + vtkMath::Pi();
	}

	Pixel phi_1, phi_2;
	phi_1 = phi;
	if ( phi > 0 ){
		phi_2 = phi - vtkMath::Pi();
	} else {
		phi_2 = phi + vtkMath::Pi();
	}

	// theta > 0
	if ( rotation->GetMatrix()->GetElement(2,1) * psi_1 > 0 ){
		psi = psi_1;
	} else {
		psi = psi_2;
	}

	if ( rotation->GetMatrix()->GetElement(1,2) * phi_1 > 0 ){
		phi = phi_1;
	} else {
		phi = phi_2;
	}

	if ( cos( psi ) * rotation->GetMatrix()->GetElement(2,0) > 0 ){
		// theta < 0
		theta = - theta;
		if ( rotation->GetMatrix()->GetElement(2,1) * psi_1 < 0 ){
			psi = psi_1;
		} else {
			psi = psi_2;
		}

		if ( rotation->GetMatrix()->GetElement(1,2) * phi_1 < 0 ){
			phi = phi_1;
		} else {
			phi = phi_2;
		}
	}

	Array< REAL, 2 > Areal( A.shape() );
	Array< REAL, 2 > Rreal( A.shape() );
	Areal = cast<REAL>(A);

	this->harmonics( Areal, - psi, - theta, - phi, Rreal );
	R = cast<Pixel>(Rreal);
}

template< class Pixel >
void nbfProjectionRotationMetric3D< Pixel > :: harmonics( Array< REAL, 2 > & A, REAL alpha, REAL beta, REAL gamma, Array< REAL, 2 > & R )
{
	int bwIn, bwOut, degOut ;

	bwIn = this->B;
	bwOut = this->B;
	degOut = this->B-1;

	// initialize SOFT structures if not done already
	if ( this->rotation_spherical_harmonics_initialized == false ){
		this->rot_sigInR = (REAL *) malloc(sizeof(REAL)*(4*bwIn*bwIn));
		this->rot_sigInI = (REAL *) malloc(sizeof(REAL)*(4*bwIn*bwIn));
		this->rot_sigOutR = (REAL *) malloc(sizeof(REAL)*(4*bwOut*bwOut));
		this->rot_sigOutI = (REAL *) malloc(sizeof(REAL)*(4*bwOut*bwOut));

		if ( bwOut > bwIn )
			this->rot_scratch = (REAL *) malloc(sizeof(REAL)*((16*bwOut*bwOut) + (48 * bwOut)));
		else
			this->rot_scratch = (REAL *) malloc(sizeof(REAL)*((16*bwIn*bwIn) + (48 * bwIn)));

		this->rot_seminaive_naive_tablespace =
			(double *) malloc(sizeof(double) *
			(Reduced_Naive_TableSize(bwIn,bwIn) +
			Reduced_SpharmonicTableSize(bwIn,bwIn)));

		this->rot_trans_seminaive_naive_tablespace2 =
			(double *) malloc(sizeof(double) *
			(Reduced_Naive_TableSize(bwOut,bwOut) +
			Reduced_SpharmonicTableSize(bwOut,bwOut)));

		this->rot_seminaive_naive_tablespace2 =
			(double *) malloc(sizeof(double) *
			(Reduced_Naive_TableSize(bwOut,bwOut) +
			Reduced_SpharmonicTableSize(bwOut,bwOut)));


		/****
		At this point, check to see if all the memory has been
		allocated. If it has not, there's no point in going further.
		****/

		if ( (this->rot_scratch == NULL) || 
			(this->rot_sigInI == NULL ) ||
			(this->rot_sigOutR == NULL ) || (this->rot_sigOutI == NULL ) ||
			(this->rot_seminaive_naive_tablespace == NULL) ||
			(this->rot_trans_seminaive_naive_tablespace2 == NULL) )
		{
			perror("Error in allocating memory");
			exit( 1 ) ;
		}


		//fprintf(stdout,"Generating seminaive_naive tables...\n");
		this->rot_seminaive_naive_table = SemiNaive_Naive_Pml_Table(bwIn, bwIn,
			this->rot_seminaive_naive_tablespace,
			this->rot_scratch);


		//fprintf(stdout,"Generating seminaive_naive tables...\n");
		this->rot_seminaive_naive_table2 = SemiNaive_Naive_Pml_Table(bwOut, bwOut,
			this->rot_seminaive_naive_tablespace2,
			this->rot_scratch);


		//fprintf(stdout,"Generating trans_seminaive_naive tables...\n");
		this->rot_trans_seminaive_naive_table2 =
			Transpose_SemiNaive_Naive_Pml_Table(this->rot_seminaive_naive_table2,
			bwOut, bwOut,
			this->rot_trans_seminaive_naive_tablespace2,
			this->rot_scratch);

		this->rotation_spherical_harmonics_initialized = true;
	}

	// convert input image into SOFT structure
	for ( int i = 0; i < A.size(); i++ ){
		*(this->rot_sigInR+i) = *( A.dataZero() + i );
		*(this->rot_sigInI+i) = 0;
	}

	//fprintf(stdout,"about to rotate ...\n");

	// rotate image in the unit sphere
	rotateFct( bwIn, bwOut, degOut,
		this->rot_sigInR, this->rot_sigInI,
		this->rot_sigOutR, this->rot_sigOutI,
		alpha, beta, gamma,
		this->rot_scratch,
		this->rot_seminaive_naive_table,
		this->rot_trans_seminaive_naive_table2 ) ;

	//fprintf(stdout,"finished rotating ...\n");

	// store result
	for ( int i = 0; i < R.size(); i++ ){
		*( R.dataZero() + i ) = this->rot_sigOutR[i];
	}
}


template< class Pixel >
void nbfProjectionRotationMetric3D< Pixel > :: rotatedSphericalImageMinMax( Array< REAL, 3 > & C, Array< REAL, 2 > & F2, Array< REAL, 2 > & Rmin, Array< REAL, 2 > & Rmax, REAL & minValue, REAL & maxValue )
{
	// look for maxima and minima of volume
	maxValue = - numeric_limits< REAL > :: max();
	minValue = numeric_limits< REAL > :: max();
	TinyVector< int, 3 > minPosition, maxPosition;

	Array< REAL, 3 > :: iterator iter = C.begin();
	while ( iter != C.end() ){
		if ( *iter > maxValue ){
			maxValue = *iter;
			maxPosition = iter.position();
		}
		if ( *iter < minValue ){
			minValue = *iter;
			minPosition = iter.position();
		}
		iter++;
	}

	// rotation of minima
	Pixel alpha = minPosition[0]/((REAL) this->bwOut) * vtkMath::Pi();
	Pixel beta = (2*minPosition[1]+1.0)/(4.*this->bwOut) * vtkMath::Pi();
	Pixel gamma = minPosition[2]/((REAL) this->bwOut) * vtkMath::Pi();
	this->harmonics( F2, alpha, beta, gamma, Rmin );

	// rotation of maxima
	alpha = maxPosition[0]/((REAL) this->bwOut) * vtkMath::Pi();
	beta = (2*maxPosition[1]+1.0)/(4.*this->bwOut) * vtkMath::Pi();
	gamma = maxPosition[2]/((REAL) this->bwOut) * vtkMath::Pi();
	this->harmonics( F2, alpha, beta, gamma, Rmax );
}


template< class Pixel >
void nbfProjectionRotationMetric3D< Pixel > :: getTrueLocalMaximaPositions( Array< REAL, 3 > & C, 
																			PriorityQueue< TinyVector< Pixel, 4 >, greaterPD< Pixel, 3 > , vector< TinyVector< Pixel, 4 > > > & pqNB,
																			Pixel tolerance )
{
	pqNB.clear();

	// Store all candidate rotations in vector
	vector< TinyVector< Pixel, 4 > > positions;
	
	// traverse SO3 correlation volumes
	Array< REAL, 3 > :: iterator iterC = C.begin();
	while ( iterC != C.end() ){
		// this signals location of local maxima in correlation volume
		if ( abs(*iterC) > 0 ){

			// compute rotation angles from SO3 coordinates
			Pixel alpha = iterC.position()[0] / ((REAL) this->bwOut ) * 180.0;
			Pixel beta = ( 2 * iterC.position()[1] + 1.0 ) / ( 4. * this->bwOut ) * 180.0;

			// WARNING - If sampling if too coarse (<=32) the smallest angle that can be represented is 1.4 degrees. And since this is the smallest angle possible, if the restriction is zero degrees all transformations will be invalid.
			if ( this->B <= 32 ){
				beta = ( 2 * iterC.position()[1] + 0.0 ) / ( 4. * this->bwOut ) * 180.0;
			}
			Pixel gamma = iterC.position()[2] / ((REAL) this->bwOut ) * 180.0;

			// build transformation to check compliance with constraints
			vtkTransform * t = vtkTransform::New();
			t->RotateZ( - gamma );
			t->RotateY( - beta );
			t->RotateZ( - alpha );
			
			// store only if inside admisible range (change fitness sign so we now look for minimizers)
			if ( this->isTransformValid( t ) == true ){
				positions.push_back( TinyVector< Pixel, 4 >( alpha, beta, gamma, -(*iterC) ) );
			}
			t->Delete();
		}
		++iterC;
	}

	this->eliminateEquivalentAngles( positions, pqNB, tolerance );
}

template< class Pixel >
void nbfProjectionRotationMetric3D< Pixel > :: eliminateEquivalentAngles( vector< TinyVector< Pixel, 4 > > & positions,
																		  PriorityQueue< TinyVector< Pixel, 4 >, greaterPD< Pixel, 3 > , vector< TinyVector< Pixel, 4 > > > & pqNB,
																		  Pixel tolerance )
{
	// Trim local minimas that are too close together in SO3 space.
	// Distance between two rotations is computed as norm between transformation matrices.

	for ( int i = 0; i < positions.size(); i++ ){

		// build trial transformation
		vtkTransform * t1 = vtkTransform::New();
		t1->RotateZ( - positions[i][2] );
		t1->RotateY( - positions[i][1] );
		t1->RotateZ( - positions[i][0] );

		// extract rotation matrix
		double matrix1[16];
		vtkMatrix4x4::DeepCopy( matrix1, t1->GetMatrix() );
		t1->Delete();

		for ( int j = 0; j < positions.size(); j++ ){
	
			// skip if this point was already discarded
			if ( positions[i][3] == - numeric_limits< Pixel > :: max() ){
				break;
			}

			// skip if point was already discarded
			if ( ( i != j ) && ( positions[j][3] > - numeric_limits< Pixel > :: max() ) ){
		
				// build trial transformation
				vtkTransform * t2 = vtkTransform::New();
				t2->RotateZ( - positions[j][2] );
				t2->RotateY( - positions[j][1] );
				t2->RotateZ( - positions[j][0] );

				// extract rotation matrix
				double matrix2[16];
				vtkMatrix4x4::DeepCopy( matrix2, t2->GetMatrix() );
				t2->Delete();

				// compute similarity between rotations
				double error = 0;
				for ( int m = 0; m < 16; m++ ){
					error += fabs( matrix1[m] - matrix2[m] );
				}

				// if error is below tolerance, we assume these are the same rotations 
				// we only keep the one with smaller distance
				if ( error < tolerance ){
					// only keep smaller of the two
					if ( positions[i][3] > positions[j][3] ){
						positions[i][3] = - numeric_limits< Pixel > :: max();
					} else {
						positions[j][3] = - numeric_limits< Pixel > :: max();
					}
				}
			}
		}
	}
	for ( int i = 0; i < positions.size(); i++ ){
		if ( positions[i][3] > - numeric_limits< Pixel > :: max() ){
			pqNB.push( positions[i] );
		}
	}
}


//template< class Pixel >
//Pixel nbfProjectionRotationMetric3D< Pixel > :: getWedgeOverlap()
//{
//		if ( ( this->wedgedInput1 != NULL ) && ( this->wedgedInput2 != NULL ) ){
//		
//			Array< Pixel, 2 > W1( this->Pp1saved.shape() ), W2( this->Pp1saved.shape() );
//			
//			this->wedgedInput1->getSphericalWedgeImage( W1 );
//			this->wedgedInput2->getSphericalWedgeImage( W2 , this->transform );
//			cout << sum(W1) / ( 1.0 * W1.size() ) << endl;
//			return sum( W1 * W2 ) / ( 1.0 * W1.size() );
//		}
//		else {
//			return 1;
//		}
//}
