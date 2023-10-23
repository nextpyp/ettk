#pragma once

#include <vtkMath.h>
#include <vtkImageData.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfImageWriter.h>
#include <em/nbfImageMetric.h>
#include <em/nbfFourierFilter.h>
#include <vtkImageConstantPad.h>
#include <nbfTimer.h>

template< class Pixel, int const Dim > class nbfFourierImageMetric;
template< class Pixel, int const Dim > class nbfFourierImageMetricCore;

#include <nbfLinearInterpolator.h>

#include <fftw3.h>

using namespace blitz;
 	
#include <random/uniform.h>

/** Compute normalized cross-correlation between two images in Fourier space.
	Input images may be in real or already in reciprocal representation.

	Correlation peak is located with sub-pixel accuracy, using parabolic fit (Frank p. 82) in each axis.

	@todo Improve peak search computation, i.e. compute two most significant peaks and 
	use difference for validation purposes. 
*/
template< class Pixel, int const Dim >
class nbfCorrelationImageMetric : public nbfImageMetric< Pixel, Dim >
{
public:

	/// Constructor with wedge object argument (neccesary to carry out FFT computations).
	nbfCorrelationImageMetric( nbfImageFilter< Pixel, Dim > * = NULL, nbfFourierFilter< Pixel, Dim > * = NULL );

	virtual ~nbfCorrelationImageMetric();

	virtual int getId(){ return NBF_IMAGE_METRIC_CORR; }

	/// Redefine
	virtual void execute();

	/// Optimal transformation.
	vtkTransform * getTransform( int index = 0 ){
		if ( this->candidateTransforms.size() > 0 ){
			double matrix[16];
			for ( int i = 0; i < 16; i++ ){
				matrix[i] = this->candidateTransforms[index](i);
			}
			vtkTransform * t = vtkTransform::New();
			t->Concatenate( matrix );
			this->transform->DeepCopy( t );
			t->Delete();
		}
		return this->transform;
	}

	/// Corrected best match.
	vtkImageData * getAligned(){ return this->aligned; }

	/// Corrected best match (fourier).
	vtkImageData * getAlignedFourier(){ return this->alignedFourier; }

	//void getCCCImage( Array< Pixel, Dim > & A ){ A.reference( this->cccImage ); }

	static TinyVector< Pixel, Dim > locateMaximaSubPixelParabolic( Array< Pixel, Dim > &, TinyVector< int, Dim > & );
	static TinyVector< Pixel, 2 > locateMaximaSubPixelCentroid( Array< Pixel, 2 > &, TinyVector< int, 2 > & );
	static TinyVector< Pixel, 3 > locateMaximaSubPixelCentroid( Array< Pixel, 3 > &, TinyVector< int, 3 > & );

	/// Redefine from nbfImageMetric
	Pixel getWedgeOverlap( vtkTransform * );

	// MPI - result Array stores all information of the alignment: distance, ccc, wedge overlap and matrix[16]
	// compute distace matrices
	void getDistances( vector< nbfWedgedSubImage3D< Pixel > > &, Array< Pixel, 4 > &, int = 0 );
	void getDistances( vector< nbfWedgedAverageImage3D< Pixel > > &, Array< Pixel, 4 > & , int = 0);
	// compute inter-vector distances
	void getDistances( vector< nbfWedgedAverageImage3D< Pixel > > &, vector< nbfWedgedSubImage3D< Pixel > > &, Array< Pixel, 4 > &, int = 0 );
	void getDistances( vector< nbfWedgedImage3D< Pixel > * > &, vector< nbfWedgedImage3D< Pixel > * > &, vector< TinyVector< int, 2 > > &, Array< Pixel, 3 > &, int = 0 );
	// void getDistances( vector< nbfWedgedImage3D< Pixel > * > &, vector< nbfWedgedImage3D< Pixel > * > &, Array< Pixel, 3 > & );
	void getDistances( vector< nbfWedgedAverageImage3D< Pixel > > &, vector< nbfWedgedAverageImage3D< Pixel > > &, Array< Pixel, 4 > &, int = 0 );
	void getDistances( vector< nbfWedgedSubImage3D< Pixel > > &, vector< nbfWedgedSubImage3D< Pixel > > &, Array< Pixel, 4 > &, int = 0 );

	void getImage( nbfWedgedAverageImage3D< Pixel > & );
	void getImageVariance( nbfWedgedAverageImage3D< Pixel > & );

	void updateAccumulatedWedgeImage( nbfWedgedAverageImage3D< Pixel > & );
	void updateSphericalWedgeImage( nbfWedgedAverageImage3D< Pixel > &, TinyVector< int, 2 > & );

	// given list of volumes compute a representation for each volume: real + imag + wedge
	void getRepresentations( stringstream &, Array< Pixel, 3 > &, int, Pixel = 1, int = 0 );
	void getDistances( stringstream &, Array< Pixel, 2 > & );

	// For CTF correction
	void get2DRepresentations( stringstream &, Array< Pixel, 3 > &, int, Pixel = 1, int = 0 );
	void get2DAverages( stringstream &, Array< Pixel, 3 > &, Array< Pixel, 3 > & );

	void finalizeMPI();
	void slaveMPI();

	void setFourierFilter( nbfFourierFilter< Pixel, Dim > * filter ){ this->fourierFilter = filter; }

	// MPI
	void makeDistanceMatrix( Array< Pixel, 4 > & );

	void refinementOn(){ this->refinement = true; }
	void refinementOff(){ this->refinement = false; }

	Pixel getGivenDistance( int index ){ this->correlationPeak = this->candidateCorrelationPeaks[index]; return this->getDistance(); }
	Pixel getGivenCorrelationPeak( int index ){	return this->candidateCorrelationPeaks[index]; }
	Pixel getGivenCorrelationScale( int index ){ return this->candidateCorrelationScales[index]; }
	Pixel getGivenWedgeOverlap( int index ){ return this->candidateWedgeOverlaps[index]; }

	void setNumberOfCandidates( int i ){ this->numberOfCandidates = i; }
	int getNumberOfCandidates(){ return this->numberOfCandidates; }

	/** Normalization of input images in Fourier domain.
	*/
	Pixel normalizeFourier( Array< complex< double >, Dim > & );
	Pixel normalizeFourierHalf( Array< complex< double >, Dim > &, bool = false );

protected:

	/// Peak search (implements parabolic fit for each axis)
	TinyVector< Pixel, Dim > locateMaxima( Array< double, Dim > & );
	void locateMaxima( Array< Pixel, Dim > &, vector< TinyVector< Pixel, Dim > > &, vector< Pixel > & );

	vtkTransform * transform;
	vtkImageData * aligned;
	vtkImageData * alignedFourier;

	// store vectors of candidates
	int numberOfCandidates;
	vector< Array< double, 1 > > candidateTransforms;
	vector< Pixel > candidateCorrelationPeaks;
	vector< Pixel > candidateCorrelationScales;
	vector< Pixel > candidateWedgeOverlaps;

	//Array< Pixel, Dim > cccImage;

	Array< complex< double >, Dim > FFT1;
	//Array< complex< Pixel >, Dim > FFT2;

	Array< complex< double >, Dim > FFT1saved;

	// MPI
	void getDistances( Array< Pixel, 2 > & );
	static const int tag_done = 1;
	static const int tag_processing = 2;
	static const int tag_averaging = 3;
	static const int tag_representation = 4;
	static const int tag_representation_distances = 5;
	static const int tag_variance = 6;
	static const int tag_update_accum_wedge = 7;
	static const int tag_update_accum_spherical_wedge = 8;
	static const int tag_representation_2D = 9;
	static const int tag_average_2D = 10;
	//static const int parameterSize = 4; // jobNumber, sizeOfObject, typeOfMetric, number of consecutive distances to compute
	static const int parameterSize = 11; // jobNumber, sizeOfObject, typeOfMetric, number of consecutive distances to compute
	static const int resultSize = 1 + 3 + 16; // index, distance, ccc, overlap, matrix[16]
	bool MPIdistanceMatrix;
	int MPIbundleNumber;

	bool refinement;
};

template< class Pixel, int const Dim  >
nbfCorrelationImageMetric< Pixel, Dim > :: nbfCorrelationImageMetric( nbfImageFilter< Pixel, Dim > * i, nbfFourierFilter< Pixel, Dim  > * f )
: nbfImageMetric< Pixel, Dim >(i,f), MPIdistanceMatrix(false)
{
	this->transform = vtkTransform::New();
	this->aligned = vtkImageData::New();
	this->alignedFourier = vtkImageData::New();

	this->refinement = true;

	this->MPIbundleNumber = 1;

	this->numberOfCandidates = 1;
}


template< class Pixel, int const Dim  >
nbfCorrelationImageMetric< Pixel, Dim > :: ~nbfCorrelationImageMetric()
{
	this->transform->Delete();
	this->aligned->Delete();
	this->alignedFourier->Delete();
}

template< class Pixel, int const Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: execute()
{
	nbfMatlabWriter w;
	w.setFileName("p.matlab");

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

		this->FFT1.resize( this->fourierFilter->blitzFFT2.shape() );
		this->FFT1 = this->fourierFilter->blitzFFT2;

		// store FFT1 to avoid recomputation
		this->FFT1saved.resize( this->FFT1.shape() );
		this->FFT1saved = this->FFT1;
		this->input1Changed = false;
	} 
	else {
		this->FFT1.resize( this->FFT1saved.shape() );
		this->FFT1 = this->FFT1saved;
	}

	// normalize input1
	this->normalizeFourier( this->FFT1 );

	if ( this->input2Changed == true ){	
		// window data in real space
		this->imageFilter->execute( this->input2 );
		this->input2Changed = false;
	}

	nbfVTKInterface::vtkToBlitzReference( this->input2, A );
	//w.write(A);
	this->fourierFilter->initializeFFTW( A.shape() );
	real( this->fourierFilter->blitzFFT1 ) = A * this->fourierFilter->shift;
	imag( this->fourierFilter->blitzFFT1 ) = 0;

	fftw_execute( this->fourierFilter->fftplan );

	// normalize input2
	this->normalizeFourier( this->fourierFilter->blitzFFT2 );

	// compute: imf1 * conj(imf2) and multiply for FFT shifting
	this->fourierFilter->blitzFFT2 = this->FFT1 * conj( this->fourierFilter->blitzFFT2 ) * cast<double>( this->fourierFilter->shift );

	// MCF - mutual correlation function
	if ( this->useMutualCorrelation == true ){
		real( this->fourierFilter->blitzFFT1 ) = pow( real( this->fourierFilter->blitzFFT2 * conj( this->fourierFilter->blitzFFT2 ) ), .25 );
		this->fourierFilter->blitzFFT2 = where( real( this->fourierFilter->blitzFFT1 ) > 0, this->fourierFilter->blitzFFT2 / real( this->fourierFilter->blitzFFT1 ), this->fourierFilter->blitzFFT2 );
	}

	// real( fftshift( ifft( imf1 .* conj(imf2) ) ) )
	fftw_execute( this->fourierFilter->ifftplan );

	// store CCC
	real( this->FFT1 ) = pow(-1.0,Dim-1.0) * real( this->fourierFilter->blitzFFT1 ) * this->fourierFilter->shift;

	// locate maxima position with sub-pixel accuracy
	Array< double, Dim > T( real( this->FFT1 ) );
	//w.write(T);
	TinyVector< double, Dim > cccPeakPosition = this->locateMaxima( T );

	if ( sum( cccPeakPosition ) == 0 ){
		this->correlationPeak = numeric_limits< Pixel > :: max();
		return;
	}

	// interpolate maxima value and assign to attribute
	nbfLinearInterpolator< double, Dim > interp( T );

	// USE NEGATIVE CORRELATION
	this->correlationPeak = - interp.interpolateSingle( cccPeakPosition ) / ( 1.0 * this->FFT1.size() );

	// transform to translation coordinates (remove offset)
	//cccPeakPosition = - ( floor( this->cccImage.shape() / 2.0 ) + 1.0 - cccPeakPosition ) + 1;
	cccPeakPosition = - ( floor( this->FFT1.shape() / 2.0 ) - cccPeakPosition );

	// reset transformation to identity
	vtkTransform * t = vtkTransform::New();
	this->transform->DeepCopy( t );
	t->Delete();

	TinyVector< Pixel, 3 > save;
	save[firstDim] = cccPeakPosition[firstDim];
	save[secondDim] = cccPeakPosition[secondDim];
	if ( Dim < 3 ){
		save[thirdDim] = 0;
	} else{
		save[thirdDim] = cccPeakPosition[thirdDim];
	}
	this->transform->Translate( save(0), save(1), save(2) );
}


template< class Pixel, int const Dim  >
TinyVector< Pixel, Dim > nbfCorrelationImageMetric< Pixel, Dim > :: locateMaxima( Array< double, Dim > & A )
{
	// peak search using parabolic axis-wise fit

	TinyVector< int, Dim > cccPeakPosition = maxIndex( A );

	if ( A.isInRange( cccPeakPosition ) == false ){
	//if ( ( cccPeakPosition[0] == numeric_limits< int > :: max() ) || ( cccPeakPosition[0] == - numeric_limits< int > :: max() ) ){
		cerr << "WARNING - cannot find maxima of CCC inside image." << endl;
		TinyVector< Pixel, Dim > t; t = 0;
		return t;
		//return ( A.shape() / 2.0 );
	}

	TinyVector< Pixel, Dim > res = cccPeakPosition;
	return res;

	//return nbfCorrelationImageMetric< double, Dim > :: locateMaximaSubPixelParabolic( A, cccPeakPosition );
	//return nbfCorrelationImageMetric< double, Dim > :: locateMaximaSubPixelCentroid( A, cccPeakPosition );
}


template< class Pixel, int const Dim  >
TinyVector< Pixel, Dim > nbfCorrelationImageMetric< Pixel, Dim > :: locateMaximaSubPixelParabolic( Array< Pixel, Dim > & A, TinyVector< int, Dim > & p )
{
	TinyVector< Pixel, Dim > centroid = p;
#if 1
	Pixel fb = A(p);
	Pixel fa, fc;

	TinyVector< int, Dim > delta;
	TinyVector< int, Dim > current;

	// traverse in each dimension
	for ( int i = 0; i < Dim; i++ ){
		fa = fc = fb;
		delta = 0;
		delta[ i ] = 1;

		current = p - delta;
		bool specialA = true;
		if ( A.isInRange( current ) ){
			if ( A(current) > - numeric_limits< Pixel > :: max() ){
				fa = A( current );
				specialA = false;
			}
		}
		current = p + delta;
		bool specialC = true;
		if ( A.isInRange( current ) ){
			if ( A(current) > - numeric_limits< Pixel > :: max() ){
				fc = A( current );
				specialC = false;
			}
		}

		if ( ( specialA == false ) && ( specialC == false ) ){
			//Pixel f = min( fa, fc ) - 1;
            centroid[i] -= .5 * ( fa - fc ) / ( 2.0 * fb - fa - fc );
			//if ( f > 0 ){
			//	centroid[i] = ( fa * ( p[i] - 1 ) + fb * p[i] + fc * ( p[i] + 1 ) ) / ( fa + fb + fc );
			//}
			//else {
			//	centroid[i] = ( ( fa - f ) * ( p[i] - 1 ) + ( fb - f ) * p[i] + ( fc - f ) * ( p[i] + 1 ) ) / ( fa + fb + fc - 3 * f );
			//}
		}
	}
#endif
	return centroid;
}

template< class Pixel, int const Dim  >
TinyVector< Pixel, 2 > nbfCorrelationImageMetric< Pixel, Dim > :: locateMaximaSubPixelCentroid( Array< Pixel, 2 > & A, TinyVector< int, 2 > & p )
{
	TinyVector< Pixel, 2 > centroid = 0;
	Pixel accum = 0;

	for ( int r = -1; r <=1; r++ ){
		for ( int c = -1; c <=1; c++ ){
			TinyVector< Pixel, 2 > pos(r,c);
			pos += p;
			if ( A.isInRange(pos) ){
				Pixel f = A(pos);
				if ( f > - numeric_limits< Pixel > :: max() ){
					centroid = centroid + f * pos ;
					accum += f;
				}
			}
		}
	}
	return centroid / accum;
}

template< class Pixel, int const Dim  >
TinyVector< Pixel, 3 > nbfCorrelationImageMetric< Pixel, Dim > :: locateMaximaSubPixelCentroid( Array< Pixel, 3 > & A, TinyVector< int, 3 > & p )
{
	TinyVector< Pixel, 3 > centroid = 0.0;
	Pixel accum = 0;

	vector< TinyVector< int, 3 > > points;
	points.push_back( p );
	TinyVector< int, 3 > pos = p + TinyVector< int, 3 >( 0, 0, 1 );
	points.push_back( pos );
	pos = p + TinyVector< int, 3 >( 0, 0, -1 );
	points.push_back( pos );
	pos = p + TinyVector< int, 3 >( 1, 0, 0 );
	points.push_back( pos );
	pos = p + TinyVector< int, 3 >( -1, 0, 0 );
	points.push_back( pos );
	pos = p + TinyVector< int, 3 >( 0, 1, 0 );
	points.push_back( pos );
	pos = p + TinyVector< int, 3 >( 0, -1, 0 );
	points.push_back( pos );

	for ( int i = 0; i < points.size(); i++ ){
		if ( A.isInRange( points[i] ) ){
			Pixel f = A( points[i] );
			if ( f >= 0 ){
				centroid = centroid + f * points[i] ;
				accum += f;
			} else {
				if ( f != - numeric_limits< Pixel > :: max() ){
					cerr << "WARNING - Centroid algorithm gets negative distance when expecting positive values. Continuing. f = " << f << endl;
				}
			}
		}
	}

	//for ( int r = -1; r <=1; r++ ){
	//	for ( int c = -1; c <=1; c++ ){
	//		for ( int d = -1; d <=1; d++ ){
	//			TinyVector< Pixel, 3 > pos(r,c,d);
	//			pos += p;
	//			if ( A.isInRange(pos) ){
	//				Pixel f = A(pos);
	//				if ( f >= 0 ){
	//					centroid = centroid + f * pos ;
	//					accum += f;
	//				} else {
	//					if ( f != - numeric_limits< Pixel > :: max() ){
	//						cerr << "WARNING - Centroid algorithm gets negative distance when expecting positive values. Continuing. f = " << f << endl;
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
	return centroid / accum;
}

template< class Pixel, int const Dim  >
void nbfCorrelationImageMetric< Pixel, Dim > :: locateMaxima( Array< Pixel, Dim > & A, vector< TinyVector< Pixel, Dim > > & local, vector< Pixel > & ccc )
{
	// peak search using parabolic axis-wise fit

	typename Array< Pixel, Dim > :: iterator iter = A.begin();
	int size = 2;
	Range R(-size,size);
	while ( iter != A.end() ){
		TinyVector< int, Dim > p = iter.position();
		Range I( max( p[0] - size, 0 ), min( p[0] + size, A.ubound(0) ) );
		Range J( max( p[1] - size, 0 ), min( p[1] + size, A.ubound(1) ) );
		Range K( max( p[2] - size, 0 ), min( p[2] + size, A.ubound(2) ) );
		if ( (*iter) == max( A( I,  J,  K ) ) ){
			local.push_back( p );
			ccc.push_back( (*iter) );
		}
		++iter;
	}
}

// Assume complex input and compute normalized version for correlation computation
// This is an in-place filter, overwrites the input
template< class Pixel, int const Dim  >
Pixel nbfCorrelationImageMetric< Pixel, Dim > :: normalizeFourier( Array< complex< double >, Dim > & input )
{
	this->fourierFilter->execute( input );
	
	// kill DC component
	//( *input.begin() ) = 0;
	TinyVector< int, Dim > center = floor( input.shape() / 2 );
	//TinyVector< int, Dim > center = floor( ( input.shape() - 1 ) / 2.0 );
	input( center ) = 0;

	// calculate normalization constant
	Pixel norm = sqrtf( mean( real( input * conj(input) ) ) );
	
	input /= norm;
	return norm;
}

// Assume complex input and compute normalized version for correlation computation
// This is an in-place filter, overwrites the input
template< class Pixel, int const Dim  >
Pixel nbfCorrelationImageMetric< Pixel, Dim > :: normalizeFourierHalf( Array< complex< double >, Dim > & input, bool bypass )
{
	//this->fourierFilter->executeHalf( input );

	// kill DC component
	TinyVector< int, Dim > center;
	if ( Dim == 2 ){
		center = TinyVector< int, Dim >( input.rows() / 2, input.ubound(secondDim) );
	}
	if ( Dim == 3 ){
		center = TinyVector< int, Dim >( input.rows() / 2, input.cols() / 2, input.ubound(thirdDim) );
	}
	input( center ) = 0;

	// calculate normalization constant
	Pixel norm;
	if ( Dim == 2 ){
		Array< complex< double >, Dim > inputv( input( Range :: all(), Range(0, input.cols() - 2 ) ) );
		norm = 2 * sum( real( inputv * conj(inputv) ) );
		Array< complex< double >, 1 > inputu( input( Range :: all(), input.ubound(secondDim) ) );
		norm = sqrt( ( norm + sum( real( inputu * conj( inputu) ) ) ) / ( 1.0 * input.rows() * 2 * ( input.cols() - 1 ) ) );
	}
	if ( Dim == 3 ){
		Array< complex< double >, Dim > inputv( input( Range :: all(), Range :: all(), Range(0, input.depth() - 2 ) ) );
		norm = 2 * sum( real( inputv * conj(inputv) ) );
		Array< complex< double >, 2 > inputu( input( Range :: all(), Range :: all(), input.ubound(thirdDim) ) );
		norm = sqrt( ( norm + sum( real( inputu * conj( inputu) ) ) ) / ( 1.0 * input.rows() * input.cols() * 2 * ( input.depth() - 1 ) ) );
	}
	if ( bypass == false ){
		input /= norm;
	}
	return norm;
}

template< class Pixel, int const Dim  >
Pixel nbfCorrelationImageMetric< Pixel, Dim > :: getWedgeOverlap( vtkTransform * t )
{
	if ( ( this->wedgedInput1 != NULL ) && ( this->wedgedInput2 != NULL ) ){

		Array< Pixel, Dim > wedge1, wedge2;

		if ( this->wedgedInput1->isWedgeEffective() ){
			wedge1.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
			this->wedgedInput1->getWedgeImageHalf( wedge1 );
		}

		if ( this->wedgedInput2->isWedgeEffective() ){
			wedge2.resize( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
			this->wedgedInput2->getWedgeImageHalf( wedge2, t );
		}

		if ( this->wedgedInput1->isWedgeEffective() || this->wedgedInput2->isWedgeEffective() ){
			// update fourier filter and compute overlap only in interest region
			this->fourierFilter->initializeFFTWhalf( this->imageFilter->getPaddingFactor() * this->wedgedInput1->getDimensions() );
			this->fourierFilter->updateFilterHalf();
		} else {
			return 1;
		}

		if ( this->wedgedInput1->isWedgeEffective() && this->wedgedInput2->isWedgeEffective() ){
			return sum( wedge1 * wedge2 * this->fourierFilter->filter ) / sum( this->fourierFilter->filter );
		} else if ( this->wedgedInput1->isWedgeEffective() ){
			return sum( wedge1 * this->fourierFilter->filter ) / sum( this->fourierFilter->filter );
		} else {
			return sum( wedge2 * this->fourierFilter->filter ) / sum( this->fourierFilter->filter );
		}
	}
}

// MPI

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: getDistances( vector< nbfWedgedSubImage3D< Pixel > > & list, Array< Pixel, 4 > & D, int type )
{
	vector< nbfWedgedImage3D< Pixel > * > list1;
	vector< nbfWedgedImage3D< Pixel > * > list2;

	// build two identical lists
	for ( int i = 0; i < list.size(); i++ ){
		list1.push_back( & ( list[i] ) );
		list2.push_back( & ( list[i] ) );
	}

	vector< TinyVector< int, 2 > > indexesVector;

	// build vector of coordinates to be computed
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = i + 1; j < D.cols(); j++ ){
			indexesVector.push_back( TinyVector< int, 2 >( i, j ) );
		}
	}

	// store results in serialized vector
	Array< Pixel, 3 > serialD( indexesVector.size(), D.depth(), this->getNumberOfCandidates() );

	int count = 0;
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = i + 1; j < D.cols(); j++ ){
			serialD( count, Range::all(), Range::all() ) = D( i, j, Range::all(), Range::all() );
			count++;
		}
	}
	
	this->getDistances( list1, list2, indexesVector, serialD, type );

	count = 0;
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = i + 1; j < D.cols(); j++ ){
			D( i, j, Range::all(), Range::all() ) = serialD( count, Range::all(), Range::all() );
			count++;
		}
	}

	this->makeDistanceMatrix(D);
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: getDistances( vector< nbfWedgedAverageImage3D< Pixel > > & list, Array< Pixel, 4 > & D, int type )
{
	vector< nbfWedgedImage3D< Pixel > * > list1;
	vector< nbfWedgedImage3D< Pixel > * > list2;

	for ( int i = 0; i < list.size(); i++ ){
		list1.push_back( & ( list[i] ) );
		list2.push_back( & ( list[i] ) );
	}

	vector< TinyVector< int, 2 > > indexesVector;

	// build vector of coordinates to be computed
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = i + 1; j < D.cols(); j++ ){
			indexesVector.push_back( TinyVector< int, 2 >( i, j ) );
		}
	}

	// store results in serialized vector
	Array< Pixel, 3 > serialD( indexesVector.size(), D.depth(), this->getNumberOfCandidates() );
	
	int count = 0;
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = i + 1; j < D.cols(); j++ ){
			serialD( count, Range::all(), Range::all() ) = D( i, j, Range::all(), Range::all() );
			count++;
		}
	}

	this->getDistances( list1, list2, indexesVector, serialD, type );

	count = 0;
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = i + 1; j < D.cols(); j++ ){
			D( i, j, Range::all(), Range::all() ) = serialD( count, Range::all(), Range::all() );
			count++;
		}
	}

	this->makeDistanceMatrix(D);
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: getDistances( vector< nbfWedgedAverageImage3D< Pixel > > & ilist1, vector< nbfWedgedSubImage3D< Pixel > > & ilist2, Array< Pixel, 4 > & D, int type )
{
	vector< nbfWedgedImage3D< Pixel > * > list1;
	vector< nbfWedgedImage3D< Pixel > * > list2;

	// build first input vector
	for ( int i = 0; i < ilist1.size(); i++ ){
		list1.push_back( & ( ilist1[i] ) );
	}

	// build second input vector
	for ( int i = 0; i < ilist2.size(); i++ ){
		list2.push_back( & ( ilist2[i] ) );
	}

	vector< TinyVector< int, 2 > > indexesVector;

	// positions to be computed (all to all)
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = 0; j < D.cols(); j++ ){
			indexesVector.push_back( TinyVector< int, 2 >( i, j ) );
		}
	}

	// store results in serialized vector
	//Array< Pixel, 3 > serialD( indexesVector.size(), this->resultSize - 1, this->getNumberOfCandidates() );
	Array< Pixel, 3 > serialD( indexesVector.size(), D.depth(), this->getNumberOfCandidates() );
	
	int count = 0;
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = 0; j < D.cols(); j++ ){
			serialD( count, Range::all(), Range::all() ) = D( i, j, Range::all(), Range::all() );
			count++;
		}
	}

	this->getDistances( list1, list2, indexesVector, serialD, type );

	count = 0;
	for ( unsigned int i = 0; i < D.rows(); i++ ){
		for ( unsigned int j = 0; j < D.cols(); j++ ){
			D( i, j, Range::all(), Range::all() ) = serialD( count, Range::all(), Range::all() );
			count++;
		}
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: getDistances( vector< nbfWedgedAverageImage3D< Pixel > > & ilist1, vector< nbfWedgedAverageImage3D< Pixel > > & ilist2, Array< Pixel, 4 > & D, int type )
{
	vector< nbfWedgedImage3D< Pixel > * > list1;
	vector< nbfWedgedImage3D< Pixel > * > list2;

	// build first input vector
	for ( int i = 0; i < ilist1.size(); i++ ){
		list1.push_back( & ( ilist1[i] ) );
	}

	// build second input vector
	for ( int i = 0; i < ilist2.size(); i++ ){
		list2.push_back( & ( ilist2[i] ) );
	}

	vector< TinyVector< int, 2 > > indexesVector;

	// positions to be computed (all to all)
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = 0; j < D.cols(); j++ ){
			indexesVector.push_back( TinyVector< int, 2 >( i, j ) );
		}
	}

	// store results in serialized vector
	//Array< Pixel, 2 > serialD( indexesVector.size(), this->resultSize - 1 );
	Array< Pixel, 3 > serialD( indexesVector.size(), D.depth(), this->getNumberOfCandidates() );
	
	int count = 0;
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = 0; j < D.cols(); j++ ){
			serialD( count, Range::all(), Range::all() ) = D( i, j, Range::all(), Range::all() );
			count++;
		}
	}

	this->getDistances( list1, list2, indexesVector, serialD, type );

	count = 0;
	for ( unsigned int i = 0; i < D.rows(); i++ ){
		for ( unsigned int j = 0; j < D.cols(); j++ ){
			D( i, j, Range::all(), Range::all() ) = serialD( count, Range::all(), Range::all() );
			count++;
		}
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: getDistances( vector< nbfWedgedSubImage3D< Pixel > > & ilist1, vector< nbfWedgedSubImage3D< Pixel > > & ilist2, Array< Pixel, 4 > & D, int type )
{
	vector< nbfWedgedImage3D< Pixel > * > list1;
	vector< nbfWedgedImage3D< Pixel > * > list2;

	// build first input vector
	for ( int i = 0; i < ilist1.size(); i++ ){
		list1.push_back( & ( ilist1[i] ) );
	}

	// build second input vector
	for ( int i = 0; i < ilist2.size(); i++ ){
		list2.push_back( & ( ilist2[i] ) );
	}

	vector< TinyVector< int, 2 > > indexesVector;

	// positions to be computed (all to all)
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = 0; j < D.cols(); j++ ){
			indexesVector.push_back( TinyVector< int, 2 >( i, j ) );
		}
	}

	// store results in serialized vector
	//Array< Pixel, 2 > serialD( indexesVector.size(), this->resultSize - 1 );
	Array< Pixel, 3 > serialD( indexesVector.size(), D.depth(), this->getNumberOfCandidates() );
	
	int count = 0;
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = 0; j < D.cols(); j++ ){
			serialD( count, Range::all(), Range::all() ) = D( i, j, Range::all(), Range::all() );
			count++;
		}
	}

	this->getDistances( list1, list2, indexesVector, serialD, type );

	count = 0;
	for ( unsigned int i = 0; i < D.rows(); i++ ){
		for ( unsigned int j = 0; j < D.cols(); j++ ){
			D( i, j, Range::all(), Range::all() ) = serialD( count, Range::all(), Range::all() );
			count++;
		}
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: getImage( nbfWedgedAverageImage3D< Pixel > & input )
{
	// check if no need to update
	if ( input.isImageUpToDate() == true ){
		return;
	} else {
		input.initialize();
	}

	int num_procs;
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	// switch to uni-processor
	if ( num_procs == 1 ){
		vtkImageData * av = vtkImageData::New();
		input.getImage( av );
		av->Delete();
	} else {

		// multiprocessor
		int * params = new int[ this->parameterSize ];

		// get actual number of volumes to process (i.e. weights > 0 )
		int effectiveVolumesToAverage = sum( where( abs( input.weights( Range::all(), 0 ) ) > 0, 1, 0 ) );
		int volumesPerProcess = blitz :: extrema :: max( ceil( effectiveVolumesToAverage / ( num_procs - 1.0 ) ), 10 );

		//cout << "EffectiveVolumesToAverage = " << effectiveVolumesToAverage << endl;
		//cout << "volumesPerProcess = " << volumesPerProcess << endl;

		// start all slave processes
		int startIndex = 0;
		for ( int i = 1; i < num_procs; i++ ){

			stringstream objectData;

			// create a copy of the average object so we can manipulate its weights
			nbfWedgedAverageImage3D< Pixel > subAverage;
			subAverage = input;

			if ( startIndex > 0 ){
				subAverage.weights( Range(fromStart,startIndex-1), Range::all() ) = 0;
			}
			int j = 0;
			while ( ( j < volumesPerProcess ) && ( startIndex < subAverage.weights.rows() ) ){
				if ( abs( subAverage.weights( startIndex, 0 ) ) > 0 ){ 
					j++;
					//nonZeroElements.push_back( startIndex );
				}
				startIndex++;
			}
			if ( startIndex < subAverage.weights.rows() ){
				subAverage.weights( Range(startIndex,toEnd), Range::all() ) = 0;
			}

			//// create subAverage with only non-zero elements
			//vector< int > nonZeroElements;
			//vector< nbfWedgedSubImage3D< Pixel > > nonZeroVolumes;
			//int j = 0;
			//while ( ( j < volumesPerProcess ) && ( startIndex < input.weights.rows() ) ){
			//	if ( abs( input.weights( startIndex, 0 ) ) > 0 ){
			//		j++;
			//		nonZeroElements.push_back( startIndex );
			//		nonZeroVolumes.push_back( input.getVolumesRO()[ startIndex ] );
			//	}
			//	startIndex++;
			//}
			//Array< Pixel, 3 > nonZeroAlignments( nonZeroElements.size(), 17, input.multipleAlignments.depth() );
			//for ( int k = 0; k < nonZeroElements.size(); k++ ){
			//	nonZeroAlignments( k, 0, Range :: all() ) = input.weights( nonZeroElements[k], Range :: all() );
			//	nonZeroAlignments( k, Range(1,toEnd), Range :: all() ) = input.multipleAlignments( nonZeroElements[k], Range :: all(), Range :: all() );
			//}
			//nbfWedgedAverageImage3D< Pixel > subAverage( nonZeroVolumes );
			//subAverage.setAlignments( nonZeroAlignments );
			////////////////////// END

			subAverage.updateState();

			// store object data
			subAverage.serialize( objectData );

			// send current size of object
			int objectSize = objectData.str().size();
			params[0] = objectSize;

			//cout << "MASTER: sending nbfWedgedAverageImage3D object of size = " << objectSize << endl;
			//cout << "MASTER: sending size = " << objectSize << "\n" << objectData.str().c_str() << endl;

			// send data size and tag information
			MPI_Send ( params, this->parameterSize, MPI_INT, i, this->tag_averaging, MPI_COMM_WORLD );

			// send object data and object type information
			MPI_Send ( (char*)objectData.str().c_str(), objectSize, MPI_CHAR, i, this->tag_averaging, MPI_COMM_WORLD );

			if ( startIndex == input.weights.rows() ){
				num_procs = i + 1;
				break;
			}
		}

		//cout << "done sending " << endl;

		int jobsDone = 0;

		// build cumulative average of individual ranges
		Array< double, 3 > average( input.getDimensions() );

#ifdef NBF_AVERAGE_IN_RECIPROCAL_SPACE
		// define geometry for real FFTs
		TinyVector< int, 3 > shape( average.rows(), average.cols(), average.depth() / 2 + 1 );
		Array< double, 3 > realAverage( shape );
		realAverage = 0;
		Array< double, 3 > imagAverage( shape );
		imagAverage = 0;
		Array< double, 3 > accumulatedFourier( shape );
		accumulatedFourier = 0;
		Array< double, 4 > C;
#else
		Array< double, 3 > processResults;
		average = 0;
#endif

		//cout << "receiving loop " << endl;

		// loop while not all jobs completed
		while ( jobsDone < num_procs - 1 ){

			MPI_Status status;

			// receive notification from slaves
			MPI_Recv( params, this->parameterSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			int source = status.MPI_SOURCE;

			// we can now read the result from disk
			stringstream processFile;
			if ( source < 10 ){
				processFile << "mpi_00" << source << ".tmp";
			} else if ( source < 100 ){
				processFile << "mpi_0" << source << ".tmp";
			} else {
				processFile << "mpi_" << source << ".tmp";
			}
			nbfMatlabReader r;
			r.setFileName( processFile.str().c_str() );

#ifdef NBF_AVERAGE_IN_RECIPROCAL_SPACE
			//cout << "reading " << processFile.str().c_str() << endl;
			r.read(C);
			realAverage += C( Range::all(), Range::all(), Range::all(), 0 );
			imagAverage += C( Range::all(), Range::all(), Range::all(), 1 );
			accumulatedFourier += C( Range::all(), Range::all(), Range::all(), 2 );
			//cout << "done reading" << endl;
#else
			r.read( processResults );
			average += processResults;
#endif
			jobsDone++;
		}

		//cout << "receiving loop done." << endl;

#ifdef NBF_AVERAGE_IN_RECIPROCAL_SPACE
		// compute accumulated wedge image
		nbfWedgedAverageImage3D< Pixel > :: fourierFilter.initializeFFTWhalf( average.shape() );

		// define view of complex representation
		Array< complex< double >, 3 > FFTre( reinterpret_cast<complex<double>*>(nbfWedgedAverageImage3D< Pixel > :: fourierFilter.blitzFFT.data()), shape, neverDeleteData );
		real( FFTre ) = where( accumulatedFourier > 0, realAverage / accumulatedFourier, 0 );
		imag( FFTre ) = where( accumulatedFourier > 0, imagAverage / accumulatedFourier, 0 );
		FFTre( FFTre.rows() / 2, FFTre.cols() / 2, FFTre.depth() - 1 ) = 0;
		fftw_execute( nbfWedgedAverageImage3D< Pixel > :: fourierFilter.ifftplanreal );
		average = nbfWedgedAverageImage3D< Pixel > :: fourierFilter.blitzFFTreal / average.size() * nbfWedgedAverageImage3D< Pixel > :: fourierFilter.shift;
#endif
		// normalize to zero mean and unit variance
		Pixel varianceInside = sqrt( mean( average * average ) );
		average /= varianceInside;

		input.setAverageImage( average );

		delete [] params;
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: getImageVariance( nbfWedgedAverageImage3D< Pixel > & input )
{
	// check if no need to update
	if ( input.isImageUpToDate() == false ){
		cerr << "ERROR - Average must be computed before variance. In " << __FILE__ << "," << __LINE__ << endl;
		return;
	}

	int num_procs;
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	// switch to uni-processor
	if ( num_procs == 1 ){
		vtkImageData * av = vtkImageData::New();
		input.getImageVariance( av );
		av->Delete();
	} else {

		// multiprocessor
		int * params = new int[ this->parameterSize ];

		// get actual number of volumes to process (i.e. weights > 0 )
		int effectiveVolumesToAverage = sum( where( abs( input.weights( Range::all(), 0 ) ) > 0, 1, 0 ) );
		int volumesPerProcess = blitz :: extrema :: max( ceil( effectiveVolumesToAverage / ( num_procs - 1.0 ) ), 25 );

		// start all slave processes
		int startIndex = 0;
		for ( int i = 1; i < num_procs; i++ ){

			stringstream objectData;

			// create a copy of the average object so we can manipulate its weights
			nbfWedgedAverageImage3D< Pixel > subAverage;
			subAverage = input;

			if ( startIndex > 0 ){
				subAverage.weights( Range(fromStart,startIndex-1), Range::all() ) = 0;
			}
			int j = 0;
			while ( ( j < volumesPerProcess ) && ( startIndex < subAverage.weights.rows() ) ){
				if ( abs( subAverage.weights( startIndex, 0 ) ) > 0 ){ 
					j++;
				}
				startIndex++;
			}
			if ( startIndex < subAverage.weights.rows() ){
				subAverage.weights( Range(startIndex,toEnd), Range::all() ) = 0;
			}

			subAverage.updateState();

			// set global average image needed for the computation of variance
			subAverage.setAverageImage( input.average );

			// store object data
			subAverage.serialize( objectData );

			// send current size of object
			int objectSize = objectData.str().size();
			params[0] = objectSize;

			// send data size and tag information
			MPI_Send ( params, this->parameterSize, MPI_INT, i, this->tag_variance, MPI_COMM_WORLD );

			// send object data and object type information
			MPI_Send ( (char*)objectData.str().c_str(), objectSize, MPI_CHAR, i, this->tag_variance, MPI_COMM_WORLD );

			if ( startIndex == input.weights.rows() ){
				num_procs = i + 1;
				break;
			}
		}

		int jobsDone = 0;

		// build cumulative variance of individual ranges
		Array< double, 3 > variance( input.getDimensions() );

		Array< double, 3 > processResults;
		variance = 0;

		// loop while not all jobs completed
		while ( jobsDone < num_procs - 1 ){

			MPI_Status status;

			// receive notification from slaves
			MPI_Recv( params, this->parameterSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			int source = status.MPI_SOURCE;

			// we can now read the result from disk
			stringstream processFile;
			if ( source < 10 ){
				processFile << "mpi_00" << source << ".tmp";
			} else if ( source < 100 ){
				processFile << "mpi_0" << source << ".tmp";
			} else {
				processFile << "mpi_" << source << ".tmp";
			}
			nbfMatlabReader r;
			r.setFileName( processFile.str().c_str() );

			r.read( processResults );
			variance += processResults;

			jobsDone++;
		}

		input.setVarianceImage( variance );

		delete [] params;
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: updateAccumulatedWedgeImage( nbfWedgedAverageImage3D< Pixel > & input )
{
	// check if no need to update
	if ( input.isWedgeUpToDate() == true ){
		return;
	}

	int num_procs;
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	// switch to uni-processor
	if ( num_procs == 1 ){
		input.updateAccumulatedWedgeImage();
	} else {

		// multiprocessor
		int * params = new int[ this->parameterSize ];

		// get actual number of volumes to process (i.e. weights > 0 )
		int effectiveVolumesToProcess = sum( where( abs( input.weights( Range::all(), 0 ) ) > 0, 1, 0 ) );
		int volumesPerProcess = blitz :: extrema :: max( ceil( effectiveVolumesToProcess / ( num_procs - 1.0 ) ), 25 );

		// start all slave processes
		int startIndex = 0;
		for ( int i = 1; i < num_procs; i++ ){

			stringstream objectData;

			// create a copy of the average object so we can manipulate its weights
			nbfWedgedAverageImage3D< Pixel > subAverage;
			subAverage = input;

			if ( startIndex > 0 ){
				subAverage.weights( Range(fromStart,startIndex-1), Range::all() ) = 0;
			}
			int j = 0;
			while ( ( j < volumesPerProcess ) && ( startIndex < subAverage.weights.rows() ) ){
				if ( abs( subAverage.weights( startIndex, 0 ) ) > 0 ){ 
					j++;
				}
				startIndex++;
			}
			if ( startIndex < subAverage.weights.rows() ){
				subAverage.weights( Range(startIndex,toEnd), Range::all() ) = 0;
			}

			subAverage.updateState();

			// store object data
			subAverage.serialize( objectData );

			// send current size of object
			int objectSize = objectData.str().size();
			params[0] = objectSize;

			// send data size and tag information
			MPI_Send ( params, this->parameterSize, MPI_INT, i, this->tag_update_accum_wedge, MPI_COMM_WORLD );

			// send object data and object type information
			MPI_Send ( (char*)objectData.str().c_str(), objectSize, MPI_CHAR, i, this->tag_update_accum_wedge, MPI_COMM_WORLD );

			if ( startIndex == input.weights.rows() ){
				num_procs = i + 1;
				break;
			}
		}

		int jobsDone = 0;

		// build cumulative wedge of individual ranges
		Array< Pixel, 3 > wedge( input.getDimensions() );
		wedge = 0;

		Array< Pixel, 3 > processResults;

		// loop while not all jobs completed
		while ( jobsDone < num_procs - 1 ){

			MPI_Status status;

			// receive notification from slaves
			MPI_Recv( params, this->parameterSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			int source = status.MPI_SOURCE;

			// we can now read the result from disk
			stringstream processFile;
			if ( source < 10 ){
				processFile << "mpi_00" << source << ".tmp";
			} else if ( source < 100 ){
				processFile << "mpi_0" << source << ".tmp";
			} else {
				processFile << "mpi_" << source << ".tmp";
			}
			nbfMatlabReader r;
			r.setFileName( processFile.str().c_str() );

			r.read( processResults );
			wedge += processResults;

			jobsDone++;
		}

		input.setAccumulatedWedgeImage( wedge );

		delete [] params;
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: updateSphericalWedgeImage( nbfWedgedAverageImage3D< Pixel > & input, TinyVector< int, 2 > & size )
{
	int num_procs;
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	// switch to uni-processor
	if ( num_procs == 1 ){
		input.updateSphericalWedgeImage( size, reinterpret_cast<nbfProjectionRotationMetric3D<Pixel>*>(this) );
	} else {

		// multiprocessor
		int * params = new int[ this->parameterSize ];

		// get actual number of volumes to process (i.e. weights > 0 )
		int effectiveVolumesToProcess = sum( where( abs( input.weights( Range::all(), 0 ) ) > 0, 1, 0 ) );
		int volumesPerProcess = blitz :: extrema :: max( ceil( effectiveVolumesToProcess / ( num_procs - 1.0 ) ), 25 );

		// start all slave processes
		int startIndex = 0;
		for ( int i = 1; i < num_procs; i++ ){

			stringstream objectData;

			// create a copy of the average object so we can manipulate its weights
			nbfWedgedAverageImage3D< Pixel > subAverage;
			subAverage = input;

			if ( startIndex > 0 ){
				subAverage.weights( Range(fromStart,startIndex-1), Range::all() ) = 0;
			}
			int j = 0;
			while ( ( j < volumesPerProcess ) && ( startIndex < subAverage.weights.rows() ) ){
				if ( abs( subAverage.weights( startIndex, 0 ) ) > 0 ){ 
					j++;
				}
				startIndex++;
			}
			if ( startIndex < subAverage.weights.rows() ){
				subAverage.weights( Range(startIndex,toEnd), Range::all() ) = 0;
			}

			subAverage.updateState();

			// store object data
			subAverage.serialize( objectData );

			// send current size of object
			int objectSize = objectData.str().size();
			params[0] = objectSize;

			// send size parameters
			params[1] = size[0];
			params[2] = size[1];

			// send data size and tag information
			MPI_Send ( params, this->parameterSize, MPI_INT, i, this->tag_update_accum_spherical_wedge, MPI_COMM_WORLD );

			// send object data and object type information
			MPI_Send ( (char*)objectData.str().c_str(), objectSize, MPI_CHAR, i, this->tag_update_accum_spherical_wedge, MPI_COMM_WORLD );

			if ( startIndex == input.weights.rows() ){
				num_procs = i + 1;
				break;
			}
		}

		int jobsDone = 0;

		// build cumulative wedge of individual ranges
		Array< Pixel, 2 > wedge( size );
		wedge = 0;

		Array< Pixel, 2 > processResults;

		// loop while not all jobs completed
		while ( jobsDone < num_procs - 1 ){

			MPI_Status status;

			// receive notification from slaves
			MPI_Recv( params, this->parameterSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			int source = status.MPI_SOURCE;

			// we can now read the result from disk
			stringstream processFile;
			if ( source < 10 ){
				processFile << "mpi_00" << source << ".tmp";
			} else if ( source < 100 ){
				processFile << "mpi_0" << source << ".tmp";
			} else {
				processFile << "mpi_" << source << ".tmp";
			}
			nbfMatlabReader r;
			r.setFileName( processFile.str().c_str() );

			r.read( processResults );
			wedge += processResults;

			jobsDone++;
		}

		input.sphericalWedgeImage.resize( wedge.shape() );
		input.sphericalWedgeImage = where( wedge > max( wedge ) / 2.0, 1, 0 );
		input.setSphericalWedgeUpToDate();

		delete [] params;
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: getRepresentations( stringstream & inputFile, Array< Pixel, 3 > & R, int useRealRepresentation, Pixel bin_factor, int fold )
{
	int num_procs;
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	vector< nbfWedgedSubImage3D< Pixel > > volumeList;
	nbfWedgedSubImage3D< Pixel > :: read( inputFile.str(), volumeList );

	// reset array
	R.free();

	// switch to uni-processor
	if ( num_procs == 1 ){

		nbfFourierImageMetric< Pixel, 3 > fmetric( this->imageFilter, this->fourierFilter );
		fmetric.setToComputeOverlapNormalizedDistances( this->overlapNormalizedDistances );
		//fmetric.setToUseMutualCorrelation( this->useMutualCorrelation );
		fmetric.useMutualCorrelationOff();
		for ( int i = 0; i < volumeList.size(); i++ ){
			fmetric.setInput1( &volumeList[i] );
			if ( useRealRepresentation > 0 ){
				Array< Pixel, 1 > C;
				fmetric.getLowDimensionRepresentationReal( C, bin_factor, fold );
				if ( R.size() == 0 ){
					cerr << "Initializing representation size to " << volumeList.size() << " x " << C.rows() << " ... ";
					R.resize( volumeList.size(), C.rows(), 1 );
					cerr << "done." << endl;
				}
				R( i, Range::all(), 0 ) = C;
			} else {
				Array< complex< Pixel >, 1 > C;
				Array< Pixel, 1 > W;
				fmetric.getLowDimensionRepresentationHalf( C, W, fold );
				if ( R.size() == 0 ){
					cerr << "Initializing representation size to " << volumeList.size() << " x " << C.rows() << " x  " << 3 << " ... ";
					R.resize( volumeList.size(), C.rows(), 3 );
					cerr << "done." << endl;
				}
				R( i, Range::all(), 0 ) = real(C);
				R( i, Range::all(), 1 ) = imag(C);
				R( i, Range::all(), 2 ) = W;
			}
		}
	} else {

		// multiprocessor
		int * params = new int[ this->parameterSize ];

		int volumesPerProcess = ceil( volumeList.size() / ( num_procs - 1.0 ) );

		// start all slave processes
		for ( int i = 1; i < num_procs; i++ ){

			// compute volume range for this process
			int startIndex = ( i - 1 ) * volumesPerProcess;
			int endIndex = blitz::extrema::min( volumeList.size() - 1, i * volumesPerProcess - 1 );

			params[0] = inputFile.str().size();
			params[1] = startIndex;
			params[2] = endIndex;
			params[3] = useRealRepresentation;
			params[4] = this->overlapNormalizedDistances;
			params[5] = 0; // this->useMutualCorrelation;
			//params[5] = this->useMutualCorrelation;
			params[6] = bin_factor;
			params[7] = fold;

			// send data size and tag information
			MPI_Send ( params, this->parameterSize, MPI_INT, i, this->tag_representation, MPI_COMM_WORLD );

			// send object data and object type information
			MPI_Send ( (char*)inputFile.str().c_str(), params[0], MPI_CHAR, i, this->tag_representation, MPI_COMM_WORLD );

			if ( endIndex == volumeList.size() - 1 ){
				num_procs = i + 1;
				break;
			}
		}

		int jobsDone = 0;

		// loop while not all jobs completed
		while ( jobsDone < num_procs - 1 ){

			MPI_Status status;

			// receive notification from slaves
			MPI_Recv( params, this->parameterSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			int source = status.MPI_SOURCE;

			// we can now read the result from disk
			stringstream currentFile;
			if ( source < 10 ){
				currentFile << "mpi_00" << source << ".tmp";
			} else if ( source < 100 ){
				currentFile << "mpi_0" << source << ".tmp";
			} else {
				currentFile << "mpi_" << source << ".tmp";
			}
			nbfMatlabReader r;
			r.setFileName( currentFile.str().c_str() );
			Array< Pixel, 3 > partialR;
			r.read( partialR );

			if ( R.size() == 0 ){
				R.resize( volumeList.size(), partialR.cols(), partialR.depth() );
			}
			R( Range( params[1], params[2] ), Range::all(), Range::all() ) = partialR;

			jobsDone++;
		}
		delete [] params;
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: get2DRepresentations( stringstream & inputFile, Array< Pixel, 3 > & R, int useRealRepresentation, Pixel bin_factor, int fold )
{
	int num_procs;
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	vector< nbfWedgedSubImage3D< Pixel > > volumeList;
	nbfWedgedSubImage3D< Pixel > :: read( inputFile.str(), volumeList );

	cout << "Get total number of images ..." << endl;

	// find grand total number of volumes
	int grandTotal = 0;
	vector< int > sizes;
	for ( int i = 0; i < volumeList.size(); i++ ){
		nbfMrcReader reader;
		reader.setFileName( volumeList[i].getFileName().c_str() );
		grandTotal += reader.getDims()[2];
		sizes.push_back( reader.getDims()[2] );
	}
	cout << "Total number of images = " << grandTotal << endl;
		
	// reset array
	R.free();

	// switch to uni-processor
	if ( num_procs == 1 ){

		nbfFourierImageMetric< Pixel, 3 > fmetric( this->imageFilter, this->fourierFilter );
		int globalIndex = 0;
		for ( int i = 0; i < volumeList.size(); i++ ){
			fmetric.setInput1( &volumeList[i] );
			for ( int j = 0; j < sizes[i]; j++ ){
				Array< complex< Pixel >, 1 > C;
				Array< Pixel, 2 > M;
				fmetric.get2DimensionRepresentationHalf( C, M, j );
				if ( R.size() == 0 ){
					cerr << "Initializing representation size to " << grandTotal << " x " << C.rows() << " x 1 ... ";
					R.resize( grandTotal, C.rows(), 1 );
					cerr << "done." << endl;
				}
				R( globalIndex++, Range::all(), 0 ) = real( C );
			}
		}
	} else {

		// multiprocessor
		int * params = new int[ this->parameterSize ];

		int volumesPerProcess = ceil( grandTotal / ( num_procs - 1.0 ) );

		// start all slave processes
		for ( int i = 1; i < num_procs; i++ ){

			// compute volume range for this process
			int startIndex = ( i - 1 ) * volumesPerProcess;
			int endIndex = blitz::extrema::min( grandTotal - 1, i * volumesPerProcess - 1 );

			params[0] = inputFile.str().size();
			params[1] = startIndex;
			params[2] = endIndex;
			params[3] = useRealRepresentation;
			params[4] = this->overlapNormalizedDistances;
			params[5] = 0; // this->useMutualCorrelation;
			//params[5] = this->useMutualCorrelation;
			params[6] = bin_factor;
			params[7] = fold;

			// send data size and tag information
			MPI_Send ( params, this->parameterSize, MPI_INT, i, this->tag_representation_2D, MPI_COMM_WORLD );

			// send object data and object type information
			MPI_Send ( (char*)inputFile.str().c_str(), params[0], MPI_CHAR, i, this->tag_representation_2D, MPI_COMM_WORLD );

			if ( endIndex == volumeList.size() - 1 ){
				num_procs = i + 1;
				break;
			}
		}

		int jobsDone = 0;

		// loop while not all jobs completed
		while ( jobsDone < num_procs - 1 ){

			MPI_Status status;

			// receive notification from slaves
			MPI_Recv( params, this->parameterSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			int source = status.MPI_SOURCE;

			// we can now read the result from disk
			stringstream currentFile;
			if ( source < 10 ){
				currentFile << "mpi_00" << source << ".tmp";
			} else if ( source < 100 ){
				currentFile << "mpi_0" << source << ".tmp";
			} else {
				currentFile << "mpi_" << source << ".tmp";
			}
			nbfMatlabReader r;
			r.setFileName( currentFile.str().c_str() );
			Array< Pixel, 3 > partialR;
			r.read( partialR );

			if ( R.size() == 0 ){
				// R.resize( volumeList.size(), partialR.cols(), partialR.depth() );
				R.resize( grandTotal, partialR.cols(), partialR.depth() );
			}
			R( Range( params[1], params[2] ), Range::all(), Range::all() ) = partialR;

			char command[400];
			sprintf( command, "rm %s", currentFile.str().c_str() );
			system( command );

			// int currentTotal = 0;
			// for ( int i = 0; i < params[1]; i++ ){
				// nbfMrcReader reader;
				// reader.setFileName( volumeList[i].getFileName().c_str() );
				// currentTotal += reader.getDims()[2];
			// }

			// R( Range( currentTotal, currentTotal + partialR.rows() - 1 ), Range::all(), Range::all() ) = partialR;

			jobsDone++;
		}
		delete [] params;
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: get2DAverages( stringstream & inputFile, Array< Pixel, 3 > & classes, Array< Pixel, 3 > & averages )
{
	int num_procs;
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	vector< nbfWedgedSubImage3D< Pixel > > volumeList;
	nbfWedgedSubImage3D< Pixel > :: read( inputFile.str(), volumeList );

	// reset array
	averages.free();

	// switch to uni-processor
	if ( num_procs == 1 ){

		// turn all filters off to generate averages
		this->fourierFilter->bandPassOff();
		this->fourierFilter->bfactorOff();
		this->imageFilter->paddingOn(2);
		nbfFourierImageMetric< Pixel, 3 > fmetric( this->imageFilter, this->fourierFilter );
		
		vector< int > sizes;
		for ( int i = 0; i < volumeList.size(); i++ ){
			nbfMrcReader reader;
			reader.setFileName( volumeList[i].getFileName().c_str() );
			sizes.push_back( reader.getDims()[2] );
		}

		for ( int i = 0; i < classes.rows(); i++ ){
			for ( int j = 0; j < classes.cols(); j++ ){
				if ( classes( i, j ) > 0 ){
					Array< Pixel, 2 > M;
					Array< complex< Pixel >, 1 > C;
					int stackIndex = 0;
					int lastStackCounter = 0;
					int stackCounter = sizes[stackIndex];
					while ( j >= stackCounter ){
						stackIndex++;
						lastStackCounter = stackCounter;
						stackCounter += sizes[stackIndex];
					}
					cout << j << ", Volume stack " << stackIndex << ", index in stack " << j-lastStackCounter << endl;
					fmetric.setInput1( &volumeList[stackIndex] );
					fmetric.get2DimensionRepresentationHalf( C, M, j-lastStackCounter, false );
					if ( averages.size() == 0 ){
						cout << "Initializing averages size to " << M.rows() << " x " << M.cols() << " x " << classes.rows();
						averages.resize( M.rows(), M.cols(), classes.rows() );
						averages = 0;
					}
					averages( Range :: all(), Range :: all(), i ) += M;
				}
			}
		}
	
	} else {

		// multiprocessor
		int * params = new int[ this->parameterSize ];

		// dump classes so that slaves can see it
		nbfMatlabWriter w;
		w.setFileName("classes.bin");
		w.write(classes);

		int jobsSubmitted = 0;
		
		// start all slave processes
		for ( int i = 1; i < num_procs; i++ ){

			params[0] = inputFile.str().size();
			params[1] = jobsSubmitted;

			// send data size and tag information
			MPI_Send ( params, this->parameterSize, MPI_INT, i, this->tag_average_2D, MPI_COMM_WORLD );

			// send object data and object type information
			MPI_Send ( (char*)inputFile.str().c_str(), params[0], MPI_CHAR, i, this->tag_average_2D, MPI_COMM_WORLD );

			// cout << "Class " << jobsSubmitted << " submitted to node "<< i << endl;
			
			jobsSubmitted++;
			
			if ( jobsSubmitted == classes.rows() ){
				// num_procs = i + 1;
				break;
			}
		}

		int jobsDone = 0;

		// loop while not all jobs completed
		while ( jobsDone < classes.rows() ){

			MPI_Status status;

			// receive notification from slaves
			MPI_Recv( params, this->parameterSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			int source = status.MPI_SOURCE;

			// cout << "Received class " << params[1] << " from node "<< source << endl;
	
			// we can now read the result from disk
			stringstream currentFile;
			if ( source < 10 ){
				currentFile << "mpi_00" << source << ".tmp";
			} else if ( source < 100 ){
				currentFile << "mpi_0" << source << ".tmp";
			} else {
				currentFile << "mpi_" << source << ".tmp";
			}
			nbfMatlabReader r;
			r.setFileName( currentFile.str().c_str() );
			Array< Pixel, 2 > partialR;
			r.read( partialR );

			// delete tmp file
			char command[200];
			sprintf( command, "rm %s", currentFile.str().c_str() );
			system( command );

			if ( averages.size() == 0 ){
				averages.resize( partialR.rows(), partialR.cols(), classes.rows() );
			}

			averages( Range :: all(), Range :: all(), params[1] ) = partialR;

			jobsDone++;
			
			if ( jobsSubmitted < classes.rows() ){
				params[0] = inputFile.str().size();
				params[1] = jobsSubmitted;

				// send data size and tag information
				MPI_Send ( params, this->parameterSize, MPI_INT, source, this->tag_average_2D, MPI_COMM_WORLD );

				// send object data and object type information
				MPI_Send ( (char*)inputFile.str().c_str(), params[0], MPI_CHAR, source, this->tag_average_2D, MPI_COMM_WORLD );
				
				// cout << "Class " << jobsSubmitted << " submitted to node "<< source << endl;
				
				jobsSubmitted++;
			}			
			
		}
		
		// delete tmp file
		char command[200];
		sprintf( command, "rm classes.bin" );
		system( command );
		
		delete [] params;
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: getDistances( stringstream & representationFile, Array< Pixel, 2 > & D )
{
	int num_procs;
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	nbfMatlabReader r;
	r.setFileName( representationFile.str().c_str() );
	Array< Pixel, 3 > R;
	r.read(R);

	// reset array
	D.resize( R.rows(), R.rows() );
	D = 0;

	// switch to uni-processor
	if ( num_procs == 1 ){

		for ( int i = 0; i < D.rows(); i++ ){
			for ( int j = i + 1; j < D.rows(); j++ ){

				if ( R.depth() == 1 ){
					
					// Using PCA dimensionality reduction
	
					// L2-norm
					//D(i,j) = sqrt( mean( pow2( R( i, Range::all(), 0 ) - R( j, Range::all(), 0 ) ) ) );
					//D(i,j) = sqrt( sum( pow2( R( i, Range::all(), 0 ) - R( j, Range::all(), 0 ) ) ) );
					D(i,j) = sqrt( sum( pow2( R( i, Range(1,toEnd), 0 ) - R( j, Range(1,toEnd), 0 ) ) ) );

				} else {

					// Using Fourier space metric
					
					Array< complex< Pixel >, 1 > c1( R.cols() );
					real( c1 ) = R( i, Range::all(), 0 );
					imag( c1 ) = R( i, Range::all(), 1 );

					Array< complex< Pixel >, 1 > c2( R.cols() );
					real( c2 ) = R( j, Range::all(), 0 );
					imag( c2 ) = R( j, Range::all(), 1 );

					Array< Pixel, 1 > wcombined1( R.cols() );
					wcombined1 = R( i, Range::all(), 2 );

					Array< Pixel, 1 > wcombined2( R.cols() );
					wcombined2 = R( j, Range::all(), 2 );

					Pixel overlap = sum( wcombined1 * wcombined2 );

					c1 = c1 * wcombined1 * wcombined2;
					c2 = c2 * wcombined1 * wcombined2;

					//// SKIP NORMALIZATION
					//c1 /= sqrtf( mean( real( c1 * conj(c1) ) ) );
					//c2 /= sqrtf( mean( real( c2 * conj(c2) ) ) );

					c1 = c1 - c2;

					if ( ( this->overlapNormalizedDistances == true ) && ( overlap > 0 ) ){
						// L1 norm
						D(i,j) = sum( real( c1 * conj(c1) ) ) / overlap;
					} else {
						D(i,j) = mean( real( c1 * conj(c1) ) );
					}
				}

				// symmetrize
				D(j,i) = D(i,j);
			}
		}

	} else {

		// multiprocessor
		int * params = new int[ this->parameterSize ];

		int volumesPerProcess = blitz::extrema::max( 1, ceil( D.rows() * ( D.rows() - 1 ) / 2 / ( num_procs - 1.0 ) ) );

		// start all slave processes
		for ( int i = 1; i < num_procs; i++ ){

			// compute volume range for this process
			int startIndex = ( i - 1 ) * volumesPerProcess;
			int endIndex = blitz::extrema::min( D.rows() * ( D.rows() - 1 ) / 2 - 1, i * volumesPerProcess - 1 );

			params[0] = representationFile.str().size();
			params[1] = startIndex;
			params[2] = endIndex;
			params[3] = this->overlapNormalizedDistances;

			// send data size and tag information
			MPI_Send ( params, this->parameterSize, MPI_INT, i, this->tag_representation_distances, MPI_COMM_WORLD );

			// send object data and object type information
			MPI_Send ( (char*)representationFile.str().c_str(), params[0], MPI_CHAR, i, this->tag_representation_distances, MPI_COMM_WORLD );

			if ( endIndex == ( D.rows() * ( D.rows() - 1 ) / 2 - 1 ) ){
				num_procs = i + 1;
				break;
			}
		}

		int jobsDone = 0;

		// loop while not all jobs completed
		while ( jobsDone < num_procs - 1 ){

			MPI_Status status;

			// receive notification from slaves
			MPI_Recv( params, this->parameterSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			int source = status.MPI_SOURCE;

			// we can now read the result from disk
			stringstream currentFile;
			if ( source < 10 ){
				currentFile << "mpi_00" << source << ".tmp";
			} else if ( source < 100 ){
				currentFile << "mpi_0" << source << ".tmp";
			} else {
				currentFile << "mpi_" << source << ".tmp";
			}
			nbfMatlabReader r;
			r.setFileName( currentFile.str().c_str() );
			Array< Pixel, 2 > partialD;
			r.read( partialD );

			int index = 0; int count = 0;
			for ( int i = 0; i < D.rows(); i++ ){
				for ( int j = i + 1; j < D.rows(); j++ ){
					if ( ( index >= params[1] ) && ( index <= params[2] ) ){
						D(i,j) = partialD( count, 0 );
						D(j,i) = D(i,j);
						count++;
					}
					index++;
				}
			}
			
			jobsDone++;
		}
		delete [] params;
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: getDistances( vector< nbfWedgedImage3D< Pixel > * > & list1, 
													          vector< nbfWedgedImage3D< Pixel > * > & list2, 
															  vector< TinyVector< int, 2 > > & indexesVector,
															  Array< Pixel, 3 > & D, int typeOfMetric )
{
	//// retrieve temporal alignment matrix from file (if any)
	//stringstream fileName1;
	//fileName1 << "tmp.distances.matlab";

	//nbfMatlabReader mreader;
	//mreader.setFileName( fileName1.str().c_str() );

	//Array< Pixel, 3 > tmpD;
	//mreader.read( tmpD );
	//if ( sum( abs( D.shape() - tmpD.shape() ) ) == 0 ){
	//	D = where( tmpD != -1, tmpD, D );
	//}

	int jobNotSubmitted = 0;
	int jobSubmitted = 1;
	int jobDone = 2;

	Array< int, 1 > queueStatus( indexesVector.size() );
	queueStatus = jobNotSubmitted;

	for ( int i = 0; i < indexesVector.size(); i++ ){
		if ( D(i,0,0) != -1 ){
			queueStatus(i) = jobDone;
		}
	}

	int distancesToCompute = sum( where( queueStatus == jobNotSubmitted, 1, 0 ) );
	cout << "Total alignments to compute = " << distancesToCompute << endl;

	int num_procs;
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	// switch to uni-processor
	if ( num_procs == 1 ){

		nbfCorrelationImageMetric< Pixel, 3 > * currentMetric;
		nbfFourierImageMetric< Pixel, 3 > refinement( this->imageFilter, this->fourierFilter );
		nbfFourierImageMetricCore< Pixel, 3 > core( this->imageFilter, this->fourierFilter );

		if ( typeOfMetric == 0 ){
			currentMetric = this;
			this->refinementOn();
		} else if ( typeOfMetric == 1 ){
			currentMetric = &refinement;
		} else if ( typeOfMetric == 2 ){
			currentMetric = &core;
		} else if ( typeOfMetric == 3 ){
			currentMetric = this;
			this->refinementOff();
		} else {
			cerr << "ERROR - Must specify valid type of metric." << endl;
		}

		// make sure we are using same constraints as this metric
		if ( ( typeOfMetric == 1 ) || ( typeOfMetric == 2 ) ){
			currentMetric->setRotationSearchRestriction( this->restrictRotationSearch );
			currentMetric->setTranslationSearchRestriction( this->restrictTranslationSearch );
			currentMetric->setToComputeOverlapNormalizedDistances( this->overlapNormalizedDistances );
			currentMetric->setToUseMutualCorrelation( this->useMutualCorrelation );
			currentMetric->setMissingWedgeCompensation( this->useMissingWedgeCompensation );
			//currentMetric->imageFilter->maskOn( this->imageFilter->maskFile );
		}

		int currentx = -1;
		int currenty = -1;

		for ( int i = 0; i < indexesVector.size(); i++ ){

			// skip if already computed
			if ( D( i, 0, 0 ) == -1 ){
			
				int x = indexesVector[i][0];
				int y = indexesVector[i][1];

				if ( x != currentx ){
					currentMetric->setInput1( list1[ x ] );
				}
				if ( y != currenty ){
					currentMetric->setInput2( list2[ y ] );
					//stringstream f;
					//list2[ y ]->serialize(f);
					//cout << f.str().c_str() << endl;
				}

				currentx = x;
				currenty = y;

				int offset = 1;
				if ( currentMetric->getNumberOfCandidates() == 1 ){
					D( i, 0, 0 ) = currentMetric->getDistance();

					if ( D.cols() > 1 + 16 ){
						D( i, 1, 0 ) = currentMetric->getCorrelationScale();
						offset = 2;
					}
					if ( D.cols() > 2 + 16 ){
						D( i, 2, 0 ) = currentMetric->getWedgeOverlap( currentMetric->getTransform() );
						offset = 3;
					}
					double matrix[16];
					vtkMatrix4x4::DeepCopy( matrix, currentMetric->getTransform()->GetMatrix() );

					for ( int k = 0; k < 16; k++ ){
						D( i, k + offset, 0 ) = matrix[k];
					}

				} else {
                    currentMetric->getDistance();
					D( i, 0, Range::all() ) = numeric_limits< Pixel > :: max();
					for ( int candidate = 0; candidate < currentMetric->candidateTransforms.size(); candidate++ ){
						D( i, 0, candidate ) = currentMetric->getGivenDistance(candidate);
						if ( D.cols() > 1 + 16 ){
							D( i, 1, candidate ) = currentMetric->getGivenCorrelationScale(candidate);
							offset = 2;
						}
						if ( D.cols() > 2 + 16 ){
							D( i, 2, candidate ) = currentMetric->getGivenWedgeOverlap(candidate);
							offset = 3;
						}
						// set transformation matrix
						for ( int k = 0; k < 16; k++ ){
							D( i, k + offset, candidate ) = currentMetric->candidateTransforms[candidate](k);
						}
					}
				}

				//nbfMatlabWriter mw;
				//mw.setFileName( fileName1.str().c_str() );
				//mw.write( D );
			}
		}
	} else {

		// estimate bundle size as average number of jobs per process
		switch ( typeOfMetric ){
			case 0:
				this->MPIbundleNumber = 5;
				break;
			case 1:
				this->MPIbundleNumber = 5;
				break;
			case 2:
				this->MPIbundleNumber = 50;
				break;
			case 3:
				this->MPIbundleNumber = 5;
				break;
		}

		this->MPIbundleNumber = blitz::extrema::min( blitz::extrema::max( 1, ceil( 1.0 * distancesToCompute / ( num_procs - 1 ) ) ), this->MPIbundleNumber );

		cout << " Using bundle size = " << this->MPIbundleNumber << endl;

		// multiprocessor
		int * params = new int[ this->parameterSize ];

		// start all slave processes
		for ( int i = 1; i < num_procs; i++ ){
			if ( sum( where( queueStatus == jobNotSubmitted, 1, 0 ) ) > 0 ){

				stringstream objectData1;
				
				int bundleMaxNumber = blitz::extrema::min( this->MPIbundleNumber, sum( where( queueStatus == jobNotSubmitted, 1, 0 ) ) );
				
				vector< int > positionsToCompute;

				for ( int bundle = 0; bundle < bundleMaxNumber; bundle++ ){
					// search for next position to compute
					int nextBundleJob = 0;
					while ( nextBundleJob < queueStatus.size() ){
						if ( queueStatus(nextBundleJob) == jobNotSubmitted ){
							break;
						}
						nextBundleJob++;
					}

					// update queue status
					queueStatus( nextBundleJob ) = jobSubmitted;
					
					positionsToCompute.push_back( nextBundleJob );

					objectData1 << nextBundleJob << endl;

					bool skipIndex = false;
					// send first object only if different than last one sent
					if ( bundle > 0	 ){
						if ( indexesVector[ nextBundleJob ][0] == indexesVector[ positionsToCompute[ positionsToCompute.size() - 2 ] ][0] ){
							skipIndex = true;
						}
					}
					objectData1 << skipIndex << endl;
					if ( skipIndex == false ){
						list1[ indexesVector[ nextBundleJob ][0] ]->serialize( objectData1 );
					}
					// send second object
					list2[ indexesVector[ positionsToCompute[ bundle ] ][1] ]->serialize( objectData1 );
				}

				// parameters for MPI
				int objectSize = objectData1.str().size();
				params[0] = 1 + D.cols();
				params[1] = objectSize;
				params[2] = typeOfMetric;
				params[3] = bundleMaxNumber;
				params[4] = this->getTranslationSearchRestriction();
				params[5] = this->getRotationSearchRestriction();
				params[6] = this->overlapNormalizedDistances;
				params[7] = this->useMutualCorrelation;
				params[8] = this->useMissingWedgeCompensation;
				params[9] = list1[ 0 ]->getTypeId();
				params[10] = list2[ 0 ]->getTypeId();

				// save parameters to file				
				stringstream processFile;
				if ( i < 10 ){
					processFile << "mpi_00" << i << ".tmp";
				} else if ( i < 100 ){
						processFile << "mpi_0" << i << ".tmp";
				} else {
					processFile << "mpi_" << i << ".tmp";
				}
				std :: ofstream dst( processFile.str().c_str(), std :: ofstream :: binary );
				dst << objectData1.rdbuf();
				dst.close();

				//cout << "Parameters saved to file " << processFile.str().c_str() << ", size = " << objectData1.str().size() << endl;

				// send parameters
				MPI_Send ( params, this->parameterSize, MPI_INT, i, this->tag_processing, MPI_COMM_WORLD );

				//// send object data
				//MPI_Send ( (char*)objectData1.str().c_str(), objectSize, MPI_CHAR, i, list1[ 0 ]->getTypeId(), MPI_COMM_WORLD );
			}
		}

		// this->MPIbundleNumber = 1;

		int count = 1;

		int alignmentsDone = 0;

		// loop while not all jobs completed
		while ( sum( where( queueStatus != jobDone, 1, 0 ) ) > 0 ){

			MPI_Status status;

			// receive notification from slaves
			MPI_Recv( params, this->parameterSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			int source = status.MPI_SOURCE;

			// we can now read the result from disk
			stringstream processFile;
			if ( source < 10 ){
				processFile << "mpi_00" << source << ".tmp";
			} else if ( source < 100 ){
				processFile << "mpi_0" << source << ".tmp";
			} else {
				processFile << "mpi_" << source << ".tmp";
			}

			nbfMatlabReader r;
			r.setFileName( processFile.str().c_str() );
			Array< Pixel, 2 > processResults;
			r.read( processResults );

			int tag = status.MPI_TAG;

			// retrieve number of alignment candidates
			int alignmentCandidates = params[2];

			int resultIndex = 1;

			for ( int bundleIndex = 0; bundleIndex < processResults(0,0); bundleIndex++ ){
				int index = processResults( 0, resultIndex++ );
				for ( int candidates = 0; candidates < alignmentCandidates; candidates++ ){
					//for ( int z = 0; z < this->resultSize - 1; z++ ){
					for ( int z = 0; z < D.cols(); z++ ){
						D( index, z, candidates ) = processResults( 0, resultIndex++ );
					}
				}

				// update queue status
				queueStatus( index ) = jobDone;

				//// print progress and save intermediate results depending on the speed of computations
				//if ( ( index > 0 ) &&
				//	( ( ( typeOfMetric == 0 ) && ( index % 100 == 0 ) ) || // this takes very long to compute, so we save/print often
				//	( ( typeOfMetric == 1 ) && ( index % 100 == 0 ) ) || // this takes very long to compute, so we save/print often
				//	( ( typeOfMetric == 2 ) && ( index % 500 == 0 ) ) || // this is the quickest, so we save sporadically
				//	( ( typeOfMetric == 3 ) && ( index % 100 == 0 ) ) ) ){ // this is quick so we save more sporadically

				//		if ( D.cols() > 16 + 2 ){
				//			cerr << "D(" << index << " of " << queueStatus.size() << ") = " << D(index,0) << ", scale=" << D(index,1) << ", overlap=" << D(index,2) << endl;
				//		} else if ( D.cols() > 16 + 1 ){
				//			cerr << "D(" << index << " of " << queueStatus.size() << ") = " << D(index,0) << ", scale=" << D(index,1) << endl;
				//		} else if ( D.cols() > 16 ){
				//			cerr << "D(" << index << " of " << queueStatus.size() << ") = " << D(index,0) << endl;
				//		}

				//		nbfMatlabWriter mw;
				//		mw.setFileName( fileName1.str().c_str() );
				//		mw.write( D );
				//	}
			}

			// update counter of aligments done
			alignmentsDone += processResults(0,0);

			if ( ( alignmentsDone > 0 ) && ( alignmentsDone % 1000 == 0 ) ){
				cout << " " << alignmentsDone << " of " << distancesToCompute << " alignments completed." << endl;
			}
			//int alignmentsCompleted = sum( where( queueStatus == jobDone, 1, 0 ) );
			//if ( ( alignmentsCompleted > 0 ) && ( alignmentsCompleted >= count * 100 ) ){
			//	nbfMatlabWriter mw;
			//	mw.setFileName( fileName1.str().c_str() );
			//	mw.write( D );
			//	if ( alignmentsCompleted < 1000 ){
			//		cerr << " " << alignmentsCompleted << " of " << distancesToCompute << " completed." << endl;
			//	} else {
			//		cerr << alignmentsCompleted << " of " << distancesToCompute << " completed." << endl;
			//	}
			//	count++;
			//}

			// send remaining jobs if not already submitted
			if ( ( tag == this->tag_done ) && ( sum( where( queueStatus == jobNotSubmitted, 1, 0 ) ) > 0 ) ){

				stringstream objectData1;

				int bundleMaxNumber = blitz::extrema::min( this->MPIbundleNumber, sum( where( queueStatus == jobNotSubmitted, 1, 0 ) ) );
				
				vector< int > positionsToCompute;

				for ( int bundle = 0; bundle < bundleMaxNumber; bundle++ ){
					// search for next position to compute
					int nextBundleJob = 0;
					while ( nextBundleJob < queueStatus.size() ){
						if ( queueStatus(nextBundleJob) == jobNotSubmitted ){
							break;
						}
						nextBundleJob++;
					}

					// update queue status
					queueStatus( nextBundleJob ) = jobSubmitted;

					positionsToCompute.push_back( nextBundleJob );

					objectData1 << nextBundleJob << endl;

					bool skipIndex = false;
					// send first object only if different than last one sent
					if ( bundle > 0	 ){
						if ( indexesVector[ nextBundleJob ][0] == indexesVector[ positionsToCompute[ positionsToCompute.size() - 2 ] ][0] ){
							skipIndex = true;
						}
					}
					objectData1 << skipIndex << endl;
					if ( skipIndex == false ){
						list1[ indexesVector[ nextBundleJob ][0] ]->serialize( objectData1 );
					}
					// send second object
					list2[ indexesVector[ positionsToCompute[ bundle ] ][1] ]->serialize( objectData1 );
				}

				// parameters for MPI
				int objectSize = objectData1.str().size();
				params[0] = 1 + D.cols();
				params[1] = objectSize;
				params[2] = typeOfMetric;
				params[3] = bundleMaxNumber;
				params[4] = this->getTranslationSearchRestriction();
				params[5] = this->getRotationSearchRestriction();
				params[6] = this->overlapNormalizedDistances;
				params[7] = this->useMutualCorrelation;
				params[8] = this->useMissingWedgeCompensation;
				params[9] = list1[ 0 ]->getTypeId();
				params[10] = list2[ 0 ]->getTypeId();

				// save parameters to file				
				stringstream processFile;
				if ( source < 10 ){
					processFile << "mpi_00" << source << ".tmp";
				} else if ( source < 100 ){
					processFile << "mpi_0" << source << ".tmp";
				} else {
					processFile << "mpi_" << source << ".tmp";
				}
				std :: ofstream dst( processFile.str().c_str(), std :: ofstream :: binary );
				dst << objectData1.rdbuf();
				dst.close();

				//cout << "Parameters saved to file " << processFile.str().c_str() << ", size = " << objectData1.str().size() << endl;

				// send parameters
				MPI_Send ( params, this->parameterSize, MPI_INT, source, this->tag_processing, MPI_COMM_WORLD );

				// send object data
				//MPI_Send ( (char*)objectData1.str().c_str(), objectSize, MPI_CHAR, source, list1[ 0 ]->getTypeId(), MPI_COMM_WORLD );
			}		
		}
		delete [] params;
	}
	//tmpD.resize( D.shape() );
	//tmpD = -1;
	//nbfMatlabWriter mw;
	//mw.setFileName( fileName1.str().c_str() );
	//mw.write( tmpD );
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: makeDistanceMatrix( Array< Pixel, 4 > & D )
{
	for ( int i = 0; i < D.rows(); i++ ){
		for ( int j = i + 1; j < D.cols(); j++ ){

			// symetrize distance, ccc and overlap values
			D( j, i, Range(0,2), Range::all() ) = D( i, j, Range(0,2), Range::all() );

			for ( int p = 0; p < D.extent(fourthDim); p++ ){

				// assign inverse transformation to symmetric entry
				double matrix[16];
				for ( int k = 0; k < 16; k++ ){
					matrix[k] = D( i, j, k + 3, p );
				}

				double matrixInverse[16];
				vtkMatrix4x4::Invert( matrix, matrixInverse );

				for ( int k = 0; k < 16; k++ ){
					D( j, i, k + 3, p ) = matrixInverse[k];
				}
			}
		}
	}

	// set diagonal to 0
	for ( int i = 0; i < D.rows(); i++ ){
		D( i, i, Range::all(), Range::all() ) = 0;
		// set transform to identity
		double matrix[16];
		vtkMatrix4x4::Identity( matrix );
		for ( int k = 0; k < 16; k++ ){
			D( i, i, k + 3, Range::all() ) = matrix[k];
		}
	}
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: finalizeMPI()
{
	int num_procs;
	MPI_Comm_size( MPI_COMM_WORLD, &num_procs );

	int * voidParams = new int[ this->parameterSize ];

	// send termination message to all nodes
	for ( int i = 1; i < num_procs; i++ ){
		MPI_Send ( voidParams, this->parameterSize, MPI_INT, i, this->tag_done, MPI_COMM_WORLD );
	}
	delete [] voidParams;
}

template< class Pixel, const int Dim >
void nbfCorrelationImageMetric< Pixel, Dim > :: slaveMPI()
{
	int source = 0;
	int tag = this->tag_processing;

	char * objectData1;
	char * objectData2;

	MPI_Status status;

	int * params = new int[ this->parameterSize ];

	nbfWedgedImage3D< Pixel > * vol1 = NULL;
	nbfWedgedImage3D< Pixel > * vol2 = NULL;

	int objectType1, objectType2;

	//  Get this processes's rank.
	int my_id;
	MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

	// we save the result to disk
	stringstream processFile;
	if ( my_id < 10 ){
		processFile << "mpi_00" << my_id << ".tmp";
	} else if ( my_id < 100 ){
		processFile << "mpi_0" << my_id << ".tmp";
	} else {
		processFile << "mpi_" << my_id << ".tmp";
	}

	while ( tag != this->tag_done ){

		// retrieve object size and tag information
		MPI_Recv( params, this->parameterSize, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );		
		tag = status.MPI_TAG;

		// DISTANCE COMPUTATIONS
		if ( tag == this->tag_processing ){
			
			int newResultSize = params[0];
			int objectSize = params[1];
			int typeOfMetric = params[2];
			int bundleMaxNumber = params[3];
			float translationRestriction = params[4];
			float rotationRestriction = params[5];
			float overlapNormalizedDistances = params[6];
			float useMutualCorrelation = params[7];
			float useMissingWedgeCompensation = params[8];
			objectType1 = params[9];
			objectType2 = params[10];

			// save parameters to file				
			std :: ifstream src( processFile.str().c_str(), std :: ofstream :: binary );

			// allocate array
			objectData1 = new char[ objectSize ];

			//cout << "About to retieve file " << processFile.str().c_str() << endl;
			stringstream inputString1;
			src >> inputString1.rdbuf();
			src.close();
			//cout << "About to retieve file " << inputString1.str().size() << endl;

			//// retrieve object data
			//MPI_Recv( objectData1, objectSize, MPI_CHAR, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			//stringstream inputString1( objectData1 );

			// store results in array
			//Array< Pixel, 2 > bundleResults( 1, 1 + params[3] * ( 1 + this->getNumberOfCandidates() * ( this->resultSize - 1 ) ) );
			Array< Pixel, 2 > bundleResults( 1, 1 + bundleMaxNumber * ( 1 + this->getNumberOfCandidates() * ( newResultSize - 1 ) ) );
			
			// first store bundle size
			bundleResults( 0, 0 ) = bundleMaxNumber;

			int bundleIndex = 1;

			nbfCorrelationImageMetric< Pixel, 3 > * currentMetric;
			nbfFourierImageMetric< Pixel, 3 > refinement( this->imageFilter, this->fourierFilter );
			nbfFourierImageMetricCore< Pixel, 3 > core( this->imageFilter, this->fourierFilter );

			if ( typeOfMetric == 0 ){
				currentMetric = this;
				this->refinementOn();
			} else if ( typeOfMetric == 1 ){
				currentMetric = &refinement;
			} else if ( typeOfMetric == 2 ){
				currentMetric = &core;
			} else if ( typeOfMetric == 3 ){
				currentMetric = this;
				this->refinementOff();
			} else {
				cerr << "ERROR - Must specify valid type of metric. In " << __FILE__ << ", " << __LINE__ << endl;
			}

			// set metric parameters
			currentMetric->setRotationSearchRestriction( rotationRestriction );
			currentMetric->setTranslationSearchRestriction( translationRestriction );
			currentMetric->setToComputeOverlapNormalizedDistances( overlapNormalizedDistances > 0 );
			currentMetric->setToUseMutualCorrelation( useMutualCorrelation > 0 );
			currentMetric->setMissingWedgeCompensation( useMissingWedgeCompensation > 0 );

			int cindex;

			for ( int bundle = 0; bundle < bundleMaxNumber; bundle++ ){

				// retrieve index of distance to compute
				inputString1 >> cindex;

				// determine if using same input as before
				bool skipIndex;
				inputString1 >> skipIndex;

				// set input1
				if ( skipIndex == false ){
					if ( vol1 != NULL ){
						delete vol1;
						vol1 = NULL;
					}
					if ( objectType1 == NBF_WEDGED_SUB_IMAGE_3D ){
						vol1 = new nbfWedgedSubImage3D< Pixel >();
						vol1->unserialize( inputString1 );
						vtkTransform * t = vtkTransform :: New();
						reinterpret_cast< nbfWedgedSubImage3D< Pixel > * >(vol1)->getTransform(t);
					} else if ( objectType1 == NBF_WEDGED_AVERAGE_IMAGE_3D ){
						vol1 = new nbfWedgedAverageImage3D< Pixel >();
						vol1->unserialize( inputString1 );
					} else {
						vol1 = NULL;
						cerr << "ERROR: Unknown object type." << endl;
					}
					currentMetric->setInput1( vol1 );
				}

				// set input2
				if ( objectType2 == NBF_WEDGED_SUB_IMAGE_3D ){
					vol2 = new nbfWedgedSubImage3D< Pixel >();
					vol2->unserialize( inputString1 );
					vtkTransform * t = vtkTransform :: New();
					reinterpret_cast< nbfWedgedSubImage3D< Pixel > * >(vol2)->getTransform(t);
				} else if ( objectType2 == NBF_WEDGED_AVERAGE_IMAGE_3D ){
					vol2 = new nbfWedgedAverageImage3D< Pixel >();
					vol2->unserialize( inputString1 );
				} else {
					vol2 = NULL;
					cerr << "ERROR: Unknown object type." << endl;
				}
				currentMetric->setInput2( vol2 );

				currentMetric->getDistance();

				// compute size of result
				//params[0] = currentMetric->candidateCorrelationPeaks.size() * ( this->resultSize - 1 );
				params[0] = currentMetric->candidateCorrelationPeaks.size() * ( newResultSize - 1 );

				// set current volume index
				params[1] = cindex;
				
				// set number of candidate alignments
				params[2] = currentMetric->candidateCorrelationPeaks.size();

				int index = 0;
				
				// store current volume index in result array
				bundleResults(0,bundleIndex++) = cindex;
				
				for ( int current = 0; current < currentMetric->candidateCorrelationPeaks.size(); current++ ){
					bundleResults(0,bundleIndex++) = currentMetric->getGivenDistance(current);
					if ( newResultSize - 1 > 1 + 16 ){
						bundleResults(0,bundleIndex++) = currentMetric->getGivenCorrelationScale(current);
					}
					if ( newResultSize - 1 > 2 + 16 ){
						bundleResults(0,bundleIndex++) = currentMetric->getGivenWedgeOverlap(current);
					}
					for ( int i = 0; i < 16; i++ ){
						bundleResults(0,bundleIndex++) = currentMetric->candidateTransforms[current](i);
					}
				}

				// if last element, save results in file and signal master we have finished and are ready to receive more
				if ( bundle == ( bundleMaxNumber - 1 ) ){
					nbfMatlabWriter w;
					w.setFileName( processFile.str().c_str() );
					w.write( bundleResults );
					MPI_Send ( params, this->parameterSize, MPI_INT, 0, this->tag_done, MPI_COMM_WORLD );
				}

				if ( vol2 != NULL ){
					delete vol2;
					vol2 = NULL;
				}
			}

			if ( vol1 != NULL ){
				delete vol1;
				vol1 = NULL;
			}

			delete [] objectData1;

		} else if ( tag == this->tag_averaging ){

			// allocate array
			int objectSize = params[0];
			char * objectData = new char[ objectSize ];

			// retrieve object data
			MPI_Recv( objectData, objectSize, MPI_CHAR, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			stringstream inputString( objectData );

			nbfWedgedAverageImage3D< Pixel > currentAverage;
			//cout << " SLAVE: unserializing..." << endl;
			//cout << " SLAVE received size = " << objectSize << "\n" << inputString.str().c_str() << endl;
			currentAverage.unserialize( inputString );
			//cout << " SLAVE: done unserializing." << endl;

			// get partial average
			Array< double, 3 > realIm, imagIm, accumulatedFourier;
			currentAverage.getImageMPI( realIm, imagIm, accumulatedFourier );

			// we save the result to disk
			stringstream currentFile;

			if ( my_id < 10 ){
				currentFile << "mpi_00" << my_id << ".tmp";
			} else if( my_id < 100 ){
				currentFile << "mpi_0" << my_id << ".tmp";
			} else {
				currentFile << "mpi_" << my_id << ".tmp";
			}
			nbfMatlabWriter w;
			w.setFileName( currentFile.str().c_str() );

#ifdef NBF_AVERAGE_IN_RECIPROCAL_SPACE
			Array< double, 4 > C( realIm.rows(), realIm.cols(), realIm.depth(), 3 );
			C( Range::all(), Range::all(), Range::all(), 0 ) = realIm;
			C( Range::all(), Range::all(), Range::all(), 1 ) = imagIm;
			C( Range::all(), Range::all(), Range::all(), 2 ) = accumulatedFourier;
			w.write(C);
			//cout << "SLAVE: C.shape() = " << C.shape() << endl;
#else
			w.write( realIm );
#endif
			// acknowledge master we've finished
			MPI_Send ( params, this->parameterSize, MPI_INT, 0, this->tag_averaging, MPI_COMM_WORLD );

			delete [] objectData;

		} else if ( tag == this->tag_variance ){

			// allocate array
			int objectSize = params[0];
			char * objectData = new char[ objectSize ];

			// retrieve object data
			MPI_Recv( objectData, objectSize, MPI_CHAR, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			stringstream inputString( objectData );

			nbfWedgedAverageImage3D< Pixel > currentAverage;
			currentAverage.unserialize( inputString );

			// get partial average
			Array< double, 3 > varImage;
			currentAverage.getImageVarianceMPI( varImage );

			// we save the result to disk
			stringstream currentFile;

			if ( my_id < 10 ){
				currentFile << "mpi_00" << my_id << ".tmp";
			} else if ( my_id < 100 ){
				currentFile << "mpi_0" << my_id << ".tmp";
			} else {
				currentFile << "mpi_" << my_id << ".tmp";
			}
			nbfMatlabWriter w;
			w.setFileName( currentFile.str().c_str() );

			w.write( varImage );

			// aknowledge master we've finished
			MPI_Send ( params, this->parameterSize, MPI_INT, 0, this->tag_variance, MPI_COMM_WORLD );

			delete [] objectData;

		} else if ( tag == this->tag_update_accum_wedge ){

			// allocate array
			int objectSize = params[0];
			char * objectData = new char[ objectSize ];

			// retrieve object data
			MPI_Recv( objectData, objectSize, MPI_CHAR, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			stringstream inputString( objectData );

			nbfWedgedAverageImage3D< Pixel > currentAverage;
			currentAverage.unserialize( inputString );

			// get partial wedge
			currentAverage.updateAccumulatedWedgeImage();

			// we save the result to disk
			stringstream currentFile;

			if ( my_id < 10 ){
				currentFile << "mpi_00" << my_id << ".tmp";
			} else if( my_id < 100 ){
				currentFile << "mpi_0" << my_id << ".tmp";
			} else {
				currentFile << "mpi_" << my_id << ".tmp";
			}
			nbfMatlabWriter w;
			w.setFileName( currentFile.str().c_str() );

			w.write( currentAverage.wedgeAccum );

			// acknowledge master we've finished
			MPI_Send ( params, this->parameterSize, MPI_INT, 0, this->tag_averaging, MPI_COMM_WORLD );

			delete [] objectData;

		} else if ( tag == this->tag_update_accum_spherical_wedge ){

			// allocate array
			int objectSize = params[0];
			char * objectData = new char[ objectSize ];

			// retrieve size parameters
			TinyVector< int, 2 > size( params[1], params[2] );

			// retrieve object data
			MPI_Recv( objectData, objectSize, MPI_CHAR, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			stringstream inputString( objectData );

			nbfWedgedAverageImage3D< Pixel > currentAverage;
			currentAverage.unserialize( inputString );

			// get partial cumulative wedge
			currentAverage.updateAccumulatedSphericalWedgeImageMPI( size, reinterpret_cast<nbfProjectionRotationMetric3D<Pixel>*>(this) );

			// we save the result to disk
			stringstream currentFile;

			if ( my_id < 10 ){
				currentFile << "mpi_00" << my_id << ".tmp";
			} else if( my_id < 100 ){
				currentFile << "mpi_0" << my_id << ".tmp";
			} else {
				currentFile << "mpi_" << my_id << ".tmp";
			}
			nbfMatlabWriter w;
			w.setFileName( currentFile.str().c_str() );

			w.write( currentAverage.sphericalWedgeAccum );

			// acknowledge master we've finished
			MPI_Send ( params, this->parameterSize, MPI_INT, 0, this->tag_averaging, MPI_COMM_WORLD );

			delete [] objectData;

		} else if ( tag == this->tag_representation ){

			// allocate array
			int objectDataSize = params[0];
			char * objectData = new char[ objectDataSize ];

			// retrieve object data
			MPI_Recv( objectData, objectDataSize, MPI_CHAR, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			stringstream inputFile( objectData );
		
			vector< nbfWedgedSubImage3D< Pixel > > volumeList;

			char * filename = new char[ objectDataSize + 1 ];
			inputFile.get( filename, objectDataSize + 1 );
			
			nbfWedgedSubImage3D< Pixel > :: read( filename, volumeList );

			Array< Pixel, 3 > R;
			
			nbfFourierImageMetric< Pixel, 3 > fmetric( this->imageFilter, this->fourierFilter );
			fmetric.setToComputeOverlapNormalizedDistances( params[4] > 0 );
			fmetric.setToUseMutualCorrelation( params[5] > 0 );
			
			int index = 0;
			
			int useRealRepresentation = params[3];

			for ( int i = params[1]; i <= params[2]; i++ ){
				fmetric.setInput1( &volumeList[i] );
				if ( useRealRepresentation > 0 ){
					Array< Pixel, 1 > C;
					fmetric.getLowDimensionRepresentationReal( C, params[6], params[7] );
					if ( R.size() == 0 ){
						if ( my_id == 1 ) cerr << "Initializing representation size to " << volumeList.size() << " x " << C.rows() << " ... ";
						R.resize( params[2] - params[1] + 1, C.rows(), 1 );
						if ( my_id == 1 ) cerr << "done." << endl;
					}
					R( index, Range::all(), 0 ) = C;
				} else {
					Array< complex< Pixel >, 1 > C;
					Array< Pixel, 1 > W;
					fmetric.getLowDimensionRepresentationHalf( C, W, params[7] );
					if ( R.size() == 0 ){
						if ( my_id == 1 ) cerr << "Initializing representation size to " << volumeList.size() << " x " << C.rows() << " x " << 3 << " ... ";
						R.resize( params[2] - params[1] + 1, C.rows(), 3 );
						if ( my_id == 1 ) cerr << "done." << endl;
					}
					R( index, Range::all(), 0 ) = real(C);
					R( index, Range::all(), 1 ) = imag(C);
					R( index, Range::all(), 2 ) = W;
				}
				index++;
			}

			// we save the result to disk
			stringstream currentFile;

			if ( my_id < 10 ){
				currentFile << "mpi_00" << my_id << ".tmp";
			} else if ( my_id < 100 ){
				currentFile << "mpi_0" << my_id << ".tmp";
			} else {
				currentFile << "mpi_" << my_id << ".tmp";
			}
			nbfMatlabWriter w;
			w.setFileName( currentFile.str().c_str() );
			w.write( R );

			delete [] filename;

			// send params unchanged (only contains starting end ending indexes)
			MPI_Send ( params, this->parameterSize, MPI_INT, 0, this->tag_representation, MPI_COMM_WORLD );
			
			delete [] objectData;

		} else if ( tag == this->tag_representation_2D ){

			// allocate array
			int objectDataSize = params[0];
			char * objectData = new char[ objectDataSize ];

			// retrieve object data
			MPI_Recv( objectData, objectDataSize, MPI_CHAR, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			stringstream inputFile( objectData );
		
			vector< nbfWedgedSubImage3D< Pixel > > volumeList;

			char * filename = new char[ objectDataSize + 1 ];
			inputFile.get( filename, objectDataSize + 1 );
	
			nbfWedgedSubImage3D< Pixel > :: read( filename, volumeList );

			Array< Pixel, 3 > R;
			
			nbfFourierImageMetric< Pixel, 3 > fmetric( this->imageFilter, this->fourierFilter );
			fmetric.setToComputeOverlapNormalizedDistances( params[4] > 0 );
			fmetric.setToUseMutualCorrelation( params[5] > 0 );
			
			int useRealRepresentation = params[3];

			int grandTotal = params[2] - params[1] + 1;
			// for ( int i = params[1]; i <= params[2]; i++ ){
				// nbfMrcReader reader;
				// reader.setFileName( volumeList[i].getFileName().c_str() );
				// grandTotal += reader.getDims()[2];
			// }
			
			vector< int > sizes;
			for ( int i = 0; i < volumeList.size(); i++ ){
				nbfMrcReader reader;
				reader.setFileName( volumeList[i].getFileName().c_str() );
				sizes.push_back( reader.getDims()[2] );
			}

			int index = 0;
			for ( int i = params[1]; i <= params[2]; i++ ){
				// determine series index
				int stackIndex = 0;
				int lastStackCounter = 0;
				int stackCounter = sizes[stackIndex];
				while ( i >= stackCounter ){
					stackIndex++;
					lastStackCounter = stackCounter;
					stackCounter += sizes[stackIndex];
				}
				fmetric.setInput1( &volumeList[stackIndex] );
				Array< complex< Pixel >, 1 > C;
				Array< Pixel, 2 > M;
				fmetric.get2DimensionRepresentationHalf( C, M, i-lastStackCounter );
				if ( R.size() == 0 ){
					if ( my_id == 1 ) cerr << "Initializing representation size to " << grandTotal << " x " << C.rows() << " x " << 1 << " ... ";
					R.resize( grandTotal, C.rows(), 1 );
					if ( my_id == 1 ) cerr << "done." << endl;
				}
				R( index++, Range::all(), 0 ) = real( C );
			}

			// we save the result to disk
			stringstream currentFile;

			if ( my_id < 10 ){
				currentFile << "mpi_00" << my_id << ".tmp";
			} else if ( my_id < 100 ){
				currentFile << "mpi_0" << my_id << ".tmp";
			} else {
				currentFile << "mpi_" << my_id << ".tmp";
			}
			nbfMatlabWriter w;
			w.setFileName( currentFile.str().c_str() );
			w.write( R );

			delete [] filename;

			// send params unchanged (only contains starting end ending indexes)
			MPI_Send ( params, this->parameterSize, MPI_INT, 0, this->tag_representation_2D, MPI_COMM_WORLD );
			
			delete [] objectData;

		} else if ( tag == this->tag_average_2D ){

			// allocate array
			int classIndex = params[1];
			
			// cout << "Slave " << my_id << " processing class " << classIndex << endl;
			
			int objectDataSize = params[0];
			char * objectData = new char[ objectDataSize ];

			// retrieve object data
			MPI_Recv( objectData, objectDataSize, MPI_CHAR, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			stringstream inputFile( objectData );
		
			vector< nbfWedgedSubImage3D< Pixel > > volumeList;

			char * filename = new char[ objectDataSize + 1 ];
			inputFile.get( filename, objectDataSize + 1 );

			nbfWedgedSubImage3D< Pixel > :: read( filename, volumeList );

			Array< Pixel, 3 > classes;
			nbfMatlabReader reader;
			reader.setFileName( "classes.bin" );
			reader.read(classes);

			Array< Pixel, 2 > R;
			
			this->fourierFilter->bandPassOff();
			this->fourierFilter->bfactorOff();
			this->imageFilter->paddingOn(2);
			
			nbfFourierImageMetric< Pixel, 3 > fmetric( this->imageFilter, this->fourierFilter );
			// fmetric.setToComputeOverlapNormalizedDistances( params[4] > 0 );
			// fmetric.setToUseMutualCorrelation( params[5] > 0 );
			
			vector< int > sizes;
			for ( int i = 0; i < volumeList.size(); i++ ){
				nbfMrcReader reader;
				reader.setFileName( volumeList[i].getFileName().c_str() );
				sizes.push_back( reader.getDims()[2] );
			}

			for ( int j = 0; j < classes.cols(); j++ ){
				if ( classes( classIndex, j ) > 0 ){
					Array< Pixel, 2 > M;
					Array< complex< Pixel >, 1 > C;
					int stackIndex = 0;
					int lastStackCounter = 0;
					int stackCounter = sizes[stackIndex];
					while ( j >= stackCounter ){
						stackIndex++;
						lastStackCounter = stackCounter;
						stackCounter += sizes[stackIndex];
					}
					// cout << "Slave " << my_id << " processing class " << classIndex << ", index " << j << ", Volume stack " << stackIndex << ", index in stack " << j-lastStackCounter << endl;
					fmetric.setInput1( &volumeList[stackIndex] );
					fmetric.get2DimensionRepresentationHalf( C, M, j-lastStackCounter, false );
					if ( R.size() == 0 ){
						R.resize( M.shape() );
						R = 0;
					}
					R += M;
				}
			}

			// we save the result to disk
			stringstream currentFile;

			if ( my_id < 10 ){
				currentFile << "mpi_00" << my_id << ".tmp";
			} else if ( my_id < 100 ){
				currentFile << "mpi_0" << my_id << ".tmp";
			} else {
				currentFile << "mpi_" << my_id << ".tmp";
			}
			nbfMatlabWriter w;
			w.setFileName( currentFile.str().c_str() );
			w.write( R );

			delete [] filename;

			// send params unchanged (only contains starting end ending indexes)
			MPI_Send ( params, this->parameterSize, MPI_INT, 0, this->tag_representation_2D, MPI_COMM_WORLD );
			
			delete [] objectData;

		} else if ( tag == this->tag_representation_distances ){
			
			// allocate array
			int objectDataSize = params[0];
			char * objectData = new char[ objectDataSize ];

			// retrieve object data
			MPI_Recv( objectData, objectDataSize, MPI_CHAR, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			stringstream inputFile( objectData );
		
			char * filename = new char[ objectDataSize + 1 ];
			inputFile.get( filename, objectDataSize + 1 );
			
			nbfMatlabReader r;
			r.setFileName( filename );
			Array< Pixel, 3 > R;
			r.read(R);

			// store results in linear array
			Array< Pixel, 2 > D( params[2] - params[1] + 1, 1 );

			//  Get this processes's rank.
			int my_id;
			MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

			int index = 0; int count = 0;
			for ( int i = 0; i < R.rows(); i++ ){
				for ( int j = i + 1; j < R.rows(); j++ ){
					if ( ( index >= params[1] ) && ( index <= params[2] ) ){

						if ( R.depth() == 1 ){
							
							// Using PCA dimensionality reduction
							
							//D(count,0) = sqrt( mean( pow2( R( i, Range::all(), 0 ) - R( j, Range::all(), 0 ) ) ) );
							D(count,0) = sqrt( sum( pow2( R( i, Range::all(), 0 ) - R( j, Range::all(), 0 ) ) ) );
						} else {

							// Using Fourier space representation

							Array< complex< Pixel >, 1 > c1( R.cols() );
							real( c1 ) = R( i, Range::all(), 0 );
							imag( c1 ) = R( i, Range::all(), 1 );

							Array< complex< Pixel >, 1 > c2( R.cols() );
							real( c2 ) = R( j, Range::all(), 0 );
							imag( c2 ) = R( j, Range::all(), 1 );

							Array< Pixel, 1 > wcombined1( R.cols() );
							wcombined1 = R( i, Range::all(), 2 );

							Array< Pixel, 1 > wcombined2( R.cols() );
							wcombined2 = R( j, Range::all(), 2 );

							Pixel overlap = sum( wcombined1 * wcombined2 );

							c1 = c1 * wcombined1 * wcombined2;
							c2 = c2 * wcombined1 * wcombined2;

							//c1 /= sqrtf( mean( real( c1 * conj(c1) ) ) );
							//c2 /= sqrtf( mean( real( c2 * conj(c2) ) ) );

							c1 = c1 - c2;

							if ( ( params[3] > 0 ) && ( overlap > 0 ) ){
								// L2 norm
								D(count,0) = sum( real( c1 * conj(c1) ) ) / overlap / overlap;
							} else {
								// NO NORMALIZATION
								D(count,0) = mean( real( c1 * conj(c1) ) );
							}
						}
						count++;
					}
					index++;
				}
			}

			// we save the result to disk
			stringstream currentFile;

			if ( my_id < 10 ){
				currentFile << "mpi_00" << my_id << ".tmp";
			} else if ( my_id < 100 ){
				currentFile << "mpi_0" << my_id << ".tmp";
			} else {
				currentFile << "mpi_" << my_id << ".tmp";
			}
			nbfMatlabWriter w;
			w.setFileName( currentFile.str().c_str() );
			w.write( D );

			delete [] filename;

			// send params unchanged (only contains starting end ending indexes)
			MPI_Send ( params, this->parameterSize, MPI_INT, 0, this->tag_representation, MPI_COMM_WORLD );
			
			delete [] objectData;
		}
	}
	delete [] params;
}
