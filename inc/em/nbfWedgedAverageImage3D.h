#pragma once

#include <vtkImageData.h>
#include <vtkTransform.h>
#include <vtkImageResample.h>
#include <vtkImageChangeInformation.h>

#include <io/nbfMrcReader.h>
#include <em/nbfImageFilter.h>
#include <em/nbfWedgedImage3D.h>

template< class Pixel > class nbfProjectionRotationMetric3D;
template< class Pixel > class nbfWedgedImage3D;
template< class Pixel > class nbfWedgedSubImage3D;

/** Interface for VTK-like input-output pipeline filters.
	Update state is kept internally so execution is only done when needed.
	User is responsible for changing the state when doing changes that affect the filter's output.
*/
template< class Pixel >
class nbfWedgedAverageImage3D : public nbfWedgedImage3D< Pixel >
{
public:

	nbfWedgedAverageImage3D(){ this->initialize(); }
	nbfWedgedAverageImage3D( vector< nbfWedgedSubImage3D< Pixel > > & v ){ this->initialize(); this->volumes = v; }
	~nbfWedgedAverageImage3D();

	// copy constructor
	nbfWedgedAverageImage3D( const nbfWedgedAverageImage3D< Pixel > & c ){ (*this) = c; }

	// WARNING - ASSUME CLASS AVERAGES NEVER HAVE A MISSING WEDGE
	void initialize(){ this->average.free(); this->sphericalWedgeImage.free(); this->wedgeAccum.free(); this->wedgeImage.free(); this->updateState(); this->smoothBeforeAveraging = 0; }

	nbfWedgedAverageImage3D< Pixel > & operator= ( const nbfWedgedAverageImage3D< Pixel > & );

	// rotates the current average without re-computing averages
	void rotate( vtkMatrix4x4 *, nbfProjectionRotationMetric3D< Pixel > * );

	/// Redefine from father
	void getImage( Array< Pixel, 3 > &, vtkTransform * = NULL, bool = false );
	virtual void getImage( vtkImageData *, vtkTransform * = NULL, bool = false );
	void getImageMPI( Array< double, 3 > &, Array< double, 3 > &, Array< double, 3 > & );

	void getImageVariance( Array< Pixel, 3 > & );
	virtual void getImageVariance( vtkImageData * );
	void getImageVarianceMPI( Array< double, 3 > & );

	void setAccumulatedWedgeImage( Array< Pixel, 3 > & A ){
		
		if ( A.size() > 0 ){
			this->wedgeAccum.resize( A.shape() ); 
			this->wedgeAccum = A; 
		}

		this->wedgeImage.resize( this->wedgeAccum.shape() );
		this->wedgeImage = where( this->wedgeAccum > max( this->wedgeAccum ) / 2.0, 1, 0 );

		//this->wedgeEffective = ( min( this->wedgeImage ) != 1 );
		this->wedgeEffective = sum( this->wedgeImage ) < .95 * this->wedgeImage.size();
		if ( this->wedgeEffective == true ){
			this->wedge.smoothWedgeImage( this->wedgeImage );
		} else {
			this->wedge.reset();
			this->wedgeAccum.free();
			this->wedgeImage.free();
		}

		this->wedgeUpToDate = true;
	}

	void setAverageImage( Array< double, 3 > & A ){
		this->average.resize( A.shape() ); 
		this->average = A;
		this->imageUpToDate = true;
		this->imageVarianceUpToDate = false;
	}

	void setVarianceImage( Array< double, 3 > & A ){
		this->variance.resize( A.shape() ); 
		this->variance = A;
		this->imageVarianceUpToDate = true;
	}

	/// Redefine from father
	void getWedgeImage( Array< Pixel, 3 > &, vtkTransform * = NULL );
	void updateAccumulatedWedgeImage();
	void getAccumulatedWedgeImage( Array< Pixel, 3 > & , vtkTransform * = NULL );

	void updateAccumulatedSphericalWedgeImageMPI( TinyVector< int, 2 > &, nbfProjectionRotationMetric3D< Pixel > * );
	void updateSphericalWedgeImage( TinyVector< int, 2 > &, nbfProjectionRotationMetric3D< Pixel > * );
	// void getSphericalWedgeImage( Array< Pixel, 2 > & A, vtkTransform * B ){ this->getSphericalWedgeImage( A, B, (nbfProjectionRotationMetric3D<Pixel>*)(NULL) ); }
	void getSphericalWedgeImage( Array< Pixel, 2 > &, vtkTransform * = NULL, nbfProjectionRotationMetric3D< Pixel > * = NULL );
	void setSphericalWedgeImage( Array< Pixel, 2 > & );
	
	void getWedgeImageHalf( Array< Pixel, 3 > &, vtkTransform * = NULL );

	void setAlignments( Array< Pixel, 3 > & );

	// reset volume and weights list
	void clear(){ this->getVolumes().clear(); this->weights = 0; }

	vector< nbfWedgedSubImage3D< Pixel > > & getVolumes(){ this->updateState(); return this->volumes; }
	vector< nbfWedgedSubImage3D< Pixel > > & getVolumesRO(){ return this->volumes; }
	
	Array< double, 3 > average;
	Array< double, 3 > variance;
	Array< Pixel, 3 > wedgeImage;
	Array< Pixel, 2 > sphericalWedgeImage;
	Array< Pixel, 2 > sphericalWedgeAccum;
	Array< Pixel, 3 > wedgeAccum;

	/// Redefine from nbfImageMetric
	TinyVector< Pixel, 3 > getDimensions();

	// MPI
	int getTypeId(){ return NBF_WEDGED_AVERAGE_IMAGE_3D; }
	void serialize( stringstream & );
	void unserialize( stringstream & );

	// weights
	Array< Pixel, 2 > weights;

	// multiple alignments
	Array< Pixel, 3 > multipleAlignments;

	// check is transform is within valid limits ( translation allowance, rotation allowance )
	bool isTransformValid( vtkTransform *, Pixel, Pixel );

	static nbfFourierFilter< Pixel, 3 > fourierFilter;

	void updateState(){ 
		this->imageUpToDate = false; 
		this->imageVarianceUpToDate = false; 
		this->wedgeUpToDate = false; 
		this->sphericalWedgeUpToDate = false;
		this->wedgeEffective = true;
	}

	bool isImageUpToDate(){ return this->imageUpToDate; }
	bool isImageVarianceUpToDate(){ return this->imageVarianceUpToDate; }
	
	bool isWedgeUpToDate(){ return this->wedgeUpToDate; }
	bool isSphericalWedgeUpToDate(){ return this->sphericalWedgeUpToDate; }

	void setSphericalWedgeUpToDate(){ this->sphericalWedgeUpToDate = true; }

	bool isWedgeEffective( TinyVector< int, 3 > & d );
	bool isWedgeEffective(){ TinyVector< int, 3 > d(this->volumes[0].getDimensions()); return this->isWedgeEffective(d); }

	void resetWedge(){ this->wedgeUpToDate = true; this->sphericalWedgeUpToDate = true; this->wedgeEffective = false; }

	Pixel smoothBeforeAveraging;

protected:

	bool imageUpToDate;
	bool imageVarianceUpToDate;
	bool wedgeUpToDate;
	bool sphericalWedgeUpToDate;
	bool wedgeEffective;

	vector< nbfWedgedSubImage3D< Pixel > > volumes;
	
	nbfWedgedSubImage3D< Pixel > volume;
};

template< class Pixel > nbfFourierFilter< Pixel, 3 > nbfWedgedAverageImage3D< Pixel > :: fourierFilter;

template< class Pixel >
nbfWedgedAverageImage3D< Pixel > :: ~nbfWedgedAverageImage3D()
{
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: getImage( Array< Pixel, 3 > & A, vtkTransform * t, bool normalize )
{
	cerr << "ERROR - NOT IMPLEMENTED: " << __FILE__ << ", " << __LINE__ << endl;
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: setAlignments( Array< Pixel, 3 > & alignments )
{
	if ( ( alignments.rows() != this->volumes.size() ) || ( alignments.cols() != 17 ) ){
		cout << "WARNING - Alignment dimensions do not match. Ignoring alignments." << endl;
		return;
	}

	// store weights
	this->weights.resize( alignments.rows(), alignments.depth() );
	this->weights = alignments( Range::all(), 0, Range::all() );

	this->multipleAlignments.resize( alignments.rows(), 16, alignments.depth() );
	this->multipleAlignments = alignments( Range::all(), Range(1,toEnd), Range::all() );

	for ( int i = 0; i < this->volumes.size(); i++ ){
		for ( int candidate = 0; candidate < this->multipleAlignments.depth(); candidate++ ){
			if ( abs( this->weights( i, candidate ) ) > 0 ){
				// set current alignment
				double matrix[16];
				for ( int k = 0; k < 16; k++ ){
					matrix[k] = this->multipleAlignments( i, k, candidate );
				}
				vtkMatrix4x4 * store = vtkMatrix4x4 :: New();
				store->DeepCopy( matrix );
				this->volumes[i].setTransform( store );
				store->Delete();
			}
		}
	}

	// force state to not-up-to-date so that average is recomputed next-time
	this->updateState();
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: getImage( vtkImageData * a, vtkTransform * t, bool normalize )
{
	if ( this->imageUpToDate == false ){

		Array< double, 3 > realFFT, imagFFT, wedgeAccum;
		this->getImageMPI( realFFT, imagFFT, wedgeAccum );

#ifdef NBF_AVERAGE_IN_RECIPROCAL_SPACE
		// define view of complex representation
		Array< complex< double >, 3 > FFTre( reinterpret_cast<complex<double>*>(nbfWedgedAverageImage3D< Pixel > :: fourierFilter.blitzFFT.data()), realFFT.shape(), neverDeleteData );
		real(FFTre) = where( wedgeAccum > 0, realFFT / wedgeAccum, 0 );
		imag(FFTre) = where( wedgeAccum > 0, imagFFT / wedgeAccum, 0 );
		//if ( normalize ){
			FFTre( FFTre.rows() / 2, FFTre.cols() / 2, FFTre.depth() - 1 ) = 0;
		//}
		fftw_execute( nbfWedgedAverageImage3D< Pixel > :: fourierFilter.ifftplanreal );
		this->average.resize( nbfWedgedAverageImage3D< Pixel > :: fourierFilter.blitzFFTreal.shape() );
		this->average = nbfWedgedAverageImage3D< Pixel > :: fourierFilter.blitzFFTreal / this->average.size() * nbfWedgedAverageImage3D< Pixel > :: fourierFilter.shift;
#else
		this->average.resize( realFFT.shape() );
		this->average = realFFT;
		// this->average /= sum( this->weights );
#endif
		//// normalize to zero mean and unit variance
		//if ( normalize ){
			Pixel varianceInside = sqrt( mean( this->average * this->average ) );
			this->average /= varianceInside;
		//}
		this->imageUpToDate = true;
		this->imageVarianceUpToDate = false;
	}

	if ( t != NULL ){
		nbfWedgedSubImage3D< Pixel > volumefixed;
		Array< Pixel, 3 > faverage( this->average.shape() );
		faverage = cast< Pixel >( this->average );
		volumefixed.setFixedImage( faverage );
		volumefixed.getImage( a, t );
	} else {
		nbfVTKInterface :: blitzToVtk( this->average, a );
	}
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: getImageMPI( Array< double, 3 > & realAverage, Array< double, 3 > & imagAverage, Array< double, 3 > & accumulatedFourier )
{
	if ( this->volumes.size() != this->weights.rows() ){
		cerr << "ERROR - Weights do not match current volumes. " << this->volumes.size() << "!=" << this->weights.rows() << ". In " __FILE__ << "," << __LINE__ << endl;
		return;
	}

	if ( this->imageUpToDate == false ){

		if ( ( this->volumes.size() > 0 ) && ( this->volumes.size() == this->weights.rows() ) ){

			TinyVector< int, 3 > size( volumes[0].getDimensions() );

#ifdef NBF_AVERAGE_IN_RECIPROCAL_SPACE
			TinyVector< int, 3 > shape( size[0], size[1], size[2] / 2 + 1 );

			// store reciprocal average
			realAverage.resize( shape );
			realAverage = 0;
			imagAverage.resize( shape );
			imagAverage = 0;
			accumulatedFourier.resize( shape );
			accumulatedFourier = 0;
			
			nbfWedgedAverageImage3D< Pixel > :: fourierFilter.initializeFFTWhalf( size );

			// define view of complex representation
			Array< complex< double >, 3 > FFTre( reinterpret_cast<complex<double>*>(nbfWedgedAverageImage3D< Pixel > :: fourierFilter.blitzFFT.data()), shape, neverDeleteData );

			// temporal array
			Array< double, 3 > B( shape );
#else
			realAverage.resize( size );
			realAverage = 0;
			imagAverage.free();
			accumulatedFourier.free();
#endif
			vtkImageData * image = vtkImageData::New();

			// traverse all volumes
			for ( int i = 0; i < this->volumes.size(); i++ ){

				for ( int candidate = 0; candidate < this->weights.cols(); candidate++ ){

					if ( abs( this->weights( i, candidate ) ) > 0 ){

						// get current volume image
						this->volume = this->volumes[i];

						// set current alignment
						double matrix[16];
						for ( int k = 0; k < 16; k++ ){
							matrix[k] = this->multipleAlignments( i, k, candidate );
						}
						vtkTransform * store = vtkTransform::New();
						store->SetMatrix( matrix );
						this->volume.setTransform( (vtkTransform*)NULL );
						this->volume.getImage( image, store, false );
						store->Delete();

						//if ( this->smoothBeforeAveraging > 0 ){
						//	vtkImageGaussianSmooth * filter = vtkImageGaussianSmooth::New();
						//	filter->SetDimensionality(3);
						//	filter->SetInput( image );
						//	filter->SetStandardDeviations( this->smoothBeforeAveraging, this->smoothBeforeAveraging, this->smoothBeforeAveraging );
						//	filter->SetRadiusFactors( this->smoothBeforeAveraging * .75, this->smoothBeforeAveraging * .75, this->smoothBeforeAveraging * .75 );
						//	image->DeepCopy( filter->GetOutput() );
						//	filter->Delete();
						//}

						Array< double, 3 > T;
						nbfVTKInterface::vtkToBlitzReference( image, T );
#ifdef NBF_AVERAGE_IN_RECIPROCAL_SPACE
						nbfWedgedAverageImage3D< Pixel > :: fourierFilter.blitzFFTreal = T * nbfWedgedAverageImage3D< Pixel > :: fourierFilter.shift;
						fftw_execute( nbfWedgedAverageImage3D< Pixel > :: fourierFilter.fftplanreal );

						B = sqrt( real( FFTre * conj( FFTre ) ) );
						accumulatedFourier += abs( this->weights(i,candidate) ) * B;

						//// frequency weighted average						
						//M = exp( - () / * );
						//accumulatedFourier += M * B;						

						realAverage += this->weights(i,candidate) * B * real( FFTre );
						imagAverage += this->weights(i,candidate) * B * imag( FFTre );
#else
						realAverage += this->weights(i,candidate) * T / sum( abs( this->weights ) );
#endif
					}
				}
			}
			image->Delete();
		}
		this->imageUpToDate = true;
	}
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: rotate( vtkMatrix4x4 * trans, nbfProjectionRotationMetric3D< Pixel > * metric )
{
	// apply to all volumes in current reference

	// update multiple alignments
	for ( int m = 0; m < this->weights.rows(); m++ ){
		for ( int n = 0; n < this->weights.cols(); n++ ){
			if ( abs( this->weights( m, n ) > 0 ) ){
				// retrieve current alingments
				double matrix[16];
				for ( int k = 0; k < 16; k++ ){
					matrix[k] = this->multipleAlignments( m, k, n );
				}
				vtkMatrix4x4 * t = vtkMatrix4x4 :: New();
				t->DeepCopy( matrix );
				vtkMatrix4x4 * s = vtkMatrix4x4 :: New();
				vtkMatrix4x4 :: Multiply4x4( t, trans, s );
				t->Delete();
				vtkMatrix4x4 :: DeepCopy( matrix, s );
				s->Delete();
				for ( int k = 0; k < 16; k++ ){
					this->multipleAlignments( m, k, n ) = matrix[k];
				}
			}
		}
	}

	vtkTransform * t1 = vtkTransform :: New();
	t1->SetMatrix( trans );

	// update average
	if ( this->imageUpToDate == true ){
		vtkImageData * data = vtkImageData :: New();
		this->getImage( data, t1 );
		Array< double, 3 > A;
		nbfVTKInterface :: vtkToBlitzReference( data, A );
		this->setAverageImage( A );
		data->Delete();
	}
	t1->Delete();

	vtkTransform * t2 = vtkTransform :: New();
	t2->SetMatrix( trans );
	// permanently rotate missing wedge
	if ( this->wedgeUpToDate == true ){
		Array< Pixel, 3 > W( this->wedgeImage.shape() );
		this->getWedgeImage( W, t2 );
		this->wedgeImage = W;
	}
	t2->Delete();

	vtkTransform * t3 = vtkTransform :: New();
	t3->SetMatrix( trans );
	// permanently rotate spherical missing wedge
	if ( this->sphericalWedgeUpToDate == true ){
		Array< Pixel, 2 > S( this->sphericalWedgeImage.shape() );
		this->getSphericalWedgeImage( S, t3, metric );
		this->sphericalWedgeImage = S;
		//nbfMatlabWriter w;
		//w.setFileName("p.matlab");
		//w.write(S);
	}
	t3->Delete();
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: getImageVariance( Array< Pixel, 3 > & A )
{
	cerr << "ERROR - NOT IMPLEMENTED: " << __FILE__ << ", " << __LINE__ << endl;
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: getImageVariance( vtkImageData * a )
{
	if ( this->imageVarianceUpToDate == false ){
		this->getImageVarianceMPI( this->variance );
	}
	nbfVTKInterface::blitzToVtk( this->variance, a ); 
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: getImageVarianceMPI( Array< double, 3 > & varianceImage )
{
	if ( this->imageUpToDate == false ){
		cerr << "ERROR - Global average needs to be computed before variance computation." << endl;
	}

	if ( this->volumes.size() != this->weights.rows() ){
		cerr << "ERROR - Weights do not match current volumes." << endl;
		return;
	}

	if ( this->imageVarianceUpToDate == false ){

		if ( ( this->volumes.size() > 0 ) && ( this->volumes.size() == this->weights.rows() ) ){

			TinyVector< int, 3 > size( volumes[0].getDimensions() );
			
			varianceImage.resize( size );
			varianceImage = 0;

			vtkImageData * image = vtkImageData::New();

			// traverse all volumes
			for ( int i = 0; i < this->volumes.size(); i++ ){

				for ( int candidate = 0; candidate < this->weights.cols(); candidate++ ){

					if ( abs( this->weights( i, candidate ) ) > 0 ){

						// get current volume image
						this->volume = this->volumes[i];

						// set current alignment
						double matrix[16];
						for ( int k = 0; k < 16; k++ ){
							matrix[k] = this->multipleAlignments( i, k, candidate );
						}
						vtkTransform * store = vtkTransform::New();
						store->SetMatrix( matrix );
						this->volume.setTransform( (vtkTransform*)NULL );
						this->volume.getImage( image, store, false );
						store->Delete();

						Array< double, 3 > T;
						nbfVTKInterface::vtkToBlitzReference( image, T );
						varianceImage += abs( this->weights(i,candidate) ) * pow2( T - this->average ) / sum( abs( this->weights ) );
					}
				}
			}
			image->Delete();
		}
		this->imageVarianceUpToDate = true;
	}
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: updateAccumulatedWedgeImage()
{
	if ( this->volumes.size() != this->weights.rows() ){
		cerr << "ERROR - Weights do not match current volumes. " << this->volumes.size() << "!=" << this->weights.rows() << ". In " __FILE__ << "," << __LINE__ << endl;
		return;
	}

	if ( this->wedgeUpToDate == false ){

			TinyVector< int, 3 > size( this->volumes[0].getDimensions() );

			Array< Pixel, 3 > accum( size ), current( size );
			accum = 0;

			// traverse all volumes
			for ( int i = 0; i < this->volumes.size(); i++ ){

				for ( int candidate = 0; candidate < this->weights.cols(); candidate++ ){

					if ( abs( this->weights( i, candidate ) ) > 0 ){

						// get current volume image
						this->volume = this->volumes[i];

						// set current alignment
						double matrix[16];
						for ( int k = 0; k < 16; k++ ){
							matrix[k] = this->multipleAlignments( i, k, candidate );
						}
						vtkTransform * store = vtkTransform::New();
						store->SetMatrix( matrix );
						this->volume.setTransform( (vtkTransform*)NULL );
						this->volume.getWedgeImage( current, store );
						store->Delete();

						accum += current;
					}
				}
			}

			this->setAccumulatedWedgeImage( accum );
	}
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: getAccumulatedWedgeImage( Array< Pixel, 3 > & A, vtkTransform * t )
{
	this->updateAccumulatedWedgeImage();
	A = this->wedgeAccum;
}

template< class Pixel >
bool nbfWedgedAverageImage3D< Pixel > :: isWedgeEffective( TinyVector< int, 3 > & shape )
{
	this->updateAccumulatedWedgeImage();
	return this->wedgeEffective;
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: getWedgeImage( Array< Pixel, 3 > & A, vtkTransform * t )
{
	TinyVector< int, 3 > shape = A.shape();
	if ( this->isWedgeEffective( shape ) ){
		A = this->wedgeImage;
	} else {
		A = 1;
		return;
	}

	// rotate wedge image (if neccessary)
	if ( t != NULL ){

		// build new transform with rotation component only
		vtkTransform * rotation = vtkTransform::New();
		rotation->SetInput( t );
		double position[3];
		t->GetPosition(position);
		position[0] = - position[0];
		position[1] = - position[1];
		position[2] = - position[2];
		rotation->PostMultiply();
		rotation->Translate(position);

		vtkImageData * data = vtkImageData :: New();
		nbfVTKInterface :: blitzToVtk( A, data );

		vtkImageChangeInformation * change = vtkImageChangeInformation :: New();
		change->SetInput( data );
		change->SetOriginTranslation( - A.rows() / 2.0 , - A.cols() / 2.0, - A.depth() / 2.0 );

		vtkImageReslice * reslice = vtkImageReslice :: New();
		reslice->SetInput( change->GetOutput() );
		reslice->SetBackgroundLevel( max(A) );
		reslice->SetInterpolationModeToCubic();
		reslice->SetResliceTransform( rotation );
		reslice->Update();

		nbfVTKInterface :: vtkToBlitz( reslice->GetOutput(), A );
		reslice->Delete();
		change->Delete();
		data->Delete();
		rotation->Delete();
	}
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: getWedgeImageHalf( Array< Pixel, 3 > & A, vtkTransform * t )
{
	this->getWedgeImage( A, t );

	Array< Pixel, 3 > B( A.rows(), A.cols(), A.depth() / 2 + 1 );
	B = A( Range :: all(), Range :: all(), Range( 0, A.depth() / 2 ) );
	A.resize( B.shape() );
	A = B;
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: updateAccumulatedSphericalWedgeImageMPI( TinyVector< int, 2 > & size, nbfProjectionRotationMetric3D< Pixel > * metric )
{
	if ( this->volumes.size() != this->weights.rows() ){
		cerr << "ERROR - Weights do not match current volumes. " << this->volumes.size() << "!=" << this->weights.rows() << ". In " __FILE__ << "," << __LINE__ << endl;
		return;
	}

	this->sphericalWedgeAccum.resize( size );
	Array< Pixel, 2 > A( size );

	// initialize to 0
	this->sphericalWedgeAccum = 0;

	// combine wedge from all volumes
	for ( int i = 0; i < this->volumes.size(); i++ ){
		for ( int candidate = 0; candidate < this->weights.cols(); candidate++ ){
			if ( abs( this->weights(i,candidate) ) > 0 ){

				double matrix[16];
				for ( int k = 0; k < 16; k++ ){
					matrix[k] = this->multipleAlignments( i, k, candidate );
				}

				vtkTransform * store = vtkTransform::New();
				store->SetMatrix( matrix );

				this->volume = this->volumes[i];
				this->volume.setTransform( (vtkTransform*)NULL );

				// get current volume image
				this->volume.getSphericalWedgeImage( A, store, metric );

				store->Delete();

				this->sphericalWedgeAccum += A;
			}
		}
	}
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: updateSphericalWedgeImage( TinyVector< int, 2 > & size, nbfProjectionRotationMetric3D< Pixel > * metric )
{
	if ( this->sphericalWedgeUpToDate == false ){

		this->updateAccumulatedSphericalWedgeImageMPI( size, metric );

		// keep a wedge if less than half the max contributors
		this->sphericalWedgeImage.resize( size );
		this->sphericalWedgeImage = where( this->sphericalWedgeAccum > max( this->sphericalWedgeAccum ) / 2.0, 1, 0 );
		// this->sphericalWedgeImage = where( this->sphericalWedgeAccum > 0, 1, 0 );
		
		this->sphericalWedgeUpToDate = true;
	}
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: getSphericalWedgeImage( Array< Pixel, 2 > & A, vtkTransform * t, nbfProjectionRotationMetric3D< Pixel > * metric )
{
	TinyVector< int, 2 > size = A.shape();
	if ( this->isWedgeEffective() ){
		this->updateSphericalWedgeImage( size, metric );

		//nbfMatlabWriter w;
		//w.setFileName("p.matlab");
		//w.write(this->sphericalWedgeImage);

		if ( t != NULL ){
			// extract angles from transformation
			metric->harmonics( this->sphericalWedgeImage, t, A );
		} else {
			//A = where( this->sphericalWedgeImage > .5, 1, 0 );
			A = this->sphericalWedgeImage;
		}
	} else {
		A = 1;
	}
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: setSphericalWedgeImage( Array< Pixel, 2 > & A )
{
	this->sphericalWedgeImage.resize( A.shape() );
	this->sphericalWedgeImage = A;
	this->sphericalWedgeUpToDate = true;
}

template< class Pixel >
nbfWedgedAverageImage3D< Pixel > & nbfWedgedAverageImage3D< Pixel > :: operator= ( const nbfWedgedAverageImage3D< Pixel > & param )
{
	this->wedge = param.wedge;

	this->volumes = param.volumes;
	this->weights.resize( param.weights.shape() );
	this->weights = param.weights;

	this->multipleAlignments.resize( param.multipleAlignments.shape() );
	this->multipleAlignments = param.multipleAlignments;

	this->volume = param.volume;

	this->imageUpToDate = param.imageUpToDate;
	this->imageVarianceUpToDate = param.imageVarianceUpToDate;
	this->wedgeUpToDate = param.wedgeUpToDate;
	this->sphericalWedgeUpToDate = param.sphericalWedgeUpToDate;
	this->wedgeEffective = param.wedgeEffective;

	this->average.resize( param.average.shape() );
	this->average = param.average;

	this->variance.resize( param.variance.shape() );
	this->variance = param.variance;

	this->wedgeImage.resize( param.wedgeImage.shape() );
	this->wedgeImage = param.wedgeImage;

	this->sphericalWedgeImage.resize( param.sphericalWedgeImage.shape() );
	this->sphericalWedgeImage = param.sphericalWedgeImage;

	this->sphericalWedgeAccum.resize( param.sphericalWedgeAccum.shape() );
	this->sphericalWedgeAccum = param.sphericalWedgeAccum;

	this->wedgeAccum.resize( param.wedgeAccum.shape() );
	this->wedgeAccum = param.wedgeAccum;

	this->smoothBeforeAveraging = param.smoothBeforeAveraging;

	return *this;
}

template< class Pixel >
TinyVector< Pixel, 3 > nbfWedgedAverageImage3D< Pixel > :: getDimensions(){
	TinyVector< Pixel, 3 > dim;
	if ( volumes.size() != 0 ){
		dim = volumes[0].getDimensions();
	} else {
		dim = this->average.shape();
	}
	return dim;
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: serialize( stringstream & output ){
	output << this->wedge.lwedge << endl;
	output << this->wedge.uwedge << endl;
	output << this->imageUpToDate << endl;
	output << this->imageVarianceUpToDate << endl;
	output << this->wedgeUpToDate << endl;
	output << this->sphericalWedgeUpToDate << endl;
	output << this->wedgeEffective << endl;
	output << this->average << endl;
	output << this->variance << endl;
	output << this->wedgeImage << endl;
	output << this->wedgeAccum << endl;
	output << this->sphericalWedgeImage << endl;
	output << this->sphericalWedgeAccum << endl;
	output << this->volumes.size() << endl;
	for ( int i = 0; i < this->volumes.size(); i++ ){
		this->volumes[i].serialize(output);
	}
	output << this->weights << endl;
	output << this->multipleAlignments << endl;
}

template< class Pixel >
void nbfWedgedAverageImage3D< Pixel > :: unserialize( stringstream & input ){
	Pixel lwedge, uwedge;
	input >> lwedge;
	input >> uwedge;
	this->wedge.set( lwedge, uwedge );
	input >> this->imageUpToDate;
	input >> this->imageVarianceUpToDate;
	input >> this->wedgeUpToDate;
	input >> this->sphericalWedgeUpToDate;
	input >> this->wedgeEffective;
	input >> this->average;
	input >> this->variance;
	input >> this->wedgeImage;
	input >> this->wedgeAccum;
	input >> this->sphericalWedgeImage;
	input >> this->sphericalWedgeAccum;
	unsigned int numVolumes;
	input >> numVolumes;
	this->volumes.clear();
	for ( int i = 0; i < numVolumes; i++ ){
		nbfWedgedSubImage3D< Pixel > tmp;
		tmp.unserialize( input );
		this->volumes.push_back( tmp );
	}
	input >> this->weights;
	input >> this->multipleAlignments;
}

template< class Pixel >
bool nbfWedgedAverageImage3D< Pixel > :: isTransformValid( vtkTransform * t, Pixel tAll, Pixel rAll )
{
	if ( this->volumes.size() > 0 ){
		return this->volumes[0].isTransformValid( t, tAll, rAll );
	} else {
		return false;
	}
}
