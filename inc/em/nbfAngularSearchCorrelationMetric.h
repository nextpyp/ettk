#pragma once

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkBMPWriter.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageShrink3D.h>
#include <vtkImageEllipsoidSource.h>
#include <vtkImageGaussianSmooth.h>

//#include <itkVTKImageToImageFilter.txx>
//#include <itkFFTWRealToComplexConjugateImageFilter.txx>

#include <io/nbfVTKInterface.h>
#include <em/nbfCorrelationImageMetric.h>
#include <nbfTimer.h>

using namespace blitz;

/** 3D volume alignment based on cross-correlation. Given two input volumes, a brute force search 
	is performed with rotated versions of the second image.
	Accounting for the wedge is supported by specifying a nbfWedge object with setting for both volumes.
	Search domain bounds and sampling rate are specified with nbfAngularSearchCorrelationMetric::setAngleSearch#.
	Search can be done on downsized volumes by specifying a scale factor.
	@todo Use more sophisticated search strategies to speed up brute force approach (manage multiple scales, 
	different sampling rates, etc).
	@see nbfCorrelationMetric
*/
template< class Pixel >
class nbfAngularSearchCorrelationMetric : public nbfCorrelationImageMetric< Pixel, 3 >
{
public:

	/// Constructor with wedge argument. A wedge object is required for all operations in reciprocal space.
	nbfAngularSearchCorrelationMetric( nbfWedge< Pixel > * );

	~nbfAngularSearchCorrelationMetric();

	/// Redefine from nbfCorrelationMetric
	void execute();

	/// Set angular search bounds and sampling rates.
	void setAngleSearchX( Pixel lbound, Pixel ubound, Pixel rate = 1.0 ){ this->lboundX = lbound; this->uboundX = ubound; this->rateX = rate; this->updated = false; }
	void setAngleSearchY( Pixel lbound, Pixel ubound, Pixel rate = 1.0 ){ this->lboundY = lbound; this->uboundY = ubound; this->rateY = rate; this->updated = false; }
	void setAngleSearchZ( Pixel lbound, Pixel ubound, Pixel rate = 1.0 ){ this->lboundZ = lbound; this->uboundZ = ubound; this->rateZ = rate; this->updated = false; }

	/// Set angular search as range around given position (typically done after first brute force pass).
	void setAngleSearchRefineX( Pixel delta, Pixel rate ){ this->lboundX = this->bestAngleX - delta; this->uboundX =  this->bestAngleX + delta; this->rateX = rate; this->updated = false; }
	void setAngleSearchRefineY( Pixel delta, Pixel rate ){ this->lboundY = this->bestAngleY - delta; this->uboundY =  this->bestAngleY + delta; this->rateY = rate; this->updated = false; }
	void setAngleSearchRefineZ( Pixel delta, Pixel rate ){ this->lboundZ = this->bestAngleZ - delta; this->uboundZ =  this->bestAngleZ + delta; this->rateZ = rate; this->updated = false; }

	/// Set scale for carrying out computations ( 2 = half, 4 = one forth, etc. ).
	void setDownSampling( Pixel );
	void setDownSampling( Pixel, Pixel, Pixel );

protected:

	/// Search bounds and sampling rates for each axis.
	Pixel lboundX, uboundX, rateX;
	Pixel lboundY, uboundY, rateY;
	Pixel lboundZ, uboundZ, rateZ;
	
	/// Downsampling factors for each axis.
	Pixel downSamplingX, downSamplingY, downSamplingZ;

	/// Store current best alignment candidate parameters.
	Pixel bestAngleX, bestAngleY, bestAngleZ;
	Pixel bestTranslationX, bestTranslationY, bestTranslationZ;
	Pixel bestCorrelationPeak;
};

template< class Pixel >
nbfAngularSearchCorrelationMetric< Pixel > :: nbfAngularSearchCorrelationMetric( nbfWedge< Pixel > * w )
: lboundX(0), uboundX(0), rateX(1),
  lboundY(0), uboundY(0), rateY(1),
  lboundZ(0), uboundZ(0), rateZ(1),
  downSamplingX(1), downSamplingY(1), downSamplingZ(1),
  nbfCorrelationMetric< Pixel, 3 >(w)
{
	this->bestAngleX = numeric_limits< Pixel > :: max();
	this->bestAngleY = numeric_limits< Pixel > :: max();
	this->bestAngleZ = numeric_limits< Pixel > :: max();
	this->bestCorrelationPeak = - numeric_limits< Pixel > :: max();
}

template< class Pixel >
nbfAngularSearchCorrelationMetric< Pixel > :: ~nbfAngularSearchCorrelationMetric()
{
}

template< class Pixel >
void nbfAngularSearchCorrelationMetric< Pixel > :: setDownSampling( Pixel x, Pixel y, Pixel z )
{
	this->downSamplingX = 1.0 / x;
	this->downSamplingY = 1.0 / y;
	this->downSamplingZ = 1.0 / z;
	this->updated = false;
}

template< class Pixel >
void nbfAngularSearchCorrelationMetric< Pixel > :: setDownSampling( Pixel p )
{
	this->setDownSampling( p, p, p );
}

template< class Pixel >
void nbfAngularSearchCorrelationMetric< Pixel > :: execute()
{
	// store copy of first fourier image
	vtkImageData * imf1c = vtkImageData::New();

	// pre-center input2 for rotation
	vtkImageChangeInformation * centerInput2 = vtkImageChangeInformation::New();
	centerInput2->SetInput( this->input2 );
	centerInput2->CenterImageOn();
	centerInput2->Update();

	// input2 reslicer (set filling strategy to mirror, i.e. fill rotation gaps assuming a mirror image)
	vtkImageResample * resliceInput2 = vtkImageResample::New();
	resliceInput2->MirrorOn();

	// continue setting slicer parameters
	resliceInput2->SetInput( centerInput2->GetOutput() );
	resliceInput2->SetInterpolationModeToCubic();
	resliceInput2->SetAxisMagnificationFactor(0,this->downSamplingX);
	resliceInput2->SetAxisMagnificationFactor(1,this->downSamplingY);
	resliceInput2->SetAxisMagnificationFactor(2,this->downSamplingZ);

	// input1 reslicer
	vtkImageResample * resliceInput1 = vtkImageResample::New();
	resliceInput1->SetInput( this->input1 );
	resliceInput1->SetInterpolationModeToCubic();
	resliceInput1->SetAxisMagnificationFactor(0,this->downSamplingX);
	resliceInput1->SetAxisMagnificationFactor(1,this->downSamplingY);
	resliceInput1->SetAxisMagnificationFactor(2,this->downSamplingZ);
	resliceInput1->Update();

	// windowing of data in real space (in-place operation)
	this->window( resliceInput1->GetOutput() );

	this->wedge->setGeometry( resliceInput1->GetOutput() );

	//typedef itk::VTKImageToImageFilter< itk::Image< Pixel, 3 > > vtkToItkType;
	//vtkToItkType * vtkToItk = vtkToItkType::New();
	//vtkToItk->SetInput( resliceInput1->GetOutput() );
	//vtkToItk->Update();

	//typedef itk::FFTWRealToComplexConjugateImageFilter< Pixel, 3 >  FFTFilterType;
	//FFTFilterType * fft1 = FFTFilterType::New();
	//fft1->SetInput( vtkToItk->GetOutput() );

	// FFT of first image
	vtkImageFFT * fft1 = vtkImageFFT::New();
	fft1->SetDimensionality(3);
	fft1->SetInput( resliceInput1->GetOutput() );
	fft1->Update();

	//vtkToItk->Delete();

	// FFT of second image
	vtkImageFFT * fft2 = vtkImageFFT::New();
	fft2->SetDimensionality(3);
	fft2->SetInput( resliceInput2->GetOutput() );
	
	TinyVector< int, 3 > translation;

	if ( this->bestAngleX == numeric_limits< Pixel > :: max() ){
		this->bestAngleX = lboundX;
	}
	if ( this->bestAngleY == numeric_limits< Pixel > :: max() ){
		this->bestAngleY = lboundY;
	}
	if ( this->bestAngleZ == numeric_limits< Pixel > :: max() ){
		this->bestAngleZ = lboundZ;
	}

	// probe transform
	vtkTransform * probeRotation = vtkTransform :: New();

	// wedge transform
	vtkTransform * wedgeRotation = vtkTransform :: New();

	nbfTimer t;

	t.start();
	// brute force search
	for ( Pixel angleX = this->lboundX; angleX <= this->uboundX; angleX += this->rateX ){
		for ( Pixel angleY = this->lboundY; angleY <= this->uboundY; angleY += this->rateY ){
			for ( Pixel angleZ = this->lboundZ; angleZ <= this->uboundZ; angleZ += this->rateZ ){

				//// avoid maxima at origin location
				//if ( ( angleX == 0 ) && ( angleY == 0 ) && ( angleZ == 0 ) ){
				//	continue;
				//}

				probeRotation->Identity();	
				probeRotation->RotateX(-angleX);
				probeRotation->RotateY(-angleY);
				probeRotation->RotateZ(-angleZ);

				// apply transform to second image
				resliceInput2->SetResliceTransform( probeRotation );
				resliceInput2->Update();

				// windowing of data in real space (in-place operation)
				this->window( resliceInput2->GetOutput() );

				// FFT of second image
				fft2->Modified();
				fft2->Update();

				// operate on copy of first image's FFT to avoid recomputation of FFT.
				imf1c->DeepCopy( fft1->GetOutput() );
				
				wedgeRotation->Identity();	
				wedgeRotation->RotateX(angleX);
				wedgeRotation->RotateY(angleY);
				wedgeRotation->RotateZ(angleZ);

				// apply transform to second wedge only
				this->wedge->wedge2Transform( wedgeRotation );

				// compute CCC
				nbfCorrelationMetric< Pixel, 3 > :: execute( imf1c, fft2->GetOutput() );

				//cout << this->getCorrelationPeak() << endl;

				// see if correlation peak is better
				if ( this->getCorrelationPeak() > this->bestCorrelationPeak ){

					// store new attributes
					this->bestCorrelationPeak = this->getCorrelationPeak();
					this->bestAngleX = angleX;
					this->bestAngleY = angleY;
					this->bestAngleZ = angleZ;

					double trans[3];
					this->transform->GetPosition(trans);
					this->bestTranslationX = trans[0] /= this->downSamplingX;
					this->bestTranslationY = trans[1] /= this->downSamplingY;
					this->bestTranslationZ = trans[2] /= this->downSamplingZ;
				}

				//// debug info
				//double trans[3];
				//this->getTransform()->GetPosition(trans);
				//cout << "theta = [" << angleX << "," << angleY << "," << angleZ <<"] , score = " << this->getCorrelationPeak() << ", t = [ " << trans[0] << ", " << trans[1] << ", " << trans[2] << endl;
			}
		}
	}

	t.stop();
	cout << "total = " << t.elapsedSeconds() << endl;

	probeRotation->Delete();

	cout << "Alignment = [" << bestAngleX << "," << bestAngleY << "," << bestAngleZ <<"] , t = [" << bestTranslationX << ", " << bestTranslationY << ", " << bestTranslationZ << " ] - ccc = " << this->bestCorrelationPeak << endl;

	// re-assign attributes to reflect best match parameters
	this->transform->Identity();
	this->transform->RotateX( -bestAngleX );
	this->transform->RotateY( -bestAngleY );
	this->transform->RotateZ( -bestAngleZ );
	this->transform->Translate( -this->bestTranslationX, -this->bestTranslationY, -this->bestTranslationZ );
	
	this->correlationPeak = this->bestCorrelationPeak;

	// change strategy to mean constant value when data is missing after performing rotation+translation.
	Array< Pixel, 3 > M;
	nbfVTKInterface::vtkToBlitzReference( this->input2, M );
	resliceInput2->SetBackgroundLevel( mean(M) );

	// compute and store aligned image (using final transformation)
	resliceInput2->SetResliceTransform( this->transform );
	resliceInput2->Modified();
	resliceInput2->Update();
	this->aligned->DeepCopy( resliceInput2->GetOutput() );

	// compute FFT of aligned image
	fft1->SetInput( this->aligned );
	fft1->Modified();
	fft1->Update();
	
	// set wedge region to exactly zero
	wedgeRotation->Identity();
	wedgeRotation->RotateX(bestAngleX);
	wedgeRotation->RotateY(bestAngleY);
	wedgeRotation->RotateZ(bestAngleZ);
	this->wedge->wedge2Transform( wedgeRotation );
	
	Array< Pixel, 3 > real;
	Array< Pixel, 3 > imag;
	nbfVTKInterface::vtkToBlitz( fft1->GetOutput(), real, imag );
	this->wedge->applyComplexWedge2( real, imag );

	// store aligned image in reciprocal space
	this->alignedFourier->DeepCopy( fft1->GetOutput() );

	// garbage collect
	wedgeRotation->Delete();

	resliceInput1->Delete();
	resliceInput2->Delete();
	centerInput2->Delete();

	fft1->Delete();
	fft2->Delete();

	imf1c->Delete();

	// set state to updated
	this->updated = true;
}