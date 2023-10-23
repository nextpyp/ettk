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

#include <io/nbfVTKInterface.h>
#include <em/nbfCorrelationImageMetric3D.h>
#include <nbfTimer.h>

using namespace blitz;

/** 3D volume alignment based on cross-correlation. Given two input volumes, a brute force search 
	is performed with rotated versions of the second image.
	Accounting for the wedge is supported by specifying a nbfWedge object with setting for both volumes.
	Search domain bounds and sampling rate are specified with nbfAngularSearchCorrelationMetricFFTW::setAngleSearch#.
	Search can be done on downsized volumes by specifying a scale factor.
	@todo Use more sophisticated search strategies to speed up brute force approach (manage multiple scales, 
	different sampling rates, etc).
	@see nbfCorrelationMetricFFTW
*/
template< class Pixel >
class nbfImageMetricEM3D
{
public:

	/// Constructor.
	nbfImageMetricEM3D( nbfImageMetric< Pixle, 3 > * );

	~nbfImageMetricEM3D();

	/// Set input EM data.
	void setInput1( nbfSubVolume< Pixel > * );
	void setInput2( nbfSubVolume< Pixel > * );

	/// Wrapper for nbfImageMetric
	void execute();

	/// Distance value between inputs (function of the correlation value).
	Pixel getDistance(){ return this->metric->getDistance(); }

	/// Get value of correlation (assuming execute() was already called).
	Pixel getCorrelationPeak(){ return this->metric->getCorrelationPeak(); }

protected:

	/// Search bounds and sampling rates for each axis.
	Pixel lboundX, uboundX;
	Pixel lboundY, uboundY;
	Pixel lboundZ, uboundZ;
	
	/// Store current best alignment candidate parameters.
	Pixel bestAngleX, bestAngleY, bestAngleZ;
	Pixel bestTranslationX, bestTranslationY, bestTranslationZ;
	Pixel bestCorrelationPeak;

	nbfSubVolume< Pixel > * input1;
	nbfSubVolume< Pixel > * input2;

	/// Store search strategy.
	Array< Pixel, 2 > strategy;

	/// Store metric to compare images
	nbfCorrelationImageMetric3D< Pixel > * metric;
};

template< class Pixel >
nbfAngularSearchMetricEM3D< Pixel > :: nbfAngularSearchMetricEM3D( nbfCorrelationimageMetric3D< Pixel > * m )
: lboundX(0), uboundX(0),
  lboundY(0), uboundY(0),
  lboundZ(0), uboundZ(0),
  nbfCorrelationImageMetric3D< Pixel >( m->getWedge() )
{
	this->metric = m;
}

template< class Pixel >
nbfAngularSearchMetricEM3D< Pixel > :: ~nbfAngularSearchMetricEM3D()
{
}

template< class Pixel >
void nbfAngularSearchMetricEM3D< Pixel > :: execute()
{
	vtkImageData * vinput1 = vtkImageData::New();
	this->input1->getVolume( vinput1 );

	vtkImageData * vinput2 = vtkImageData::New();
	this->input2->getVolume( vinput2 );

	this->metric->setInput1( vinput1 );
	this->metric->setInput2( vinput2 );

	this->metric->execute();

	TinyVector< int, 3 > translation;

	// initialize internal parameters
	this->bestAngleX = numeric_limits< Pixel > :: max();
	this->bestAngleY = numeric_limits< Pixel > :: max();
	this->bestAngleZ = numeric_limits< Pixel > :: max();
	this->bestCorrelationPeak = - numeric_limits< Pixel > :: max();

	this->bestAngleX = lboundX;
	this->bestAngleY = lboundY;
	this->bestAngleZ = lboundZ;
	
	// probe transform
	vtkTransform * probeRotation = vtkTransform :: New();

	// wedge transform
	vtkTransform * wedgeRotation = vtkTransform :: New();

	Pixel lboundXcurrent = this->lboundX;
	Pixel uboundXcurrent = this->uboundX;
	Pixel lboundYcurrent = this->lboundY;
	Pixel uboundYcurrent = this->uboundY;
	Pixel lboundZcurrent = this->lboundZ;
	Pixel uboundZcurrent = this->uboundZ;

	nbfWedgeFilter3D< Pixel > wedge;

	for ( int pass = 0; pass < this->strategy.rows(); pass++ ){

		Pixel scale = this->strategy(pass,firstDim);
		Pixel delta = this->strategy(pass,secondDim);

		double downSampling = 1.0 / pow( 2.0, scale - 1.0 );

		resliceInput2->SetAxisMagnificationFactor(0,downSampling);
		resliceInput2->SetAxisMagnificationFactor(1,downSampling);
		resliceInput2->SetAxisMagnificationFactor(2,downSampling);

		resliceInput1->SetAxisMagnificationFactor(0,downSampling);
		resliceInput1->SetAxisMagnificationFactor(1,downSampling);
		resliceInput1->SetAxisMagnificationFactor(2,downSampling);
		resliceInput1->Update();
	
		wedge.setDimensions( resliceInput1->GetOutput() );
		wedge.wedge1On( this->input1 );

		this->metric->setInput1( resliceInput1->GetOutput() );

		// brute force search
		for ( Pixel angleX = lboundXcurrent; angleX <= uboundXcurrent; angleX += delta ){
			for ( Pixel angleY = lboundYcurrent; angleY <= uboundYcurrent; angleY += delta ){
				for ( Pixel angleZ = lboundZcurrent; angleZ <= uboundZcurrent; angleZ += delta ){

					probeRotation->Identity();	
					probeRotation->RotateX(-angleX);
					probeRotation->RotateY(-angleY);
					probeRotation->RotateZ(-angleZ);

					// apply transform to second image
					resliceInput2->SetResliceTransform( probeRotation );
					resliceInput2->Update();

					this->metric->setInput2( resliceInput2->GetOutput() );

					wedgeRotation->Identity();	
					wedgeRotation->RotateX(angleX);
					wedgeRotation->RotateY(angleY);
					wedgeRotation->RotateZ(angleZ);

					// apply transform to second wedge only
					// this->wedge->wedge2Transform( wedgeRotation );

					// compute CCC
					this->metric->execute();

					// see if correlation peak is better
					if ( this->metric->getCorrelationPeak() > this->bestCorrelationPeak ){

						// store new attributes
						this->bestCorrelationPeak = this->metric->getCorrelationPeak();
						this->bestAngleX = angleX;
						this->bestAngleY = angleY;
						this->bestAngleZ = angleZ;

						double trans[3];
						this->metric->getTransform()->GetPosition(trans);
						this->bestTranslationX = trans[0] /= downSampling;
						this->bestTranslationY = trans[1] /= downSampling;
						this->bestTranslationZ = trans[2] /= downSampling;
					}

					//// debug info
					//double trans[3];
					//correlationMetric.getTransform()->GetPosition(trans);
					//cout << "theta = [" << angleX << "," << angleY << "," << angleZ <<"] , score = " << correlationMetric.getCorrelationPeak() << ", t = [ " << trans[0] << ", " << trans[1] << ", " << trans[2] << "]" << endl;
				}
			}
		}
		
		// if not done, set new search bounds
		if ( pass < this->strategy.ubound(firstDim) ){
			Pixel nextDelta = this->strategy( pass + 1, secondDim );
			lboundXcurrent = max( this->bestAngleX - nextDelta, this->lboundX );
			uboundXcurrent = min( this->bestAngleX + nextDelta, this->uboundX );
			lboundYcurrent = max( this->bestAngleY - nextDelta, this->lboundY );
			uboundYcurrent = min( this->bestAngleY + nextDelta, this->uboundY );
			lboundZcurrent = max( this->bestAngleZ - nextDelta, this->lboundZ );
			uboundZcurrent = min( this->bestAngleZ + nextDelta, this->uboundZ );
		}
	}

	//t.stop();
	//cout << "total = " << t.elapsedSeconds() << endl;

	probeRotation->Delete();

	cout << "Alignment = [" << bestAngleX << "," << bestAngleY << "," << bestAngleZ <<"] , t = [" << bestTranslationX << ", " << bestTranslationY << ", " << bestTranslationZ << " ], " << this->bestCorrelationPeak << endl;

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
	vtkImageFFT * fft1 = vtkImageFFT::New();
	fft1->SetDimensionality(3);	
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

	// set state to updated
	this->updated = true;
}