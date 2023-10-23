/*
 * Phantom.cpp
 *
 *  Created on: Jan 31, 2010
 *      Author: fefo
 */

#include "mpi.h"


#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkImageReslice.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkImageReader.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageReslice.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkImageNonMaximumSuppression.h>
#include <vtkImageMathematics.h>
#include <vtkImageContinuousDilate3D.h>
#include <vtkImageMathematics.h>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageCast.h>
#include <vtkImageShrink3D.h>
#include <vtkImageEllipsoidSource.h>
#include <vtkSphereSource.h>
#include <io/nbfVTKInterface.h>
#include <io/nbfMrcReader.h>
#include <io/nbfMrcWriter.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMrcWriter.h>
#include <bs/nbfBordStrategyMirror.h>

#include <em/nbfImageMetric.h>
#include <em/nbfFourierImageMetric.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfWedgedAverageImage3D.h>

#include <em/nbfUnsupervisedLoopClustering.h>

#include <nbfCylindricalDomain3.h>

#include <mxml.h>

#define PIXEL float

int main( int argc, char *argv[] ) {

	nbfImageFilter< PIXEL, 3 > imfilter;
	nbfFourierFilter< PIXEL, 3 > fffilter;

	nbfFourierImageMetric< PIXEL, 3 > fmetric( &imfilter, &fffilter );
	vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
	nbfWedgedSubImage3D< PIXEL > :: read( argv[1], volumeList );

	Array< complex< PIXEL >, 1 > C;
	Array< PIXEL, 2 > M;
	Array< PIXEL, 3 > Ri;
	for ( int i = 0; i < volumeList.size(); i++ ){
	
		cout << "Processing " << volumeList[i].getFileName() << endl;
		
		fmetric.setInput1( &volumeList[i] );
		fmetric.get2DimensionRepresentationHalf( C, M, i );
		
		nbfPolarDomain< PIXEL, 2 > polar;
		TinyVector< PIXEL , 2 > center( (M.rows()+1)/2, M.cols()-2 );
		polar.setCenter( center );
		// polar.setMinRho( M.cols() * this->metric->fourierFilter->getBandLowCut() );
		// polar.setMaxRho( M.cols() * this->metric->fourierFilter->getBandHighCut()  );
		polar.setMaxRho( M.cols() );
		polar.setResRho( 1 * M.cols() );
		polar.setResTheta( 360 );
		Array< PIXEL, 2 > P;
		Array< bool, 2 > B;
		polar.cartesian2polar( M, P, B );
		Array< PIXEL, 1 > Ps( P.rows() );
		secondIndex j;
		Ps = sum( P( Range :: all(), Range( Ps.cols()/2, toEnd ) ), j );	
		nbfMatlabWriter m;

		stringstream outputFile;
		int size = volumeList[i].getFileName().size();
		outputFile << volumeList[i].getFileName().erase(size-4,5) << "_psd.matlab";
		m.setFileName( outputFile.str().c_str() );
		m.write(Ps);
	}
	return 0;
}
