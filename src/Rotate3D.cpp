#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkMath.h>
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
#include <vtkImageCast.h>
#include <vtkImageConstantPad.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <nbfVeselnessFilter.h>
#include <bs/nbfBordStrategyMirror.h>

//#include <em/nbfAngularSearchCorrelationMetric.h>
//#include <em/nbfAngularSearchCorrelationMetricFFTW.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfAngularSearchMetric3D.h>

#define PIXEL double

void main( int argc, char ** argv )
{
	if ( argc != 7 ){
		cout << "Usage: " << argv[0] << " input rotZ rotY rotZ flip output" << endl;
		exit(0);
	}

	char * inFile = argv[1];
	double angleZ1 = atof(argv[2]);
	double angleY = atof(argv[3]);
	double angleZ2 = atof(argv[4]);
	double flip = atof(argv[5]);
	char * outFile = argv[6];

	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();

	// read 3D image data
	reader->SetFileName( inFile );
	reader->Update();

	///////////////////

	vtkImageChangeInformation * center = vtkImageChangeInformation::New();
	center->CenterImageOn();
	center->SetInput( reader->GetOutput() );
	center->Update();

	vtkImageReslice * reslice = vtkImageReslice::New();
	//reslice->WrapOff();
	//reslice->MirrorOff();

	Array< PIXEL, 3 > A;
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), A );
	reslice->SetBackgroundLevel( mean(A) );
	reslice->SetInput( center->GetOutput() );
	reslice->SetInterpolationModeToCubic();

	vtkTransform * transform = vtkTransform::New();
	transform->RotateZ( -angleZ1 );		// -gamma
	transform->RotateY( -angleY  );		// -beta
	transform->RotateZ( -angleZ2 );		// -alpha
	if ( flip != 0 ){
		cout << "flipping...";
		transform->RotateZ(180);
		transform->RotateY(180);
		cout << "done." << endl;
	}

	reslice->SetResliceTransform( transform );
	reslice->Update();

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName( outFile );
	writer->SetInput( reslice->GetOutput() );
	writer->Write();
}
