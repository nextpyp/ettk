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
#include <vtkCylindricalTransform.h>
#include <vtkImageMedian3D.h>
#include <vtkImageButterworthHighPass.h>
#include <vtkImageNoiseSource.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>

#include <nbfTimer.h>
#include <em/nbfDensityLocalMaximaClustering.h>
#include <em/nbfAngularSearchCorrelationMetric.h>
#include <em/nbfIterativeAveraging.h>

#include <vtkPNGReader.h>

#include <random/normal.h>

#define PIXEL float

void main( int argc, char ** argv  )
{
	if ( argc != 13 ){
		cout << "Usage: input angleX angleY angleZ wedgeLow wedgeHigh wedgeAngleX wedgeAngleY wedgeAngleZ noiseMin noiseMax output" << endl;
		exit(0);
	}
	
	char * fileIn = argv[1];

	double angleX = atof(argv[2]);
	double angleY = atof(argv[3]);
	double angleZ = atof(argv[4]);

	double wedgeLow = atof(argv[5]);
	double wedgeHigh = atof(argv[6]);
	double wedgeAngleX = atof(argv[7]);
	double wedgeAngleY = atof(argv[8]);
	double wedgeAngleZ = atof(argv[9]);

	double noiseMin = atof(argv[10]);
	double noiseMax = atof(argv[11]);

	char * fileOut = argv[12];

	// read input data
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName( fileIn );
	reader->Update();

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToDouble();
	cast->SetInput( reader->GetOutput() );

	vtkImageChangeInformation * center = vtkImageChangeInformation::New();
	center->CenterImageOn();
	center->SetInput( cast->GetOutput() );
	center->Update();

	vtkImageReslice * reslice = vtkImageReslice::New();
	reslice->WrapOff();
	reslice->MirrorOff();
	reslice->SetBackgroundLevel(0);
	reslice->SetInput( center->GetOutput() );
	reslice->SetInterpolationModeToCubic();

	vtkTransform * transform = vtkTransform::New();
	transform->RotateX( angleX );
	transform->RotateY( angleY );
	transform->RotateZ( angleZ );
	transform->Translate(0,0,0);
	reslice->SetResliceTransform( transform );
	reslice->Update();

	vtkMath::RandomSeed( (unsigned int)time(0) );

	vtkImageNoiseSource * noise = vtkImageNoiseSource::New();
	noise->SetMinimum( noiseMin );
	noise->SetMaximum( noiseMax );
	int extent[6];
	reader->GetOutput()->GetExtent(extent);
	noise->SetWholeExtent(extent[0],extent[1],extent[2],extent[3],extent[4],extent[5]);
	noise->Update();

	vtkImageMathematics * sum = vtkImageMathematics::New();
	sum->SetOperationToAdd();
	sum->SetInput1( noise->GetOutput() );
	sum->SetInput2( reslice->GetOutput() );
	sum->Update();

	// apply wedge filter
	nbfWedge< double > wedge;
	wedge.setGeometry( reader->GetOutput() );
	vtkTransform * t = vtkTransform::New();
	t->RotateX( wedgeAngleX );
	t->RotateY( wedgeAngleY );
	t->RotateZ( wedgeAngleZ );
	wedge.wedge1On( wedgeLow, wedgeHigh, t );
	vtkImageData * filtered = vtkImageData::New();
	wedge.applyReal( sum->GetOutput(), filtered );

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetInput( filtered );
	writer->SetFileName( fileOut );
	writer->Write();
	writer->Delete();

	center->Delete();
	reslice->Delete();
	sum->Delete();
	noise->Delete();
	t->Delete();
	transform->Delete();
	reader->Delete();
}