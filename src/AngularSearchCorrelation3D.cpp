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

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <nbfCylindricalDomain3.h>
#include <bs/nbfBordStrategyMirror.h>

#include <nbfTimer.h>
#include <em/nbfWedge.h>
#include <em/nbfAngularSearchCorrelationMetric.h>

void main( int argv, char ** argc )
{
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName(argc[1]);
	reader->Update();

	vtkImageData * image = vtkImageData :: New();
	image->DeepCopy( reader->GetOutput() );

	// apply wedge filter
	nbfWedge< float > wedge;
	wedge.setWedge(60);
	wedge.setGeometry( image );

	vtkTransform * t1 = vtkTransform::New();
	vtkTransform * t2 = vtkTransform::New();
	wedge.composeWedges(t1,t2);

	vtkImageData * output1 = vtkImageData::New();
	wedge.applyWedgeReal( image, output1 );

	int dims[3];
	image->GetDimensions(dims);
	cout << dims[0] << ", " << dims[1] << ", " << dims[2] << endl;

	// rotated volume
	vtkImageChangeInformation * change = vtkImageChangeInformation::New();
	change->CenterImageOn();
	change->SetInput( image );
	change->Update();

	vtkImageReslice * reslice1 = vtkImageReslice::New();
	reslice1->SetInput( change->GetOutput() );
	reslice1->WrapOff();
	reslice1->SetBackgroundLevel(0);

	vtkTransform * transform1 = vtkTransform::New();
	transform1->RotateX(0);
	transform1->RotateY(0);
	transform1->RotateZ(9);
	transform1->Translate(0,0,0);
	reslice1->SetResliceTransform( transform1 );
	reslice1->Update();

	vtkImageData * output2 = vtkImageData::New();
	wedge.applyWedgeReal( reslice1->GetOutput(), output2 );

	nbfAngularSearchCorrelationMetric< float > ccc;
	ccc.setWedge( &wedge );
	ccc.setAngleSearchX(0,0,1);
	ccc.setAngleSearchY(0,0,1);
	ccc.setAngleSearchZ(0,11,1);
	ccc.setInput1( output1 );
	ccc.setInput2( output2 );

	nbfTimer t;
	t.start();
	ccc.execute();
	t.stop();
	cout << "running time = " << t.elapsedSeconds() << endl;

	 cout << "distance = " << ccc.getDistance() << endl;

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetInput( output1 );
	writer->SetFileName( argc[2] );
	writer->Write();
	writer->Delete();

}