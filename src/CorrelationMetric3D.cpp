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

#include <em/nbfAngularSearchCorrelationMetric.h>
#include <em/nbfWedge.h>

void main( int argv, char ** argc )
{
	//Array< double, 2 > M;
	//nbfMatlabReader mr;
	//mr.setFileName( argc[1]);
	//mr.read( M );

	//cout << M.shape() << endl;
	
	nbfVTKInterface converter;

	nbfMatlabWriter mwriter;
	mwriter.setFileName("ipath");

	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName(argc[1]);
	reader->Update();

	int dims[3];
	reader->GetOutput()->GetDimensions(dims);
	cout << dims[0] << ", " << dims[1] << ", " << dims[2] << endl;

	vtkImageData * image = vtkImageData :: New();
	image->DeepCopy( reader->GetOutput() );
	//converter.blitzToVtk( M, image );

	reader->SetFileName(argc[2]);
	reader->Update();	

	// apply wedge filter
	nbfWedge< float > wedge;
	wedge.setWedge(50);
	wedge.setGeometry( reader->GetOutput() );
	vtkTransform * t1 = vtkTransform::New();
	//t1->RotateZ( -90 );
	vtkTransform * t2 = vtkTransform::New();
	t2->RotateX( 70 );
	//wedge.composeWedges(t1,t2);
	//vtkImageData * filtered = vtkImageData::New();
	//wedge.applyWedgeReal( image, filtered );
	//image->DeepCopy( filtered );
	//wedge.applyWedgeReal( reader->GetOutput(), filtered );

	nbfAngularSearchCorrelationMetric< float > ccc;
	ccc.setWedge( &wedge );
	ccc.setAngleSearchX(0,0);
	ccc.setAngleSearchY(0,0);
	ccc.setAngleSearchZ(0,0);
	ccc.setInput1( image );
	ccc.setInput2( reader->GetOutput(), t2 );
	ccc.execute();
	cout << ccc.getCorrelationPeak() << endl;
	cout << ccc.getDistance() << endl;

	//Array< float, 2 > A;
	//converter.vtkToBlitz( ccc.getCorrelationMatrix(), A );

	//vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	//writer->SetFileName( argc[2] );
	//writer->SetInput( ccc.getCorrelationMatrix() );
	//writer->Write();

//	mwriter.write(A);

}