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
	if ( argc != 15 ){
		cout << "Usage: input1 input2 output" << endl;
		exit(0);
	}

	char * inFile1 = argv[1];
	char * inFile2 = argv[2];
	double wedgeL1 = atof(argv[3]);
	double wedgeU1 = atof(argv[4]);
	double wedgeAngle1X = atof(argv[5]);
	double wedgeAngle1Y = atof(argv[6]);
	double wedgeAngle1Z = atof(argv[7]);
	double wedgeL2 = atof(argv[8]);
	double wedgeU2 = atof(argv[9]);
	double wedgeAngle2X = atof(argv[10]);
	double wedgeAngle2Y = atof(argv[11]);
	double wedgeAngle2Z = atof(argv[12]);
	char * outFile = argv[13];
	char * outFileFourier = argv[14];

	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();

	// read 3D image data
	reader->SetFileName( inFile1 );
	reader->Update();

	///////////////////

	vtkImageChangeInformation * center = vtkImageChangeInformation::New();
	center->CenterImageOn();
	center->SetInput( reader->GetOutput() );
	center->Update();

	vtkImageReslice * reslice = vtkImageReslice::New();
	//reslice->WrapOff();
	//reslice->MirrorOff();
	//reslice->SetBackgroundLevel(0);
	reslice->SetInput( center->GetOutput() );
	reslice->SetInterpolationModeToCubic();

	vtkTransform * transform = vtkTransform::New();
	//transform->RotateZ( 45 );
	transform->RotateZ( 53.4 );		// -gamma 0
	transform->RotateY( -0  );		// -beta 37.265
	transform->RotateZ( -0 );		// -alpha 0
	//transform->RotateZ(180);
	//transform->RotateY(180);

	reslice->SetResliceTransform( transform );
	reslice->Update();

	// apply in-place window of image data
	nbfImageMetric< PIXEL, 3 > window;
	vtkImageData * data = vtkImageData::New();
	data->DeepCopy( reader->GetOutput() );
	window.applyWindow( data );

	vtkImageConstantPad * pad = vtkImageConstantPad::New();
	pad->SetConstant(0);
	int extent[6];
	reader->GetOutput()->GetExtent(extent);
	pad->SetOutputWholeExtent( extent[0], 2*extent[1]+1, extent[2], 2*extent[3]+1, extent[4], 2*extent[5]+1 );
	pad->SetInput( reader->GetOutput() );

	Array< PIXEL, 3 > R;
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), R );
	pad->SetConstant( mean(R) );
	pad->Update();

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToDouble();
	cast->SetInput( pad->GetOutput() );
	cast->Update();

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName( "pad.vtk" );
	writer->SetInput( cast->GetOutput() );
	writer->Write();

	////////////////////

	vtkImageData * image1 = vtkImageData::New();
	image1->DeepCopy( cast->GetOutput() );
	// image1->DeepCopy( pad->GetOutput() );

	// read second 3D image
	reader->SetFileName( inFile2 );
	reader->Update();
	//pad->SetInput( reslice->GetOutput() );
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), R );
	pad->SetConstant( mean(R) );
	pad->Update();

	cast->Modified();
	cast->Update();

	vtkImageData * image2 = reader->GetOutput();
	image2->DeepCopy( cast->GetOutput() );
	//image2->DeepCopy( reslice->GetOutput() );

	//image2->DeepCopy( pad->GetOutput() );

	//// define wedge filter
	//nbfWedge< double > wedge;
	//wedge.setGeometry( image1 );

	//vtkTransform * t1 = vtkTransform::New();
	//t1->RotateX( wedgeAngle1X );
	//t1->RotateY( wedgeAngle1Y );
	//t1->RotateZ( wedgeAngle1Z );
	//wedge.wedge1On(wedgeL1,wedgeU1,t1);

	//vtkTransform * t2 = vtkTransform::New();
	//t2->RotateX( wedgeAngle2X );
	//t2->RotateY( wedgeAngle2Y );
	//t2->RotateZ( wedgeAngle2Z );
	//wedge.wedge2On(wedgeL2,wedgeU2,t2);

	////wedge.lowPassOn(.2,.2,.2,4);
	////wedge.highPassOn(.16,.16,.16,4);
	//wedge.bandPassOn(.2,.5,.01);
	////wedge.bandPassOn(0.01,0.5,.01);

	//Array< double, 2 > strategy(1,2);
	//strategy = 1, 1;

	//cout << strategy << endl;

	//nbfAngularSearchCorrelationMetric< double > metric( &wedge );
	//metric.setAngleSearchX(0,0,2);
	//metric.setAngleSearchY(0,0,2);
	//metric.setAngleSearchZ(0,0,2);
	//metric.setInput1( image1 );
	//metric.setInput2( image2 );
	//metric.setDownSampling(1);
	//
	//metric.execute();
#if 0
	nbfAngularSearchCorrelationMetricFFTW< double, 3 > metric( &wedge );
	metric.setStrategy( strategy );
	metric.setAngleSearchX(0,0);
	metric.setAngleSearchY(0,0);
	metric.setAngleSearchZ(0,0);
	metric.setInput1( image1 );
	metric.setInput2( image2 );
	
	metric.execute();
	//cout << metric.getCorrelationPeak() << endl;
#else

	nbfWedgeFilter3D< PIXEL > w;
	//w.bandPassOn(.1,.5,.01);
	//w.bandPassOn(.1,1,.01);

	vtkTransform * t1 = vtkTransform::New();
	t1->RotateX( wedgeAngle1X );
	t1->RotateY( wedgeAngle1Y );
	t1->RotateZ( wedgeAngle1Z );
	w.wedge1On(wedgeL1,wedgeU1,t1);

	vtkTransform * t2 = vtkTransform::New();
	t2->RotateX( wedgeAngle2X );
	t2->RotateY( wedgeAngle2Y );
	t2->RotateZ( wedgeAngle2Z );
	w.wedge2On(wedgeL2,wedgeU2,t2);

	//w.applyReal( image1, image2 );
	//writer->SetInput( image2 );
	//writer->Write();

	nbfProjectionRotationMetric3D< PIXEL > m(&w);
	m.windowOff();
	//m.medianFilterOn();
	
	nbfAngularSearchMetric3D< PIXEL > metrika(&m);
	metrika.setAngleSearchX(0,0);
	metrika.setAngleSearchY(0,0);

	Array< PIXEL, 2 > A(1,2);
	A = 1, 1;
	cout << A << endl;

	metrika.setStrategy(A);
	metrika.setInput1( image1 );
	metrika.setInput2( image2 );
	cout << metrika.getDistance() << endl;

#endif

	//return;

	//vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	//writer->SetFileName( outFile );
	//writer->SetInput( metric.getAligned() );
	//writer->Write();
	////cout << "File: " << outFile << " written." << endl;

	//writer->SetFileName( outFileFourier );
	//writer->SetInput( metric.getAlignedFourier() );
	//writer->Write();

	//writer->Delete();
	image1->Delete();
	image2->Delete();
	reader->Delete();
}
