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

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>

#include <nbfTimer.h>
#include <em/nbfDensityLocalMaximaClustering.h>
#include <em/nbfAngularSearchCorrelationMetric.h>
#include <em/nbfIterativeAveraging.h>

#include <vtkPNGReader.h>

#include <random/normal.h>

#define PIXEL float

void main( int argv, char ** argc )
{
	//vtkPNGReader * reader = vtkPNGReader::New();

	//nbfVTKInterface converter;
	//Array< PIXEL, 2 > A;

	//vector< vtkImageData * > input;

	//for ( int i = 1; i <= atoi(argc[2]); i++ ){
	//	stringstream fileName;
	//	if ( i < 10 ){
	//		fileName << argc[1] << "_0" << i << ".png";
	//	}
	//	else{
	//		fileName << argc[1] << "_" << i << ".png";
	//	}
	//	cout << fileName.str() << endl;
	//	reader->SetFileName(fileName.str().c_str());
	//	reader->Update();
	//	vtkImageData * image = vtkImageData::New();
	//	image->DeepCopy( reader->GetOutput() );
	//	input.push_back( image );
	//}

	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName( argc[1] );
	reader->Update();

	vtkImageData * image = vtkImageData::New();
	image->DeepCopy( reader->GetOutput() );

	//reader->SetFileName( argc[2] );
	//reader->Update();

	//vtkImageFFT * fft = vtkImageFFT::New();
	//fft->SetDimensionality(3);
	//fft->SetInput( reader->GetOutput() );
	//fft->Update();
	//vtkImageRFFT * rfft = vtkImageRFFT::New();
	//rfft->SetDimensionality(3);
	//rfft->SetInput( fft->GetOutput() );
	//rfft->Update();

	//vtkImageExtractComponents * real = vtkImageExtractComponents::New();
	//real->SetInput( rfft->GetOutput() );
	//real->SetComponents(0);
	//real->Update();
	//
	//vtkImageMathematics * subs = vtkImageMathematics::New();
	//subs->SetOperationToSubtract();
	//subs->SetInput1( real->GetOutput() );
	//subs->SetInput2( reader->GetOutput() );
	//subs->Update();

	//vtkImageMathematics * abs = vtkImageMathematics::New();
	//abs->SetOperationToAbsoluteValue();
	//abs->SetInput1( subs->GetOutput() );
	//abs->Update();

	int dims[3];
	reader->GetOutput()->GetDimensions(dims);
	cout << dims[0] << ", " << dims[1] << ", " << dims[2] << endl;
	
	vtkImageChangeInformation * center = vtkImageChangeInformation::New();
	center->CenterImageOn();
	center->SetInput( reader->GetOutput() );
	center->Update();

	ranlib::Uniform< PIXEL > normalGen;
	normalGen.seed((unsigned int)time(0));

	vector< vtkImageData * > input;

	vtkImageReslice * reslice = vtkImageReslice::New();
	reslice->WrapOff();
	reslice->MirrorOff();
	reslice->SetBackgroundLevel(0);
	reslice->SetInput( center->GetOutput() );
	reslice->SetInterpolationModeToCubic();

	vtkTransform * transform = vtkTransform::New();
	transform->RotateX( 0 );
	transform->RotateY( 0 );
	transform->RotateZ( 0.25 );
	transform->Translate(0,.0,.0);
	reslice->SetResliceTransform( transform );
	reslice->Update();

	//// apply wedge filter
	//nbfWedge< double > wedge;
	//wedge.setGeometry( image );
	//vtkTransform * t1 = vtkTransform::New();
	//t1->RotateX( 0 );
	//wedge.wedge1On(50,t1);
	//vtkImageData * filtered1 = vtkImageData::New();
	//wedge.applyReal( image, filtered1 );

	//vtkTransform * t2 = vtkTransform::New();
	//t2->RotateX( 30 );
	//wedge.wedge1On(50,t2);
	//vtkImageData * filtered2 = vtkImageData::New();
	//wedge.applyReal( reslice->GetOutput(), filtered2 );

	//wedge.wedge1On(50,t1);
	//wedge.wedge2On(50,t2);

	nbfWedge< double > w;
	w.setGeometry( image );
	//w.lowPassOn(.25,.25,.25);
	//w.highPassOn(.15,.15,.15);

	//vtkImageMedian3D * median = vtkImageMedian3D::New();
	//median->SetKernelSize(5,5,5);
	//median->SetInput( image );
	//median->Update();
	//image->DeepCopy( median->GetOutput() );
	//median->SetInput( reslice->GetOutput() );
	//median->Update();
	//reslice->GetOutput()->DeepCopy( median->GetOutput() );

	nbfAngularSearchCorrelationMetric< double > metric( &w );
	metric.setAngleSearchX(0,0,1);
	metric.setAngleSearchY(0,0,1);
	metric.setAngleSearchZ(-1,1,1);
	metric.setInput1( image );
	metric.setInput2( reslice->GetOutput() );
	metric.setDownSampling(2);
	nbfTimer t;
	t.start();
	cout << metric.getDistance() << endl;
	t.stop();
	cout << "t = " << t.elapsedSeconds() << endl;
	cout << metric.getCorrelationPeak() << endl;

	vtkTransform * trans = metric.getTransform();
	cout << trans->GetOrientation()[0] << ", " << trans->GetOrientation()[1] << ", " << trans->GetOrientation()[2] << endl;
	metric.setAngleSearchRefineX( .5, .5 );
	metric.setAngleSearchRefineY( .5, .5 );
	metric.setAngleSearchRefineZ( .5, .5 );
	metric.setDownSampling(1);
	t.start();
	cout << metric.getDistance() << endl;
	t.stop();
	cout << "t = " << t.elapsedSeconds() << endl;
	cout << metric.getCorrelationPeak() << endl;

	metric.setAngleSearchRefineX( .25, .25 );
	metric.setAngleSearchRefineY( .25, .25 );
	metric.setAngleSearchRefineZ( .25, .25 );
	metric.setDownSampling(1);
	t.start();
	cout << metric.getDistance() << endl;
	t.stop();
	cout << "t = " << t.elapsedSeconds() << endl;
	cout << metric.getCorrelationPeak() << endl;

	//trans = metric.getTransform();
	//metric.setAngleSearchX( trans->GetOrientation()[0] - 1, trans->GetOrientation()[0] + 1, .5 );
	//metric.setAngleSearchY( trans->GetOrientation()[1] - 1, trans->GetOrientation()[1] + 1, .5 );
	//metric.setAngleSearchZ( trans->GetOrientation()[2] - 1, trans->GetOrientation()[2] + 1, .5 );
	//metric.setDownSampling(1);
	//t.start();
	//cout << metric.getDistance() << endl;
	//t.stop();
	//cout << "t = " << t.elapsedSeconds() << endl;
	//cout << metric.getCorrelationPeak() << endl;

	//return;

	//vtkImageFFT * fft = vtkImageFFT::New();
	//fft->SetDimensionality(3);
	//fft->SetInput( reader->GetOutput() );
	//fft->Update();
	//vtkImageButterworthHighPass * filter = vtkImageButterworthHighPass::New();
	//filter->SetCutOff(.1,.1,.1);
	//filter->SetOrder(2);
	//filter->SetInput( fft->GetOutput() );
	//filter->Update();
	//vtkImageRFFT * rfft = vtkImageRFFT::New();
	//rfft->SetDimensionality(3);
	//rfft->SetInput( filter->GetOutput() );
	//rfft->Update();
	//vtkImageExtractComponents * real = vtkImageExtractComponents::New();
	//real->SetInput( rfft->GetOutput() );
	//real->SetComponents(0);
	//real->Update();
	//filtered->DeepCopy( real->GetOutput() );


	//vtkTransform * t1 = vtkTransform::New();

	//// add first volume
	//input.push_back( reader->GetOutput() );
	//for ( int i = 0; i < atoi( argc[2] ); i++ ){
	//	vtkTransform * transform = vtkTransform::New();
	//	int angle = floor( ( normalGen.random() - .5 ) * 10 );
	//	cout << angle << endl;
	//	transform->RotateZ( angle );
	//	reslice->SetResliceTransform( transform );
	//	reslice->Update();
	//	vtkImageData * current = vtkImageData::New();
	//	current->DeepCopy( reslice->GetOutput() );
	//	input.push_back( current );
	//}

	//nbfDensityLocalMaximaClustering< PIXEL > cluster;
	//cluster.setInput( input );
	//nbfAngularSearchCorrelationMetric< PIXEL > metric;
	//metric.setAngleSearchZ(-10,10,1);

	//vtkImageButterworthHighPass * highpass = vtkImageButterworthHighPass::New();
	//highpass->SetCutOff(.25,.25,.25);
	////metric.setHighPassFilter( highpass );

	//vtkImageButterworthLowPass * lowpass = vtkImageButterworthLowPass::New();
	//highpass->SetCutOff(.75,.75,.75);
	////metric.setLowPassFilter( lowpass );

	//cluster.setMetric( &metric );
	//cluster.setDistanceRadius(.12);

	////metric.setInput1( reader->GetOutput() );
	////metric.setInput2( reslice->GetOutput() );
	////cout << metric.getDistance() << endl;

	//nbfTimer t;
	//t.start();
	//// vector< vector< int > > classes = cluster.execute();
	//t.stop();
	//cout << "running time = " << t.elapsedSeconds() << endl;

	////cout << "clusters = " << classes.size() << endl;

	////for ( int i = 0; i < classes.size(); i++ ){
	////	cout << "class " << i << ", elements = ";
	////	for ( int j = 0; j < classes[i].size(); j++ ){
	////		cout << classes[i][j] << ", ";
	////	}
	////	cout << endl;
	////}

	////vector< vtkImageData * > oneCluster;
	////for ( int i = 0; i < classes[0].size(); i++ ){
	////	oneCluster.push_back( input[ classes[0][i] ] );
	////}

	//nbfIterativeAveraging< PIXEL > average;
	//average.setInput( input );
	//average.setMetric( &metric );

	//vtkImageData * averageImage = vtkImageData::New();
	//t.start();
	//average.execute(1, averageImage);
	//t.stop();
	//cout << "running time = " << t.elapsedSeconds() << endl;

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetInput( reslice->GetOutput() );
	writer->SetFileName( argc[3] );
	writer->Write();
	writer->Delete();

	//center->Delete();
	//reslice->Delete();

}