#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <nbfMaximalFlow.h>
#include <io/nbfMrcReader.h>
#include <io/nbfVTKInterface.h>
#include <vtkImageCast.h>
#include <vtkImageData.h>

#include <vtkImageMathematics.h>
#include <vtkMath.h>
#include <vtkImageFFT.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageRFFT.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkImageExtractComponents.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>
#include <vtkContourFilter.h>
#include <vtkImageButterworthHighPass.h>

#include <vtkImageContinuousDilate3D.h>
#include <vtkImageThreshold.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkImageThreshold.h>
#include <vtkImageGradientMagnitude.h>

#define PIXEL float

void main( int argv, char ** argc )
{
#if 0 // high pass filtering
	nbfMrcReader reader;
	reader.setFileName( argc[1] );
	vtkImageData * image = vtkImageData::New();
	reader.read(image);

	vtkImageFFT * fft = vtkImageFFT::New();
	fft->SetDimensionality(3);
	fft->SetInput( image );
	
	vtkImageButterworthHighPass * highpass = vtkImageButterworthHighPass::New();
	highpass->SetInput( fft->GetOutput() );
	highpass->SetOrder(2);
	highpass->SetXCutOff(.01);
	highpass->SetYCutOff(.01);
	// highpass->SetZCutOff(.1);

	vtkImageRFFT * ifft = vtkImageRFFT::New();
	ifft->SetDimensionality(3);
	ifft->SetInput( highpass->GetOutput() );

	vtkImageExtractComponents * real = vtkImageExtractComponents::New();
	real->SetComponents(0);
	real->SetInput( ifft->GetOutput() );
	real->Update();
#elif 0 // compute metric from image

	// read image volume
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();	
	reader->SetFileName( argc[1] );
	reader->Update();

	vtkImageThreshold * threshold = vtkImageThreshold::New();
	threshold->ThresholdBetween(-10.0,23.0);
	threshold->SetInput( reader->GetOutput() );

	vtkImageGradientMagnitude * gradient = vtkImageGradientMagnitude::New();
	gradient->SetInput( threshold->GetOutput() );
	gradient->Update();

	vtkImageThreshold * threshold2 = vtkImageThreshold::New();
	threshold2->SetInput( gradient->GetOutput() );
	threshold2->ReplaceInOn();
	threshold2->ReplaceOutOn();
	threshold2->SetInValue(1.0);
	threshold2->SetOutValue(0.0);
	threshold2->ThresholdBetween(9,22);
	threshold2->Update();
#elif 0

	// read image volume
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();	
	reader->SetFileName( argc[1] );
	reader->Update();

	vtkImageData * g = vtkImageData::New();
	g->DeepCopy( reader->GetOutput() );

	nbfVTKInterface converter;
	Array< PIXEL, 3 > G;
	converter.vtkToBlitz( g, G );

	G = 1 - G + 1e-1;

	cout << G.shape() << endl;
	cout << min(G) << ", " << max(G) << endl;

	// read source geometry
	reader->SetFileName( argc[2] );
	reader->Update();

	Array< PIXEL, 3 > P;
	converter.vtkToBlitz( reader->GetOutput(), P );
	cout << P.shape() << endl;
	P = where( P < 5, 1, 0 );


	//// read sink geometry
	reader->SetFileName( argc[3] );
	reader->Update();

	vtkImageThreshold * thres = vtkImageThreshold::New();
	thres->SetInput( reader->GetOutput() );
	thres->ThresholdByLower(5.0);
	thres->ReplaceInOn();
	thres->ReplaceOutOn();
	thres->SetInValue(1);
	thres->SetOutValue(0);
	thres->Update();
	
	vtkImageContinuousDilate3D * dilate = vtkImageContinuousDilate3D::New();
	dilate->SetInput( thres->GetOutput() );
	dilate->SetKernelSize(5,5,5);
	dilate->Update();

	Array< PIXEL, 3 > P1;
	converter.vtkToBlitz( dilate->GetOutput(), P1 );
	cout << "P1 = " << P1.shape() << endl;
	////

	// use the inside of a tube of given radius
	P = where( P1 > 0.01, -1, P );
	
	// for outer membrane
	//P( Range(127,232), Range(68,127), P.lbound(thirdDim) ) = 1;
	
	nbfMaximalFlow< PIXEL > flow;

	//// set sink geometry to volume sides (outer membrane)
	//int pad = 0;
	//P( Range(0,50), Range::all(), Range::all() ) = -1;
	//P( Range(P.ubound(0)-pad,P.ubound(0)), Range::all(), Range::all() ) = -1;
	//P( Range::all(), Range(0,25), Range::all() ) = -1;
	//P( Range::all(), Range(169,P.ubound(1)), Range::all() ) = -1;
	//P( Range(182,P.ubound(firstDim)), Range::all(), P.ubound(thirdDim) ) = -1;

	nbfTimer t;

	t.start();
	flow.execute(P,G,atoi(argc[4]));
	t.stop();
	cout << "Time = " << t.elapsedSeconds() << endl;

	Array< PIXEL, 3 > tmp( P.shape() );
	tmp = P;
	vtkImageData * weighted = vtkImageData::New();
	converter.blitzToVtk( tmp, weighted );
#else

	// read image volume
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();	
	reader->SetFileName( argc[1] );
	reader->Update();

	vtkImageGaussianSmooth * smooth = vtkImageGaussianSmooth::New();
	smooth->SetInput( reader->GetOutput() );
	smooth->SetDimensionality(3);
	smooth->SetRadiusFactors(1.0,1.0,1.0);
	smooth->Update();

#endif

	vtkStructuredPointsWriter * vwriter = vtkStructuredPointsWriter::New();
	vwriter->SetFileName( argc[2] );
	vwriter->SetInput( smooth->GetOutput() );
	vwriter->Write();
}