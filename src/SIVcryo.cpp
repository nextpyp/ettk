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
#include <vtkPolyDataWriter.h>
#include <vtkStructuredPoints.h>
#include <vtkContourFilter.h>
#include <vtkImageButterworthHighPass.h>
#include <vtkContourFilter.h>

#include <vtkImageContinuousDilate3D.h>
#include <vtkImageContinuousErode3D.h>
#include <vtkImageThreshold.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkImageThreshold.h>
#include <vtkImageGradientMagnitude.h>

#define PIXEL float

void main( int argv, char ** argc )
{
#if 0 // compute metric from image

	vtkImageData * g = vtkImageData::New();

	nbfMrcReader reader;
	reader.setFileName( argc[1] );
	reader.read( g );

	float factor = 2.0;
	vtkImageResample * resample = vtkImageResample::New();
	resample->SetInput( g );
	resample->SetAxisMagnificationFactor(0,1.0/factor);
	resample->SetAxisMagnificationFactor(1,1.0/factor);
	resample->SetAxisMagnificationFactor(2,1.0/factor);
	resample->SetInterpolationModeToLinear();
	resample->Update();

	g->DeepCopy( resample->GetOutput() );

	nbfVTKInterface converter;
	Array< PIXEL, 3 > G;
	converter.vtkToBlitz( g, G );

	G = where( G < 8691, 8691, G );
	G = where( G > 12276, 12276, G );
	G = G - min(G);
	G = G / max(G);
	G = G + 1e-1;

	converter.blitzToVtk( G, g );

	cout << G.shape() << endl;
	cout << min(G) << ", " << max(G) << endl;

	Array< PIXEL, 3 > P( G.shape() );
	P = 0;

	// set source geometry
	P( Range(188/factor,300/factor), Range(200/factor,344/factor), Range(110/2,120/2) ) = 1;

	// set sink geometry to volume sides (outer membrane)
	int pad = 10;
	P( Range(0,pad), Range::all(), Range::all() ) = -1;
	P( Range(P.ubound(0)-pad,P.ubound(0)), Range::all(), Range::all() ) = -1;
	P( Range::all(), Range(0,pad), Range::all() ) = -1;
	P( Range::all(), Range(P.ubound(1)-pad,P.ubound(1)), Range::all() ) = -1;
	P( Range::all(), Range::all(), P.lbound(thirdDim) ) = -1;
	P( Range::all(), Range::all(), P.ubound(thirdDim) ) = -1;
	
	nbfMaximalFlow< PIXEL > flow;

	nbfTimer t;

	t.start();
	flow.execute(P,G,atoi(argc[2]));
	t.stop();
	cout << "Time = " << t.elapsedSeconds() << endl;

	Array< PIXEL, 3 > tmp( P.shape() );
	tmp = P;
	vtkImageData * weighted = vtkImageData::New();
	converter.blitzToVtk( tmp, weighted );

	converter.blitzToVtk( P, weighted );

	resample->SetInput( weighted );
	resample->SetAxisMagnificationFactor(0,factor);
	resample->SetAxisMagnificationFactor(1,factor);
	resample->SetAxisMagnificationFactor(2,factor);
	resample->Update();

#elif 0
	// inner core

	// read input image
	vtkImageData * g = vtkImageData::New();

	nbfMrcReader reader;
	reader.setFileName( argc[1] );
	reader.read( g );
	
	float factor = 2.0;
	vtkImageResample * resample = vtkImageResample::New();
	resample->SetInput( g );
	resample->SetAxisMagnificationFactor(0,1.0/factor);
	resample->SetAxisMagnificationFactor(1,1.0/factor);
	resample->SetAxisMagnificationFactor(2,1.0/factor);
	resample->SetInterpolationModeToLinear();
	resample->Update();

	g->DeepCopy( resample->GetOutput() );

	nbfVTKInterface converter;
	Array< PIXEL, 3 > G;
	converter.vtkToBlitz( g, G );

	G = where( G < 8691, 8691, G );
	G = where( G > 12276, 12276, G );
	G = G - min(G);
	G = G / max(G);
	G = G + 1e-1;

	cout << G.shape() << endl;

	// presure storage
	Array< PIXEL, 3 > P( G.shape() );
	P = 0;

	// read source from file
	vtkStructuredPointsReader * vreader = vtkStructuredPointsReader::New();	
	vreader->SetFileName( argc[2] );
	vreader->Update();

	resample->SetInput( vreader->GetOutput() );
	resample->Update();

	Array< PIXEL, 3 > S;
	nbfVTKInterface::vtkToBlitz( resample->GetOutput(), S );

	cout << S.shape() << endl;

	P = where( S < 27, 1, P );

	cout << "P = " << min(P) << ", " << max(P) << endl;

	cout << "source done." << endl;

	// read membrane and use as sink
	vtkStructuredPointsReader * vreader2 = vtkStructuredPointsReader::New();	
	vreader2->SetFileName( argc[3] );
	vreader2->Update();

	vtkImageResample * resample2 = vtkImageResample::New();
	resample2->SetInput( vreader2->GetOutput() );
	resample2->SetAxisMagnificationFactor(0,1.0/factor);
	resample2->SetAxisMagnificationFactor(1,1.0/factor);
	resample2->SetAxisMagnificationFactor(2,1.0/factor);
	resample2->SetInterpolationModeToLinear();
	resample2->Update();

	vtkImageContinuousErode3D * dilate = vtkImageContinuousErode3D::New();
	dilate->SetInput( resample2->GetOutput() );
	dilate->SetKernelSize(7,7,7);
	dilate->Update();

	Array< PIXEL, 3 > S1;
	nbfVTKInterface::vtkToBlitz( dilate->GetOutput(), S1 );

	cout << S1.shape() << endl;

	// set sink to membrane
	P = where( S1 < 1, -1, P );

	cout << "P = " << min(P) << ", " << max(P) << endl;
	cout << "sink done." << endl;

	nbfMaximalFlow< PIXEL > flow;

	nbfTimer t;

	t.start();
	//flow.execute(P,G,atoi(argc[4]));
	t.stop();
	cout << "Time = " << t.elapsedSeconds() << endl;

	Array< PIXEL, 3 > tmp( P.shape() );
	tmp = P;
	vtkImageData * weighted = vtkImageData::New();
	converter.blitzToVtk( tmp, weighted );

#elif 0
	// inner core

	// read input image
	vtkImageData * g = vtkImageData::New();

	nbfMrcReader reader;
	reader.setFileName( argc[1] );
	reader.read( g );
	
	float factor = 2.0;
	vtkImageResample * resample = vtkImageResample::New();
	resample->SetInput( g );
	resample->SetAxisMagnificationFactor(0,1.0/factor);
	resample->SetAxisMagnificationFactor(1,1.0/factor);
	resample->SetAxisMagnificationFactor(2,1.0/factor);
	resample->SetInterpolationModeToLinear();
	resample->Update();

	g->DeepCopy( resample->GetOutput() );

	nbfVTKInterface converter;
	Array< PIXEL, 3 > G;
	converter.vtkToBlitz( g, G );

	G = where( G < 8691, 8691, G );
	G = where( G > 12276, 12276, G );
	G = G - min(G);
	G = G / max(G);
	G = G + 1e-1;

	cout << G.shape() << endl;

	// read presure from file
	vtkStructuredPointsReader * vreader = vtkStructuredPointsReader::New();	
	vreader->SetFileName( argc[2] );
	vreader->Update();

	Array< PIXEL, 3 > P;
	nbfVTKInterface::vtkToBlitz( vreader->GetOutput(), P );

	cout << P.shape() << endl;
	cout << "P = " << min(P) << ", " << max(P) << endl;
	P( Range(0,154/2), Range(355/2,P.ubound(secondDim)), Range(40/2,200/2) ) = -1;

	// modify source
	P( Range(295/2,316/2), Range(270/2,318/2), Range(115/2,125/2) ) = 1;

	nbfMaximalFlow< PIXEL > flow;

	nbfTimer t;

	t.start();
	flow.execute(P,G,atoi(argc[3]));
	t.stop();
	cout << "Time = " << t.elapsedSeconds() << endl;

	Array< PIXEL, 3 > tmp( P.shape() );
	tmp = P;
	vtkImageData * weighted = vtkImageData::New();
	converter.blitzToVtk( tmp, weighted );

	resample->SetInput( weighted );
	resample->SetAxisMagnificationFactor(0,factor);
	resample->SetAxisMagnificationFactor(1,factor);
	resample->SetAxisMagnificationFactor(2,factor);
	resample->Update();

	weighted->DeepCopy( resample->GetOutput() );

#elif 0

	// segment spikes by threshold in a band outside the membrane
	// read membrane
	// read input image
	vtkImageData * g = vtkImageData::New();

	nbfMrcReader reader;
	reader.setFileName( argc[1] );
	reader.read( g );

	// read membrane
	vtkStructuredPointsReader * vreader = vtkStructuredPointsReader::New();	
	vreader->SetFileName( argc[2] );
	vreader->Update();

	float factor = 4.0;
	vtkImageResample * vresample = vtkImageResample::New();
	vresample->SetInput( vreader->GetOutput() );
	vresample->SetAxisMagnificationFactor(0,1.0/factor);
	vresample->SetAxisMagnificationFactor(1,1.0/factor);
	vresample->SetAxisMagnificationFactor(2,1.0/factor);
	vresample->SetInterpolationModeToLinear();
	vresample->Update();

	Array< PIXEL, 3 > L;
	nbfVTKInterface::vtkToBlitz( vresample->GetOutput(), L );
	Array< PIXEL, 3 > Lcopy( L.shape() );
	Lcopy = L;

	vtkImageContinuousDilate3D * dilate = vtkImageContinuousDilate3D::New();
	dilate->SetInput( vresample->GetOutput() );
	dilate->SetKernelSize(15,15,15);
	dilate->Update();

	Array< PIXEL, 3 > Ldilate;
	nbfVTKInterface::vtkToBlitzReference( dilate->GetOutput(), Ldilate );
	Ldilate = where( ( Lcopy < 1 ) && ( Ldilate > 1 ), 1, 0 ); // this selects the outside

	// upsample
	vtkImageResample * resample = vtkImageResample::New();
	resample->SetInterpolationModeToLinear();
	resample->SetInput( dilate->GetOutput() );
	resample->SetAxisMagnificationFactor(0,factor);
	resample->SetAxisMagnificationFactor(1,factor);
	resample->SetAxisMagnificationFactor(2,factor);
	resample->Modified();
	resample->Update();

	Array< PIXEL, 3 > region;
	nbfVTKInterface::vtkToBlitz( resample->GetOutput(), region );

	cout << region.shape() << endl;

	Array< PIXEL, 3 > input;
	nbfVTKInterface::vtkToBlitz( g, input );

	cout << input.shape() << endl;
	Array< PIXEL, 3 > vistaInput( input( Range::all(), Range::all(), Range(0,input.depth()-2) ) );
	cout << vistaInput.shape() << endl;
	
	vistaInput = where( region > 0, vistaInput, max(vistaInput) );
	vistaInput( Range::all(), Range::all(), Range(0,70) ) = max(vistaInput);
	vistaInput( Range::all(), Range::all(), Range(150,vistaInput.ubound(2)) ) = max(vistaInput);

	vtkImageData * weighted = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( vistaInput, weighted );

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToShort();
	cast->SetInput( weighted );
	cast->Update();

	weighted->DeepCopy( cast->GetOutput() );

#else

	vtkImageData * g = vtkImageData::New();
	nbfMrcReader reader;
	reader.setFileName( argc[1] );
	reader.read( g );

	//// read image volume
	//vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();	
	//reader->SetFileName( argc[1] );
	//reader->Update();

	//Array< double, 3 > A;
	//nbfVTKInterface::vtkToBlitzReference( reader->GetOutput(), A );
	//A = where( A < 3, 0 , 1 );

	//vtkImageGaussianSmooth * smooth = vtkImageGaussianSmooth::New();
	//smooth->SetInput( reader->GetOutput() );
	//smooth->SetDimensionality(3);
	//smooth->SetRadiusFactors(10.0,10.0,10.0);
	//smooth->Update();

	//vtkContourFilter * contour = vtkContourFilter::New();
	//contour->SetValue(0,.5);
	//contour->SetInput( smooth->GetOutput() );
	//contour->ComputeNormalsOn();
	//contour->Update();

	//vtkPolyDataWriter * w = vtkPolyDataWriter::New();
	//w->SetFileName( argc[4] );
	//w->SetInput( contour->GetOutput() );
	//w->Write();
	//return;

#endif

	vtkStructuredPointsWriter * vwriter = vtkStructuredPointsWriter::New();
	vwriter->SetFileName( argc[4] );
	vwriter->SetInput( g );
	vwriter->Write();
}