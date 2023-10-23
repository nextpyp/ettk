#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfImageWriter.h>
#include <nbfMaximalFlow.h>
#include <bs/nbfBordStrategyMirror.h>

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

#define PIXEL float

void main( int argv, char ** argc )
{
	vtkStructuredPointsReader * readerV = vtkStructuredPointsReader::New();
	readerV->SetFileName( argc[1] );
	readerV->Update();

	vtkImageData * g = vtkImageData::New();
	g->DeepCopy( readerV->GetOutput() );

	nbfVTKInterface converter;
	Array< short, 3 > imageShort;
	converter.vtkToBlitz( g, imageShort );

	Array< PIXEL, 3 > image;
	image.resize( imageShort.shape() );
	image = cast<PIXEL>( imageShort );
	cout << image.shape() << endl;

	//reader.read(image);
	image = image - min(image);
	image = image / max(image) + 1e-1;

	cout << image.shape() << endl;
	cout << min(image) << ", " << max(image) << endl;

	Timer t;

	//Array< PIXEL, 2 > P( image.shape() );
	//P = 0;
	//P( 0, Range::all() ) = -1;
	//P( image.ubound(0), Range::all() ) = -1;
	//P( Range::all(), 0 ) = -1;
	//P( Range::all(), image.ubound(1) ) = -1;
	//P( 80, 30 ) = 1;

	Array< PIXEL, 3 > P( image.shape() );
	P = 0;
	P( 0, Range::all(), Range::all() ) = -1;
	P( P.ubound(0), Range::all(), Range::all() ) = -1;
	P( Range::all(), 0, Range::all() ) = -1;
	P( Range::all(), P.ubound(1), Range::all() ) = -1;
	P( Range::all(), Range::all(), 0 ) = -1;
	P( Range::all(), Range::all(), P.ubound(2) ) = -1;
	
	Range I( P.rows() / 2 - 12, P.rows() / 2 + 12 );
	Range J( P.cols() / 2 - 12, P.cols() / 2 + 12 );
	Range K( P.depth() / 2 - 3, P.depth() / 2 + 3 );
	P( I, J, K ) = 1;

	vtkImageData * weighted = vtkImageData::New();
	nbfMaximalFlow< PIXEL > flow;

	t.start();
	flow.execute(P,image,atoi(argc[3]));
	t.stop();
	cout << "Time = " << t.elapsedSeconds() << endl;

	converter.blitzToVtk( P, weighted );
	
	vtkStructuredPointsWriter * vwriter = vtkStructuredPointsWriter::New();
	vwriter->SetFileName( argc[2] );
	vwriter->SetInput( weighted );
	vwriter->Write();
}