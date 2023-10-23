#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfImageWriter.h>
#include <io/nbfMrcReader.h>
#include <io/nbfVTKInterface.h>
#include <nbfMaximalFlow.h>
#include <vtkImageCast.h>

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

#include <vtkImageContinuousDilate3D.h>
#include <vtkImageThreshold.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>

#define PIXEL float

void main( int argv, char ** argc )
{
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	
	// read image volume
	reader->SetFileName( argc[1] );
	reader->Update();

	vtkImageData * g = vtkImageData::New();
	g->DeepCopy( reader->GetOutput() );

	nbfVTKInterface converter;
	Array< short, 3 > imageShort;
	converter.vtkToBlitz( g, imageShort );

	Array< PIXEL, 3 > G;
	G.resize( imageShort.shape() );
	G = cast<PIXEL>( imageShort );
	cout << G.shape() << endl;

	G = where( G < 586, 586, G );
	G = where( G > 794, 794, G );

	G = G - min(G);
	G = G / max(G) + 1e-1;

	cout << G.shape() << endl;
	cout << min(G) << ", " << max(G) << endl;

	// read source geometry
	reader->SetFileName( argc[2] );
	reader->Update();

	Array< PIXEL, 3 > P;
	converter.vtkToBlitz( reader->GetOutput(), P );
	cout << P.shape() << endl;

	// use the inside of a tube of given radius
	P = where( P < 20, 1, 0 );
	// P( Range::all(), Range::all(), Range(P.ubound(2)-2,P.ubound(2)) ) = -1;
	
	nbfMaximalFlow< PIXEL > flow;

	Range I(0,P.ubound(0),2);
	Range J(0,P.ubound(1),2);
	Range K(0,P.ubound(2),2);
	Array< PIXEL, 3 > Psub = P(I,J,K);
	Array< PIXEL, 3 > Gsub = G(I,J,K);

	// read outer membrane
	reader->SetFileName( argc[3] );
	reader->Update();
	vtkImageData * p = vtkImageData::New();
	p->DeepCopy( reader->GetOutput() );

	vtkImageThreshold * thres = vtkImageThreshold::New();
	thres->SetInput( p );
	thres->ThresholdByLower(3.0);
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
	Psub = where( P1 == 1, -1, Psub );

	//Gsub = where( S < 5, -1, Gsub );

	//nbfMatlabWriter w;
	//w.setFileName("test.blitz");
	//w.write(Psub);

	//return;

	int pad = 0;
	Psub( Range(0,pad), Range::all(), Range::all() ) = -1;
	Psub( Range(Psub.ubound(0)-pad,Psub.ubound(0)), Range::all(), Range::all() ) = -1;
	Psub( Range::all(), Range(0,pad), Range::all() ) = -1;
	Psub( Range::all(), Range(Psub.ubound(1)-pad,Psub.ubound(1)), Range::all() ) = -1;
	Psub( Range::all(), Range::all(), Range(0,pad) ) = -1;

	cout << "Psub = " << Psub.shape() << endl;
	cout << "Psub = " << min(Psub) << ", " << max(Psub) << endl;

	Timer t;

	t.start();
	flow.execute(Psub,Gsub,atoi(argc[4]));
	t.stop();
	cout << "Time = " << t.elapsedSeconds() << endl;

	Array< PIXEL, 3 > tmp( Psub.shape() );
	tmp = Psub;
	vtkImageData * weighted = vtkImageData::New();
	converter.blitzToVtk( tmp, weighted );

	vtkStructuredPointsWriter * vwriter = vtkStructuredPointsWriter::New();
	vwriter->SetFileName( argc[5] );
	vwriter->SetInput( weighted );
	vwriter->Write();
}