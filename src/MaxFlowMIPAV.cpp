#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMrcReader.h>
#include <io/nbfVTKInterface.h>
#include <nbfMaximalFlow.h>
#include <nbfDifferentials.h>

#include <vtkImageCast.h>

#include <vtkImageMathematics.h>
#include <vtkMath.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>

#include <vtkPoints.h>
#include <vtkImageData.h>

#define PIXEL float

void main( int argv, char ** argc )
{

	if ( argc != 4 ){
		cout << "USAGE: MaxFlowMIPAV.exe" << endl;
	}

	// read image volume
	nbfMrcReader mrc;
	mrc.setFileName( argc[1] );
	Array< float, 1 > a,b;

	vtkImageData * g = vtkImageData::New();
	mrc.read(g,a,b);

	Array< short, 3 > imageShort;
	nbfVTKInterface::vtkToBlitz( g, imageShort );

	Array< PIXEL, 3 > G;
	G.resize( imageShort.shape() );
	G = cast<PIXEL>( imageShort );
	cout << G.shape() << endl;
#if 0
	PIXEL m = mean(G);
	PIXEL var = 1500;
	G = where( G < m - var, m - var, G );
	G = where( G > m + var, m + var, G );
#endif
	G = G - min(G);
	G = G / max(G) + 1e-1;

	cout << G.shape() << endl;
	cout << min(G) << ", " << max(G) << endl;

	vtkImageData * source = vtkImageData::New();

	// read source geometry
	mrc.setFileName( argc[2] );
	mrc.setBigEndian(true);
	mrc.read( source, a, b );

	Array< PIXEL, 3 > sourceA;
	nbfVTKInterface::vtkToBlitz( source, sourceA );
	cout << "source = [ " << min(sourceA) << ", " << max(sourceA) << "]" << endl;

	//vtkImageData * sink = vtkImageData::New();
	//mrc.setFileName( argc[3] );
	//mrc.read( sink, a, b );

	//Array< PIXEL, 3 > sinkA;
	//nbfVTKInterface::vtkToBlitz( sink, sinkA );
	//cout << "sink = [ " << min(sinkA) << ", " << max(sinkA) << "]" << endl;

	Array< PIXEL, 3 > P( G.shape() );
	
	P = sourceA;
	//P = where( sinkA == 1, -1, P );
#if 0
	P( P.lbound(firstDim), Range::all(), Range::all() ) = -1;
	P( P.ubound(firstDim), Range::all(), Range::all() ) = -1;
	P( Range::all(), P.lbound(secondDim), Range::all() ) = -1;
	P( Range::all(), P.ubound(secondDim), Range::all() ) = -1;
#endif
	nbfMaximalFlow< PIXEL > flow;

	Range I(0,P.ubound(0),1);
	Range J(0,P.ubound(1),1);
	Range K(0,P.ubound(2),1);
	Array< PIXEL, 3 > Psub = P(I,J,K);
	Array< PIXEL, 3 > Gsub = G(I,J,K);

	nbfTimer t;

	t.start();
	flow.execute(Psub,Gsub,atoi(argc[4]));
	t.stop();
	cout << "Time = " << t.elapsedSeconds() << endl;

#if 0
	// curvature smooth
	Array< float, 3 > V( Psub.shape() );
	BordStrategyMirrorDouble< float, 3 > bsForV( V, 1 );

	// parameters
	float timeStep = .25;	// evolution time step
	int iterations = 3;	// evolution iterations
	for ( int i = 0; i < iterations; i++ ){
		V = Psub;
		bsForV.refresh();
		Psub = Psub + timeStep * curvature3D(V);
	}
#endif
	Psub = Psub - mean(Psub);

	/// store in file
	Array< PIXEL, 3 > tmp( Psub.shape() );
	tmp = Psub.reverse(secondDim);
	vtkImageData * weighted = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( tmp, weighted );

	vtkStructuredPointsWriter * vwriter = vtkStructuredPointsWriter::New();
	vwriter->SetFileName( argc[5] );
	vwriter->SetInput( weighted );
	vwriter->Write();
}