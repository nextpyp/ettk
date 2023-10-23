#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <Timer.hh>
#include <MacrosFlujos.hh>
#include <IO/FlujosIO.hh>
#include <NarrowBand/nbfMatlabReader.h>
#include <NarrowBand/nbfMatlabWriter.h>

#include <NarrowBand/nbfBordStrategyConst.h>
#include <NarrowBand/nbfBordStrategyMirror.h>
#include <NarrowBand/nbfDifferentials.h>
#include <NarrowBand/nbfNarrowBand.h>

#define PIXEL float
#define DIM 2

void main( int argv, char ** argc ){

	Array< float, 2 > A;
	nbfMatlabReader reader;
	reader.setFileName( argc[1] );
	reader.read( A );

	Array< float, 2 > G;
	reader.setFileName( argc[2] );
	reader.read( G );

	BordStrategyMirrorDouble< float, 2 > BSforA( A, 2 );

	// tmp array for computations
	Array< float, 2 > R( A.shape() );

	Timer tEvolution;
	tEvolution.start();

	// parameters
	float timeStep = atof( argc[4] );	// evolution time step
	int iterations = atoi( argc[5] );	// evolution iterations

	for ( int i = 0; i < iterations; i++ ){
		R = diffusion2D(G,A);
		A = A + timeStep * R;
		BSforA.refresh();
	}

	tEvolution.stop();
	cout << "Tiempo total del algoritmo (secs)= " << tEvolution.elapsedSeconds() << endl;

	nbfMatlabWriter writer;
	writer.setFileName( argc[3] );
	writer.write(A);	
}
