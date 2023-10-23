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

	// tmp arrays for computations
	Array< float, 2 > S( A.shape() ), R( A.shape() );

	Timer tEvolution, tReinit;
	tEvolution.start();
	double totalReinit = 0;

	Array< float, 2 > F( A.shape() ), Fx( A.shape() ), Fy( A.shape() );
	F = 1;
	Fx = central12n(G,firstDim);
	Fy = central12n(G,secondDim);

	// parameters
	float timeStep = atof( argc[4] );	// evolution time step
	int iterations = atoi( argc[5] );	// evolution iterations
	int gapMorel = atoi( argc[6] );		// re-initialize every number of iterations
	float stopMorel = atof( argc[7] );	// stop condition for morel PDE
	float stepMorel = atof( argc[8] );	// morel PDE time step

	for ( int i = 0; i < iterations; i++ ){
	
		//R = G * ( curvature2D(A) - 0.25 * expansion2D(A,F) * ( 1 - 1 / ( 1 + exp( -((double)i-4000)/200 ) ) ) ) + 1 * advection2D(A,Fx,Fy);
		//R = ( G - 1 ) * curvature2D(A) - 0.0 * expansion2D(A,F) + 1.0 * advection2D(A,Fx,Fy);
		R = G * curvature2D(A);
		A = A + timeStep * R;
		BSforA.refresh();

		// reinitialization
		if ( i % gapMorel == 0 ){

			tReinit.start();

			// compute distances
			S = sussmanSign( A );
			do {
				R = morelSussman2D(A,S);
				A = A + stepMorel * R;
				BSforA.refresh();
			} while( max(abs(R)) > stopMorel );

			tReinit.stop();
			totalReinit += tReinit.elapsedSeconds();
		}
	}

	tEvolution.stop();
	cout << "Tiempo total del algoritmo (secs)= " << tEvolution.elapsedSeconds() << endl;
	cout << "Time reinitializing (secs)= " << totalReinit << " (" << totalReinit / tEvolution.elapsedSeconds() * 100 << "%)" << endl;

	nbfMatlabWriter writer;
	writer.setFileName( argc[3] );
	writer.write(A);	
}
