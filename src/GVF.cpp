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

	Array< float, 2 > d1;
	reader.setFileName( argc[3] );
	reader.read( d1 );

	Array< float, 2 > d2;
	reader.setFileName( argc[4] );
	reader.read( d2 );

	Array< float, 2 > d3;
	reader.setFileName( argc[5] );
	reader.read( d3 );

	Array< float, 2 > d4;
	reader.setFileName( argc[6] );
	reader.read( d4 );

	BordStrategyMirrorDouble< float, 2 > BSforA( A, 2 );

	// tmp arrays for computations
	Array< float, 2 > S( A.shape() ), R( A.shape() );

	Timer tEvolution, tReinit;
	tEvolution.start();
	double totalReinit = 0;

	Array< float, 2 > d1x( A.shape() ), d1y( A.shape() );
	Array< float, 2 > d2x( A.shape() ), d2y( A.shape() );
	Array< float, 2 > d3x( A.shape() ), d3y( A.shape() );
	Array< float, 2 > d4x( A.shape() ), d4y( A.shape() );

	BordStrategyMirrorDouble< float, 2 > BSford1( d1, 2 );
	d1x = central12n(d1,firstDim);
	d1y = central12n(d1,secondDim);

	// normalization
	d1 = sqrt( pow2(d1x) + pow2(d1y) );
	d1x = where( d1 > 0, d1x / d1, d1x );
	d1y = where( d1 > 0, d1y / d1, d1y );

	BordStrategyMirrorDouble< float, 2 > BSford2( d2, 2 );
	d2x = central12n(d2,firstDim);
	d2y = central12n(d2,secondDim);
	
	// normalization
	d2 = sqrt( pow2(d2x) + pow2(d2y) );
	d2x = where( d2 > 0, d2x / d2, d2x );
	d2y = where( d2 > 0, d2y / d2, d2y );

	BordStrategyMirrorDouble< float, 2 > BSford3( d3, 2 );
	d3x = central12n(d3,firstDim);
	d3y = central12n(d3,secondDim);

	// normalization
	d3 = sqrt( pow2(d3x) + pow2(d3y) );
	d3x = where( d3 > 0, d3x / d3, d3x );
	d3y = where( d3 > 0, d3y / d3, d3y );

	BordStrategyMirrorDouble< float, 2 > BSford4( d4, 2 );
	d4x = central12n(d4,firstDim);
	d4y = central12n(d4,secondDim);

	// normalization
	d4 = sqrt( pow2(d4x) + pow2(d4y) );
	d4x = where( d4 > 0, d4x / d4, d4x );
	d4y = where( d4 > 0, d4y / d4, d4y );

	// parameters
	float timeStep = atof( argc[7] );	// evolution time step
	int iterations = atoi( argc[8] );	// evolution iterations
	int gapMorel = atoi( argc[9] );		// re-initialize every number of iterations
	float stopMorel = atof( argc[10] );	// stop condition for morel PDE
	float stepMorel = atof( argc[11] );	// morel PDE time step

	Array< float, 2 > F1( A.shape() ), F2( A.shape() ), F3( A.shape() ), F4( A.shape() );
	
	for ( int i = 0; i < iterations; i++ ){
	
		F1 = advection2D(A,d1x,d1y);
		F2 = advection2D(A,d2x,d2y);
		F3 = advection2D(A,d3x,d3y);
		F4 = advection2D(A,d4x,d4y);
		F1 = where( F1 >= F2, F1, F2 );
		F2 = where( F1 >= F3, F1, F3 );
		F3 = where( F2 >= F4, F2, F4 );

		//R = 1.0 * curvature2D(A) + 1.0 * where( F1 >= F2, F1, F2 );
		R = 1.0 * curvature2D(A) + 1.0 * F3;
//		R = 1.0 * curvature2D(A) + 0.25 * F1;
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
	writer.setFileName( argc[12] );
	writer.write(A);	
}
