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
#define DIM 3

void main( int argv, char ** argc ){

	Array< float, DIM > A;
	nbfMatlabReader reader;
	reader.setFileName( argc[1] );
	reader.read( A );

	Array< float, DIM > G;
	reader.setFileName( argc[2] );
	reader.read( G );

	BordStrategyMirrorDouble< float, DIM > BSforA( A, 2 );

	// tmp arrays for computations
	Array< float, DIM > S( A.shape() ), R( A.shape() );

	// Narrow band
	//  std :: vector< TinyVector<int,2> > I;
	//	float th = 3.0;
	//	nbfNarrowBand :: build( A, I, th );
	//	cout << I.size() << " of " << A.rows() * A.cols() << endl;

	Timer tEvolution, tReinit;
	tEvolution.start();
	double totalReinit = 0;

	// parameters
	float timeStep = atof( argc[3] );	// evolution time step
	int iterations = atoi( argc[4] );	// evolution iterations
	int gapMorel = atoi( argc[5] );		// re-initialize every number of iterations
	float stopMorel = atof( argc[6] );	// stop condition for morel PDE
	float stepMorel = atof( argc[7] );	// morel PDE time step

	for ( int i = 0; i < iterations; i++ ){
	
		R = G * curvature3D(A);
		A = A + timeStep * R;
		BSforA.refresh();

		// reinitialization
//		if ( i > 0 & i % gapMorel == 0 ){
		if ( i % gapMorel == 0 ){

			tReinit.start();

			// compute distances
			S = sussmanSign( A );
			do {
				R = morelSussman3D(A,S);
				A = A + stepMorel * R;
				BSforA.refresh();
			} while( max(abs(R)) > stopMorel );
#if 0
			// get maximum distance in narrowband
			S = 0;
			S[I] = abs(A);
			float limit = max(S);

			// if too big recompute narrowband
			if ( limit > 1.5 * th ){
				cout << "Recomputing narrowband, limit = " << limit << endl;

				// compute distances
				S = sussmanSign( A );
				do {
					R = morelSussman2D(A,S);
					A = A + 1e-2 * R;
					BSforA.refresh();
				} while( max(abs(R)) > stopMorel );

				// rebuild narrowband
				nbfNarrowBand :: build( A, I, th );
			}
			else {
				// only update values in narrow band
                S[I] = sussmanSign( A );
				R = 0;
				do {
					R[I] = morelSussman2D(A,S);
					A[I] = A + 1e-2 * R;
				} while( max(abs(R)) > stopMorel );
			}
#endif
			tReinit.stop();
			totalReinit += tReinit.elapsedSeconds();
		}
	}

	tEvolution.stop();
	cout << "Tiempo total del algoritmo (secs)= " << tEvolution.elapsedSeconds() << endl;
	cout << "Time reinitializing (secs)= " << totalReinit << " (" << totalReinit / tEvolution.elapsedSeconds() * 100 << "%)" << endl;

	nbfMatlabWriter writer;
	writer.setFileName( argc[8] );
	writer.write(A);	
}
