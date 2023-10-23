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

	// initial phi function
	Array< float, 3 > phi;
	nbfMatlabReader reader;
	reader.setFileName( argc[1] );
	reader.read( phi );

	// driving distance function
	Array< float, 3 > G;
	reader.setFileName( argc[2] );
	reader.read( G );

	BordStrategyMirrorDouble< float, 3 > BSforPhi( phi, 2 );
	BSforPhi.refresh();

	// tmp arrays for computations
	Array< float, 3 > R( phi.shape() );

	Timer tEvolution, tReinit;
	tEvolution.start();
	double totalReinit = 0;

	Array< float, 3 > Gx( phi.shape() ), Gy( phi.shape() ), Gz( phi.shape() );

	BordStrategyMirrorDouble< float, 3 > BSforG( G, 1 );
	BSforG.refresh();

	Gx = - central12n(G,firstDim);
	Gy = - central12n(G,secondDim);
	Gz = - central12n(G,thirdDim);

	// normalization
	R = sqrt( pow2(Gx) + pow2(Gy) + pow2(Gz) );
	Gx = where( R > 0, Gx / R, Gx );
	Gy = where( R > 0, Gy / R, Gy );
	Gz = where( R > 0, Gz / R, Gz );

	// parameters
	float timeStep = atof( argc[3] );	// evolution time step
	int iterations = atoi( argc[4] );	// evolution iterations

	for ( int i = 0; i < iterations; i++ ){
	
		R = 0.5 * curvature3D(phi) - 1.0 * advection3D(phi,Gx,Gy,Gz);

		phi = phi + timeStep * R;
		
		BSforPhi.refresh();

		cout << "i = " << i << " of " << iterations << endl;
	}

	tEvolution.stop();
	cout << "Elapsed time (secs)= " << tEvolution.elapsedSeconds() << endl;

	nbfMatlabWriter writer;
	writer.setFileName( argc[5] );
	writer.write(phi);
}
