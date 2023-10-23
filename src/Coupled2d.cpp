#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>

#include <vtkMath.h>
#include <nbfDifferentials.h>
#include <bs/nbfBordStrategyMirror.h>

#define PIXEL float
#define DIM 2

void main( int argv, char ** argc ){

	nbfMatlabReader reader;

	Array< PIXEL, DIM > I1;
	reader.setFileName( "I1" );
	reader.read( I1 );
	cout << I1.shape() << endl;

	Array< PIXEL, DIM > phi1;
	reader.setFileName( "phi1" );
	reader.read( phi1 );
	cout << phi1.shape() << endl;

	Array< PIXEL, DIM > I2;
	reader.setFileName( "I2" );
	reader.read( I2 );
	cout << I2.shape() << endl;

	Array< PIXEL, DIM > phi2;
	reader.setFileName( "phi2" );
	reader.read( phi2 );
	cout << phi2.shape() << endl;

	BordStrategyMirrorDouble< PIXEL, DIM > BSfor1( phi1, 1 );
	BordStrategyMirrorDouble< PIXEL, DIM > BSfor2( phi2, 1 );
	BSfor1.refresh();
	BSfor2.refresh();

	// tmp arrays for computations
	Array< PIXEL, DIM > R1( I1.shape() ), R2( I2.shape() );
	Array< PIXEL, DIM > S1( I1.shape() ), S2( I2.shape() );

	// parameters
	float timeStep = .25;	// evolution time step
	int iterations = 150;	// evolution iterations
	int gapMorel = 5;		// re-initialize every number of iterations
	float stepMorel = .2;	// morel PDE time step

	TinyVector< int, 2 > center(81,19);

	nbfMatlabWriter writer;
	writer.setFileName("ipath");

	for ( int i = 0; i < iterations; i++ ){
	
		R1 = where(I1<.99,0,1) * curvature2D(phi1);
		R2 = where(I2<.99,0,1) * curvature2D(phi2);

		phi1 = phi1 + timeStep * R1;
		BSfor1.refresh();

		phi2 = phi2 + timeStep * R2;
		BSfor2.refresh();

		writer.write(phi1);
		writer.write(phi2);

		// force boundary condition
		Range I( center[firstDim] - 3, center[firstDim] + 3 );
		phi1( I, Range::all() ) = ( phi1( I, Range::all() ) + phi2( I, Range::all() ) ) / 2.0;
		phi2( I, Range::all() ) = phi1( I, Range::all() );

		writer.write(phi1);
		writer.write(phi2);

		// reinitialization
		if ( i % gapMorel == 0 ){

			// compute distances
			S1 = sussmanSign( phi1 );
			double stop = 1;
			do {
				R1 = morelSussman2D(phi1,S1);
				phi1 = phi1 + stepMorel * R1;
				BSfor1.refresh();
				stop = abs( max(R1) );
			} while( stop > 1 );

			S2 = sussmanSign( phi2 );
			do {
				R2 = morelSussman2D(phi2,S2);
				phi2 = phi2 + stepMorel * R2;
				BSfor2.refresh();
				stop = abs( max(R2) );
			} while( stop > 1 );
		}

		writer.write(phi1);
		writer.write(phi2);

	}

	writer.setFileName( "phi1f" );
	writer.write(phi1);	
	writer.setFileName( "phi2f" );
	writer.write(phi2);	
}
