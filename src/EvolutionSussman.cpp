#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vector>

#include <Timer.hh>
#include <MacrosFlujos.hh>
#include <IO/FlujosIO.hh>
#include <NarrowBand/nbfMatlabReader.h>
#include <NarrowBand/nbfMatlabWriter.h>

#include <NarrowBand/nbfBordStrategyConst.h>
#include <NarrowBand/nbfBordStrategyMirror.h>
#include <NarrowBand/nbfDifferentials.h>

#include <Timer.hh>

#define PIXEL float
#define DIM 2

void main( int argv, char ** argc ){

	if ( argv != 5 ){
    cout << "Usage: sussman.exe input.array output.array stop step" << endl;
    cout << "input.array  : input file name (array format)" << endl;
    cout << "output.array : output file name (array format)" << endl;
    cout << "stop         : stop condition (1e-2)" << endl;
    cout << "step         : time step (.5)" << endl;
    exit(0);
	}

	Array< float, 2 > A;
	nbfMatlabReader reader;
	reader.setFileName( argc[1] );
	reader.read( A );

	BordStrategyMirrorDouble< PIXEL, 2 > BSforA( A, 2 );

	// parameters
	float stopMorel = atof( argc[3] );	// stop condition for morel PDE
	float stepMorel = atof( argc[4] );	// morel PDE time step

	Timer time;
	time.start();

	Array< float, 2 > S( A.shape() ), R( A.shape() );

	// compute distances

	S = sussmanSign( A );
	do {
		R = morelSussman2D(A,S);
		A = A + stepMorel * R;
		BSforA.refresh();
	} while( max(abs(R)) > stopMorel );
	
	time.stop();
	cout << "Tiempo total del algoritmo (secs)= " << time.elapsedSeconds() << endl;

	nbfMatlabWriter writer;
	writer.setFileName( argc[2] );
	writer.write(A);
}
