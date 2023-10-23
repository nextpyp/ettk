#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <fstream>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <bs/nbfBordStrategyConst.h>

#include <nbfGaussianFilter.h>
#include <nbfTimer.h>

#define PIXEL float
#define DIM 3

BZ_DECLARE_STENCIL_OPERATOR1(miLaplacian2D,A)															
return -4.0 * A(0,0)
    + A(-1,0) + A(1,0)
    + A(0,-1) + A(0,1);
BZ_END_STENCIL_OPERATOR																	

BZ_ET_STENCIL(miLaplacian2D, P_numtype)

void main( int argv, char ** argc )
{
	//if ( argv != 3 ){
	//	cout << "Usage: FILE input.array output.array" << endl;
	//	cout << "\t input.array : input file name (from matlab)" << endl;
	//	cout << "\t output.array : output file name (to be read in matlab)" << endl;
	//	exit(0);
	//}

	Array< PIXEL, 3 > A, B, C;
#if 0
	nbfMatlabReader reader;
	reader.setFileName( argc[1] );
	if ( reader.read( A ) ){
		cout << "Error reading file.\n";
		exit(0);
	}
#else
	A.resize(100,100,100);
	A = tensor::i + tensor::j + tensor::k;
#endif
	B.resize(A.shape());
	C.resize(A.shape());

	// bord strategy with appropiate band size
	BordStrategyConst<PIXEL,3> bsForA(A,5,0);
	bsForA.refresh();

	nbfGaussianFilter< PIXEL, 3 > filter( A );

	Timer t;
	t.start();
	filter.execute(B);
	t.stop();

	cout << "B = [" << min(B) << "," << max(B) << "]" << endl;
	cout << "Time = "  << t.elapsedSeconds() << " seconds." << endl;

	t.start();
	for ( int i = 0; i < 100; i++ ){
		B = Laplacian2D(A);
	}
	t.stop();
	cout << "Time = "  << t.elapsedSeconds() << " seconds." << endl;

	t.start();
	for ( int i = 0; i < 100; i++ ){
		C = miLaplacian2D(A);
	}
	t.stop();
	cout << "Time = "  << t.elapsedSeconds() << " seconds." << endl;

	cout << "diff = " << max( fabs(B-C) ) << endl;

	nbfMatlabWriter writer;
	writer.setFileName( argc[2] );
	writer.write(B);
}
