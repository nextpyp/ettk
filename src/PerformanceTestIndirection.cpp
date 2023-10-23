// Performance Test using modalities of Indirection

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>
#include <blitz/array/domain.h>

using namespace blitz;

#include <vector>

#include <Timer.hh>

#include <NarrowBand/nbfBordStrategyConst.h>
#include <NarrowBand/nbfDifferentials.h>

#define PIXEL float
#define DIM 2

BZ_DECLARE_STENCIL_OPERATOR1(nablaForward2D,A)
return sqrt( pow2( backward11n(A,firstDim) ) + 
	         pow2( forward11n(A,firstDim) ) );
BZ_END_STENCIL_OPERATOR

BZ_ET_STENCIL(nablaForward2D,P_numtype)

void main( int argv, char ** argc ){

    int N = 4;
	Array< PIXEL,2 > B(N,N), C(N,N);
	B = 1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12,
		13, 14, 15, 16;

	BordStrategyConst< PIXEL, 2 > BSforB( B, 1, 0 );
	BSforB.refresh();

	// build point set (one point)
	typedef TinyVector<int,2> coord;
	std::vector<coord> I;
	I.push_back( coord( 1,1 ) );
	I.push_back( coord( 1,2 ) );
	I.push_back( coord( 1,3 ) );

	// build point set (one point)
	std::vector< RectDomain<2> > J;
	J.push_back( strip( coord( 1,1 ), secondDim, 3 ) );

	C = 0;
	
	int times = 1e6;

	// indirection using STENCIL macro (FASTEST)
	Timer time;
	time.start();
	for ( int i = 0; i < times; i++ ){
		C[I] = nablaForward2D( B );
	}
	time.stop();
	cout << "Tiempo total del algoritmo (secs)= " << time.elapsedSeconds() << endl;
	cout << C << endl;

	// indirection using explicit expression (INTERMEDIATE)
	time.start();
	for ( int i = 0; i < times; i++ ){
		C[I] = sqrt( pow2( backward11n(B,firstDim) ) +
			         pow2( forward11n(B,firstDim) ) );	
	}
	time.stop();
	cout << "Tiempo total del algoritmo (secs)= " << time.elapsedSeconds() << endl;
	cout << C << endl;

	// full matrix operation (SLOWEST)
	time.start();
	for ( int i = 0; i < times; i++ ){
		C = sqrt( pow2( backward11n(B,firstDim) ) +
			      pow2( forward11n(B,firstDim) ) );
	}
	time.stop();
	cout << "Tiempo total del algoritmo (secs)= " << time.elapsedSeconds() << endl;
	cout << C << endl;
}
