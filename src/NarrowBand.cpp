#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <fstream>
#include <vector>

//#include <IO/nbfMrcReader.h>
#include <Timer.hh>
//#include <IO/FlujosIO.hh>

#include <NarrowBand/BordStrategyConst.h>
#include <NarrowBand/BordStrategyMirror.h>
#include <NarrowBand/Differentials.h>

#define PIXEL float
#define DIM 2

void main( int argv, char ** argc ){

	//Array< short, 3 > A;
	//nbfMrcReader reader;
	//reader.setFileName( argc[1] );
	//reader.setSubSample( 2, 2, 2 );
	//if ( reader.read( A ) ){
	//	cout << "Error reading file.\n";
	//}

	//SaveToFile(A,"full.blitz");

	//Array< short, 2 > B = A( Range::all(), 0, Range::all() );

	//B = ( B - min(B) );
	//B = B / max(B) * 255;
	//cout << min(B) << ", " << max(B) << endl;

	//WriteArray( B, "aver.bmp" );

	//ofstream writeField( "aver.blitz" );
	//if ( ~writeField.bad() ){
	//	writeField << data << endl;
	//	writeField.close();
	//}
	
    int N = 4;
//	Array< PIXEL, 3 > A( N, N, N );
	Array< PIXEL,2 > B(N,N), A(N,N), C(N,N);
	B = 1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12,
		13, 14, 15, 16;

	BordStrategyConst< PIXEL, 2 > BSforB( B, 1, 0 );
	BSforB.refresh();


	A = nablaForward2D( B );
	C = sqrt( pow2( where( backward11n(B,firstDim) > 0, backward11n(B,firstDim), 0) ) +
		      pow2( where( forward11n(B,firstDim) < 0, forward11n(B,firstDim), 0) ) );
	cout << A << endl;
	cout << C << endl;

	typedef TinyVector<int,2> coord;
	std::vector<coord> I;
	I.push_back( coord( 1,1 ) );
	//I.push_back( coord( 1,2 ) );
	//I.push_back( coord( 2,3 ) );
	//I.push_back( coord( 3,2 ) );
	C = 0;

	Timer time;
	time.start();
	for ( int i = 0; i < 10000; i++ ){
		//C[I] = sqrt( pow2( where( B > 0, B, 0) ) +
		//	         pow2( where( B < 0, B, 0) ) );
//		C[I] = sqrt( pow2( where( backward11n(B,firstDim) > 0, backward11n(B,firstDim), 0) ) +
//			         pow2( where( forward11n(B,firstDim) < 0, forward11n(B,firstDim), 0) ) );
		C[I] = nablaForward2D( B );
	}
	time.stop();
	cout << "Tiempo total del algoritmo (secs)= " << time.elapsedSeconds() << endl;
	cout << C << endl;

	time.start();
	for ( int i = 0; i < 10000; i++ ){
		//C = sqrt( pow2( where( B > 0, B, 0) ) +
		//	         pow2( where( B < 0, B, 0) ) );
		C[I] = sqrt( pow2( where( backward11n(B,firstDim) > 0, backward11n(B,firstDim), 0) ) +
			         pow2( where( forward11n(B,firstDim) < 0, forward11n(B,firstDim), 0) ) );
	}
	time.stop();
	cout << "Tiempo total del algoritmo (secs)= " << time.elapsedSeconds() << endl;
	cout << C << endl;

	time.start();
	for ( int i = 0; i < 10000; i++ ){
		//C = sqrt( pow2( where( B > 0, B, 0) ) +
		//	         pow2( where( B < 0, B, 0) ) );
		C = sqrt( pow2( where( backward11n(B,firstDim) > 0, backward11n(B,firstDim), 0) ) +
			         pow2( where( forward11n(B,firstDim) < 0, forward11n(B,firstDim), 0) ) );
	}
	time.stop();
	cout << "Tiempo total del algoritmo (secs)= " << time.elapsedSeconds() << endl;
	cout << C << endl;
	//for ( int i = 0; i < N; i++ ){
	//	A( Range::all(), Range::all(), i ) = B;
	//}

	//  for ( int i = 0; i < A.depth(); i++ ){
	//cout << A( Range::all(), Range::all(), i ) << endl;
    // }


	//// A = blitz::central12( A, firstDim );

	//// BSforA.refresh();

}
