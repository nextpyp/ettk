#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <cut/nbfCutFastMarching.h>

#define PIXEL float
#define DIM 3

void main( int argv, char ** argc )
{

  TinyVector< int, DIM > center(25,25,25);

  Array< PIXEL, DIM > input;

  input.resize( 50, 50, 50 );
  input = 1;
  //firstIndex i;
  //secondIndex j;
  //thirdIndex k;
  //input = pow2(i-center(firstDim)) + pow2(j-center(secondDim)) + pow2(k-center(thirdDim));
  //input = 1 / sqrt(input);

//nbfMatlabReader reader;
  //reader.setFileName( argc[1] );
  //if ( reader.read( input ) ){
	 // cout << "Error reading file.\n";
	 // exit(0);
  //}

  nbfCutFastMarching3D< PIXEL > fastMarching( input, center );

#if 1
  vector< TinyVector< int, 3 > > positions;
  vector< PIXEL > distances;

  for ( int k = 0; k < input.depth(); k++ ){
	  for ( int i = center[firstDim] + 1; i < input.rows(); i++ ){
			  TinyVector< int, 3 > last( i, center[secondDim], k );
			  positions.push_back( last );
			  distances.push_back( 0 );
			  //cout << "alive = [" << i << ", " << center[secondDim] << "," << k << "]" << endl;
	  }
  }

  fastMarching.setAliveSet( positions, distances );
#else
  TinyVector< int, DIM > min1(25,49,25);
  PIXEL ini = 0;

  fastMarching.setAliveSet( min1, ini );
#endif

  Array< PIXEL, DIM > distance( input.shape() );

  Timer tEvolution, tReinit;

  cout << "Running fast marching...\n";
  tEvolution.start();
  fastMarching.execute( distance );
  tEvolution.stop();
  cout << "Done.\n";

  cout << "Tiempo total del algoritmo (secs)= " << tEvolution.elapsedSeconds() << endl;

  //TinyVector< int, 2 > fin( input.ubound(firstDim),input.ubound(secondDim));
  //vector< TinyVector< PIXEL, DIM > > points;
  //fastMarching.getPath( fin, points, input );

  //cout << points.size() << endl;

  nbfMatlabWriter writer;
  writer.setFileName( argc[1] );
  writer.write( distance );
}