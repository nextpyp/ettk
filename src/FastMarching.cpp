#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

//#include <MacrosFlujos.hh>
//#include <IO/FlujosIO.hh>
#include <nbfTimer.h>

#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <fm/nbfFastMarching3Dccc.h>

#define PIXEL double
#define DIM 3

void main( int argv, char ** argc )
{
//if ( argv != 4 ){
//    cout << "Usage: Distancer.exe inputFile.array outputFile.array initialPoint" << endl;
//    cout << "\t inputFile.array : input file name (from matlab)" << endl;
//    cout << "\t outputFile.array : output file name (to be read in matlab)" << endl;
//    cout << "\t initial point : select one of the four corners 1-4" << endl;
//    exit(0);
//  }

  Array< PIXEL, DIM > input;

  input.resize( 30, 30, 30 );
  firstIndex i;
  secondIndex j;
  thirdIndex k;
  input = i * j * k;
//  input( Range(24,26), Range::all(), Range::all() ) = .5;
	
//nbfMatlabReader reader;
  //reader.setFileName( argc[1] );
  //if ( reader.read( input ) ){
	 // cout << "Error reading file.\n";
	 // exit(0);
  //}
  
  nbfFastMarching3Dccc< PIXEL > fastMarching( input );

#if 0
  vector< TinyVector< int, DIM > > min1;
  vector< PIXEL > ini;

  for ( int i = 0; i < input.depth(); i++ ){
	  //for ( int j = 0; j < input.depth(); j++ ){
		  min1.push_back( TinyVector< int, DIM >(0,0,i) );
		  ini.push_back(0);
	  //}
  }			
#else
  TinyVector< int, DIM > min1(24,24,24);
  PIXEL ini = 0;
#endif

  Array< PIXEL, 3 > T(10,10,10);
  T = sqrt( 1.0 * pow2(i-5) + pow2(j-5) + pow2(k-5) );

  fastMarching.setAliveSet( min1, ini );
  fastMarching.setTemplate( T );

  Array< PIXEL, DIM > distance( input.shape() );

  nbfTimer tEvolution, tReinit;

  cout << "Running fast marching..." << endl;
  tEvolution.start();
  fastMarching.execute( distance );
  tEvolution.stop();
  cout << "Done." << endl;

  cout << "Tiempo total del algoritmo (secs)= " << tEvolution.elapsedSeconds() << endl;

  //TinyVector< int, 2 > fin( input.ubound(firstDim),input.ubound(secondDim));
  //vector< TinyVector< PIXEL, DIM > > points;
  //fastMarching.getPath( fin, points, input );

  //cout << points.size() << endl;

  nbfMatlabWriter writer;
  writer.setFileName( "distance3d.array" );
  writer.write( distance );
}
