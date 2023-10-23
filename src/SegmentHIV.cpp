#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <MacrosFlujos.hh>
#include <IO/FlujosIO.hh>
#include <Timer.hh>

#include <NarrowBand/nbfMatlabReader.h>
#include <NarrowBand/nbfMatlabWriter.h>
#include <NarrowBand/nbfFastMarching.h>
#include <NarrowBand/nbfGeodesicPath.h>

#define PIXEL float
#define DIM 2

void main( int argv, char ** argc )
{
if ( argv != 4 ){
    cout << "Usage: Distancer.exe inputFile.array outputFile.array initialPoint" << endl;
    cout << "\t inputFile.array : input file name (from matlab)" << endl;
    cout << "\t outputFile.array : output file name (to be read in matlab)" << endl;
    cout << "\t initial point : select one of the four corners 1-4" << endl;
    exit(0);
  }

  Array< PIXEL, DIM > input;

nbfMatlabReader reader;
  reader.setFileName( argc[1] );
  if ( reader.read( input ) ){
	  cout << "Error reading file.\n";
	  exit(0);
  }
  
  Array< PIXEL, DIM > imaPeso( input.shape() );
  imaPeso = input;
  
nbfFastMarching< PIXEL, DIM > fastMarching( imaPeso );
  vector< TinyVector< int, DIM > > positions;
  vector< PIXEL > distances;

  // get infinity value from input image
  float infty = max(input);

  // set to Alive all points in the first row of the image with distance < \infty.
  for ( int j = 0; j < input.cols(); j++ ){
	  if ( input( input.lbound(firstDim), j ) < infty ){
		  TinyVector< int, DIM > last( input.lbound(firstDim), j );
		  positions.push_back( last );
		  distances.push_back( 0 );
	  }
  }
  fastMarching.setAliveSet( positions, distances );

  // stop when we reach the other side
  fastMarching.setStopBorder( firstDim );
  fastMarching.setFreezing(10);

  Array< PIXEL, DIM > distance( input.shape() );

  Timer tEvolution;
  tEvolution.start();
  fastMarching.execute( distance );
  tEvolution.stop();

  cout << "Tiempo total del algoritmo (secs)= " << tEvolution.elapsedSeconds() << endl;

  // find end path point as the first to reach the other side.
  int endPoint = ( minIndex( distance( input.lbound(firstDim), Range::all() ) ) )(firstDim);  
  TinyVector< int, DIM > fin(input.lbound(firstDim), endPoint );

  //fastMarching.getPath( fin, input );

  tEvolution.start();
nbfGeodesicPath< PIXEL > geodesic( input );
  vector< TinyVector< int, 2 > > path;
  geodesic.execute( path, input );
  tEvolution.stop();
  cout << "Tiempo total del algoritmo (secs)= " << tEvolution.elapsedSeconds() << endl;

nbfMatlabWriter writer;
  writer.setFileName( argc[2] );
  writer.write(input);
}
