#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <io/nbfMatlabReader.h>
#include <io/nbfImageReader.h>
#include <io/nbfMatlabWriter.h>
#include <fm/nbfFastFastMarching2D.h>
//#include <cut/nbfCutGeodesics.h>
#include <random/normal.h>
#include <vtkMath.h>


#define PIXEL float

void main( int argv, char ** argc )
{
	Array< PIXEL, 2 > input;

	TinyVector< float, 2 > center( 64,82 );
	//TinyVector< int, 2 > center( 80,80 );

	nbfMatlabWriter writer;
	writer.setFileName("d.array");
#if 0
//  nbfMatlabReader reader;
  nbfImageReader reader;
  reader.setFileName( argc[1] );
  if ( reader.read( input ) ){
	  cout << "Error reading file.\n";
	  exit(0);
  }
  cout << input.shape() << endl;
  input = input - min(input);
  input = input / max(input);
  //input = 1;
  input = input + 1e-3;
#else
  int n = atoi( argc[1] );
  input.resize(n,n);
  //input = 1;

  center = input.shape() / 2;

  ranlib::Normal<float> normalGen(.5,.1);
  Array<float,2>::iterator iter = input.begin();
  while ( iter != input.end() ){
	  (*iter) = floor( normalGen.random() * 255 );   
	  iter++;
  }
  input = input - min(input);
  input = input / max(input) / n;
  input = input + 1e-10;

#endif

#if 0
  vector< TinyVector< int, 2 > > min1;
  vector< PIXEL > ini;
  for ( int i = center[firstDim] + 1; i < input.rows(); i++ ){
	  min1.push_back( TinyVector< int, 2 >(i,center(secondDim)) );
	  ini.push_back(0);
  }			
#else
  TinyVector< int, 2 > min1 = input.shape() / 2;
  //min1(firstDim) = 0;
  //min1(secondDim) = 0;
  //min1(secondDim) = min1(secondDim) + 1;
  //cout << "init = " << min1 << endl;
  PIXEL ini = 0;
#endif

  //nbfMatlabWriter writer;
  //writer.setFileName("d.array");
  //writer.write(input);

  //nbfFastFastMarching2D< PIXEL > fastMarching( input, center );
  nbfFastFastMarching2D< PIXEL > fastMarching( input, atoi( argc[2] ) );
  fastMarching.setAliveSet( min1, ini );

  Array< PIXEL, 2 > distance( input.shape() );
  
  //Timer t;
  //t.start();
  fastMarching.execute( distance );
  //t.stop();
  //cout << t.elapsedSeconds() << endl;
  //writer.write(distance);

  cout << max(distance) << endl;
  return;

  //nbfCutGeodesics< PIXEL > geodesic( input, center );
  //Array< PIXEL, 2 > ipath( input.shape() );
  //TinyVector< int, 2 > start( 86, 82 );
  //vector< TinyVector< int, 2 > > path;
  ////geodesic.getForwardPath( start, path );
  //geodesic.getNewCircularPath( input, start, ipath );

  //ipath = -1;
  //for ( int i = 0; i < path.size(); i++ )
  //{
	 // ipath(path[i]) = 1;
  //}

  //writer.setFileName("i.array");
  //writer.write(ipath);

}