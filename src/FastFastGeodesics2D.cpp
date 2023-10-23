#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfImageReader.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfImageWriter.h>
#include <vtkMath.h>
#include <nbfGaussianFilter.h>
#include <nbfFastGeodesicPath.h>
#include <nbfGeodesicPath.h>
#include <fm/nbfFastFastMarching2D.h>
#include <fm/nbfFastMarching2D.h>

#include <random/normal.h>

#define PIXEL float

void main( int argv, char ** argc )
{
	Timer t;
	t.start();

#if 1
	Array< PIXEL, 2 > ge;//( 2000, 2000 );
	//ge = 1;

	//ranlib::NormalUnit<PIXEL> normalGen;
	//Array< PIXEL, 2 >::iterator iter = ge.begin();
	//while ( iter != ge.end() ){
	//	(*iter) = normalGen.random();
	//	iter++;
	//}
	//ge = ge - min(ge) + 1e-3;

	//PIXEL * myMatrix;
	//int size = 10;
	//myMatrix = new PIXEL[size*size];
	//for ( int i = 0; i < 10; i++ ){
	//	for ( int j = 0; j < 10; j++ ){
	//		myMatrix[ i * size + j ] = 1;
	//	}
	//}

	//cout << myMatrix[] << endl;

	//PIXEL (*aver)[2000] = new PIXEL [2000][2000];

	//Timer te;
	//te.start();
	//for ( int i = 0; i < 1e8; i++ ){
	//	ge(1000,1000) = ge(1000,1000) + ge(1001,1000) + ge(1000,1001) + ge(1001,1001);
	//}
	//te.stop();
	//cout << te.elapsedSeconds() << endl;

	//te.start();
	//for ( int i = 0; i < 1e8; i++ ){
	//	*aver[1000,1000] = *aver[1000,1000] + *aver[1001,1000] + *aver[1000,1001] + *aver[1001,1001];
	//}
	//te.stop();
	//cout << te.elapsedSeconds() << endl;

	//delete [] aver;

	nbfImageReader reader;
	reader.setFileName( argc[1] );
	if ( reader.read( ge ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	//ge = ge - min(ge) + 1e-3;
	//ge = ge / max(ge);
	ge = where( ge < .5, 1e-1, 2 );
	////ge = ge / max(ge);
	

	nbfMatlabWriter wm;
	wm.setFileName("distance.array");

#else
	// 1. V - input image
	Array< PIXEL, 3 > V;

	nbfMatlabReader reader;
	reader.setFileName( argc[1] );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	V = V - min(V);
	V = V / max(V);
	V = 1 - V + 1e-2;

	cout << V.shape() << " = [ " << min(V) << "," << max(V) << "]" << endl;

	Array< PIXEL, 2 > BB( V( Range::all(), Range::all(), 0 ) );
	BB = BB - min(BB) + 1e-3;
	BB = BB / max(BB);

	wm.write( BB );
#endif
	
	// MAZE EXAMPLE - starting and end points
	TinyVector<int,2> start(372,560);
	TinyVector<int,2> end(362,250);
	//TinyVector<int,2> end2(387,862);

	// LAKE EXAMPLE - starting and end points
	// ge = where( ge < .5, 1e-1, 2 );
	//TinyVector<int,2> end(554,1021);
	//TinyVector<int,2> start(592,350);
	//TinyVector<int,2> start(501,742);
	//TinyVector<int,2> start(694,607);

	// UNIFORM EXAMPLE
	//TinyVector<int,2> start(0,0);
	//TinyVector<int,2> end = ge.ubound();

	//nbfFastFastMarching2D< PIXEL > fm2( ge, 1000 );
	nbfFastMarching2D< PIXEL > fm2( ge );
	vector< TinyVector< int, 2 > > aP;
	aP.push_back( start );
	vector< PIXEL > aD;
	aD.push_back(0.0);
	fm2.setAliveSet(aP,aD);
	//fm2.setStopPoint( end );

	Array< PIXEL, 2 > dist( ge.shape() );

	t.stop();
	cout << "Initialization = " << t.elapsedSeconds() << endl;

	t.start();
	fm2.execute(dist);
	t.stop();
	cout << "Distance computation = " << t.elapsedSeconds() << ", size = " << ge.rows() << "x" << ge.cols() << endl;

	//return;

	wm.write(ge);
	wm.write(dist);

	return;

	nbfGeodesicPath< PIXEL > g(dist);

	vector< TinyVector< PIXEL, 2 > > path;
	t.start();
	g.getSimplePath( dist, end, path );
	t.stop();
	cout << "t = " << t.elapsedSeconds() << endl;
	cout << path.size() << endl;
	//g.getImplicitPath(path,ge);

	dist = 1;
	Array< PIXEL, 2 > line( path.size(), 2 );
	for ( int i = 0; i < path.size(); i++ ){
		line(i,0) = path[i](firstDim);
		line(i,1) = path[i](secondDim);
		dist( (int)line(i,0), (int)line(i,1) ) = 0;
		dist( (int)line(i,0)-1, (int)line(i,1) ) = 0;
		dist( (int)line(i,0), (int)line(i,1)-1 ) = 0;
		dist( (int)line(i,0)-1, (int)line(i,1)-1 ) = 0;
		dist( (int)line(i,0)+1, (int)line(i,1) ) = 0;
		dist( (int)line(i,0), (int)line(i,1)+1 ) = 0;
		dist( (int)line(i,0)+1, (int)line(i,1)+1 ) = 0;
	}

	wm.write(line);
	wm.write(dist);

	return;
}