#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>
#include <vector>

//#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
//#include <io/nbfImageWriter.h>
//#include <vtkMath.h>
//#include <nbfVeselnessFilter.h>
//#include <nbfGaussianFilter.h>
//#include <nbfPolarDomain.h>
//#include <nbfMinimalSurface.h>
//#include <nbfGeodesicPath3D.h>
#include <fm/nbfFastFastMarching2D.h>

#define PIXEL float

void main( int argv, char ** argc )
{
	// 1. V - input image
	Array< PIXEL, 2 > V;

	nbfMatlabReader reader;
	reader.setFileName( argc[1] );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	V = V - min(V);
	V = V / max(V) + 1e-6;
	
	vector< TinyVector< int, 2 > > min1;
	vector< PIXEL > ini;
	min1.push_back( TinyVector< int, 2 >( 100, 100 ) );
	ini.push_back(0);
	
	nbfFastFastMarching2D< PIXEL > fastMarching( V );
	fastMarching.setAliveSet( min1, ini );

	Array< PIXEL, 2 > distance( V.shape() );

	nbfTimer t;
	t.start();
	fastMarching.execute( distance );
	t.stop();
	cout << t.elapsedSeconds() << endl;

	//cout << V.shape() << " = [ " << min(V) << "," << max(V) << "]" << endl;

	//Array< PIXEL, 2 > P;
	//Array< bool, 2 > B;
	//nbfPolarDomain< PIXEL, 2 > polar;

	//float rho = 260;
	//polar.setCenter( TinyVector< PIXEL, 2 >( 181, 237 ) );
	//polar.setMaxRho( rho );
	//polar.setResRho( rho );
	//polar.setResTheta( 360 );
	//polar.cartesian2polar( V, P, B );
	//P( Range(fromStart,25), Range::all() ) = numeric_limits< float >::max();

	nbfMatlabWriter wm;
	wm.setFileName("ipath");
	wm.write(distance);

	//nbfGeodesicPath< PIXEL > geodesic( P );
	//vector< TinyVector< PIXEL, 2 > > path;
	//Timer t;
	//t.start();
	//geodesic.getCircularPath( P, path );
	//t.stop();
	//cout << "t = " << t.elapsedSeconds() << endl;
	//geodesic.getImplicitPath(path,P);
	//polar.polar2cartesian( P, V, B );

	//wm.setFileName("fm8.array");
	//wm.write(V);
}