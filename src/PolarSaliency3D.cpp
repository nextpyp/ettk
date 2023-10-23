#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfImageWriter.h>
#include <vtkMath.h>
#include <nbfVeselnessFilter.h>
#include <nbfGaussianFilter.h>
#include <nbfPolarDomain.h>
#include <nbfMinimalSurface.h>

void main( int argv, char ** argc )
{
	if ( argv != 8 ){
		cout << "Parameters:" << endl;
		cout << "\t 1 : input file" << endl;
		cout << "\t 2 : centerX" << endl;
		cout << "\t 3 : centerY" << endl;
		cout << "\t 4 : centerZ" << endl;
		cout << "\t 5 : rho" << endl;
		cout << "\t 6 : laplacian iterations (=20)" << endl;
		cout << "\t 7 : output file (.vtk)" << endl;
		exit(0);
	}

	nbfMatlabReader reader;
	nbfMatlabWriter writer;

	Timer t;

	// 1. V - input image
	Array< float, 3 > V;

	reader.setFileName( argc[1] );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	// B2
	int centerX = atoi(argc[2]);
	int centerY = atoi(argc[3]);
	int centerZ = atoi(argc[4]);
	TinyVector< int, 3 > center( centerX, centerY, centerZ );

	BordStrategyMirrorSimple< float,3 > bsForV( V, 1 );
	bsForV.refresh();

	//t.start();

	Array< float, 3 > W( V.shape() );

	// gaussian smoothing
	for ( int i = 0; i < atoi( argc[6] ); i++ ){
		W = Laplacian3D(V);
		V = V + .1 * W;
		bsForV.refresh();
	}

	//t.stop();
	//cout << "\t" << argc[6] << " steps of laplacian diffusion = "  << t.elapsedSeconds() << " seconds." << endl;

	Array< float, 3 > Wx( W.shape() );
	Array< float, 3 > Wy( W.shape() );
	Array< float, 3 > Wz( W.shape() );

	nbfVeselnessFilter< float, 3 > vesel( V );

	t.start();
	vesel.execute(Wx,Wy,Wz);
	t.stop();
	cout << "\tSaliency measure = "  << t.elapsedSeconds() << " seconds." << endl;

	cout << "[" << min(Wx) << "," << max(Wx) << "]" << endl;
	cout << "[" << min(Wy) << "," << max(Wy) << "]" << endl;
	cout << "[" << min(Wz) << "," << max(Wz) << "]" << endl;

	Wx = 1 - sqrt( sqrt( sqrt( sqrt(Wx) ) ) ) + 1e-6;
	Wy = 1 - sqrt( sqrt( sqrt( sqrt(Wy) ) ) ) + 1e-6;
	Wz = 1 - sqrt( sqrt( sqrt( sqrt(Wz) ) ) ) + 1e-6;

	cout << "[" << min(Wx) << "," << max(Wx) << "]" << endl;
	cout << "[" << min(Wy) << "," << max(Wy) << "]" << endl;
	cout << "[" << min(Wz) << "," << max(Wz) << "]" << endl;

//	writer.setFileName( argc[7] );
	writer.setFileName( "ipath" );
	//writer.write(Wx);
	//writer.write(Wy);
	writer.write(Wz);

	//// shrink range to 0-1 and compress values
	//W = 1 - sqrt( sqrt( sqrt( sqrt(W) ) ) );
	//W = where( A == true, W, numeric_limits< float > :: max() );
	//// W = where( A == true, W, 1 );

	//cout << "\tW = " << W.shape() << endl;
	//cout << "\tW = [" << min(W) << "," << max(W) << "]" << endl;

	////polar.polar2cartesian(W,V);

	nbfMinimalSurface< float > ms(W);

	t.start();
	ms.search(Wx,Wy,Wz,center,W);
	t.stop();
	cout << "\tminimal surface = "  << t.elapsedSeconds() << " seconds." << endl;

	writer.write(W);
}