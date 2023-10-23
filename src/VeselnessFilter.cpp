#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <nbfDifferentials.h>
#include <nbfMatlabReader.h>
#include <nbfMatlabWriter.h>
#include <nbfFastMarching.h>
#include <nbfGeodesicPath.h>
#include <nbfLinearInterpolator.h>
#include <nbfPolarDomain.h>
#include <vtkMath.h>
#include <nbfVeselnessFilter.h>
#include <nbfGaussianFilter.h>

#define PIXEL float
#define DIM 3

void main( int argv, char ** argc )
{
	if ( argv != 6 ){
		cout << "Usage: FILE input.array output.array" << endl;
		cout << "\t input.array : input file name (from matlab)" << endl;
		cout << "\t output.array : output file name (to be read in matlab)" << endl;
		exit(0);
	}

	nbfMatlabReader reader;
	nbfMatlabWriter writer;

	Timer t;

	Array< PIXEL, 3 > V;
	reader.setFileName( argc[1] );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}
	cout << "V = " << V.shape() << endl;

	BordStrategyMirrorSimple< PIXEL,3 > bsForV( V, 1 );
	bsForV.refresh();

	Array< PIXEL, 3 > W( V.shape() );

	t.start();
#if 1
	for ( int i = 1; i < 30; i++ ){
		W = Laplacian3D(V);
		V = V + .1 * W;
		bsForV.refresh();
	}
#else
	nbfGaussianFilter< PIXEL, 3 > gauss( V );
	for ( int i = 0; i < 5; i++ ){
        gauss.execute(W);
		V = W;
		bsForV.refresh();
	}
#endif
	t.stop();
	cout << "gauss blur "  << t.elapsedSeconds() << " seconds." << endl;

	Array< PIXEL, 3 > Vx( V.shape() ), Vy( V.shape() ), Vz( V.shape() );

	nbfVeselnessFilter< PIXEL, 3 > vesel( V );
	t.start();
	vesel.execute(W,Vx,Vy,Vz);
	t.stop();
	cout << "vesel measure "  << t.elapsedSeconds() << " seconds." << endl;

	//// binarize
	//PIXEL th = nbfArrayFilter<PIXEL,3>::graythresh(W);

	//W = where( W > th, 0, 1 );

	cout << "W = " << W.shape() << endl;
	cout << "W = [" << min(W) << "," << max(W) << "]"<< endl;

	writer.setFileName( argc[2] );
	writer.write(W);

	Array< PIXEL, 3 > D( W.shape() ), E( W.shape() );
	Array< PIXEL, 3 > Ex( W.shape() ), Ey( W.shape() ), Ez( W.shape() );

	BordStrategyMirrorDouble< PIXEL,3 > bsForEx( Vx, 1 );
	bsForEx.refresh();
	BordStrategyMirrorDouble< PIXEL,3 > bsForEy( Vy, 1 );
	bsForEy.refresh();
	BordStrategyMirrorDouble< PIXEL,3 > bsForEz( Vz, 1 );
	bsForEz.refresh();

	D = where( W < 1e-2, 1 , 0 );

	BordStrategyMirrorDouble< PIXEL,3 > bsForW( W, 1 );
	bsForW.refresh();
#if 1
	for ( int i = 1; i < 20 + 1; i++ ){
		Ex = Laplacian3D(Vx);
		Ey = Laplacian3D(Vy);
		Ez = Laplacian3D(Vz);
		Vx = Vx + .1 * Ex * D;
		Vy = Vy + .1 * Ey * D;
		Vz = Vz + .1 * Ez * D;
		if ( (i>0) & (i % 5 == 0) ){
			E = sqrt( Vx * Vx + Vy * Vy + Vz * Vz );
			Vx = where( E > 0, Vx / E, Vx );
			Vy = where( E > 0, Vy / E, Vy );
			Vz = where( E > 0, Vz / E, Vz );
			cout << "normalize" << endl;
		}
		bsForEx.refresh();
		bsForEy.refresh();
		bsForEz.refresh();
		cout << i << endl;
	}
#endif

	// tensor diffusion
	t.start();
	
	for ( int i = 0; i < 0; i++ ){
		//Ex = Laplacian3D(W);
		Ex = central22n(W,firstDim) * ( 1 - Vx ) +
			central22n(W,secondDim) * ( 1 - Vy ) +
			central22n(W,thirdDim) * ( 1 - Vz ) -
			2 * mixed22n(W,firstDim,secondDim) * Vx * Vy -
			2 * mixed22n(W,firstDim,thirdDim) * Vx * Vz -
			2 * mixed22n(W,secondDim,thirdDim) * Vy * Vz;
		cout << "Update = [" << min(Ex) << "," << max(Ex) << "]"<< endl;
		W = W + .1 * Ex;
		bsForW.refresh();
	}
	t.stop();
	cout << "diffusion "  << t.elapsedSeconds() << " seconds." << endl;

	// binarize
	cout << "W = [" << min(W) << "," << max(W) << "]"<< endl;
	//W = W - min(W);
	//W = W / max(W);

	//th = nbfArrayFilter<PIXEL,3>::graythresh(W);

	//W = where( W > th, 0, 1 );

	writer.setFileName( argc[3] );
	writer.write(W);

	writer.setFileName( argc[4] );
	writer.write(Vx);

	writer.setFileName( argc[5] );
	writer.write(Vy);

}
