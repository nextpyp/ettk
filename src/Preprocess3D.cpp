#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <nbfBordStrategyMirror.h>
#include <nbfMatlabReader.h>
#include <nbfMatlabWriter.h>
#include <vtkMath.h>
#include <nbfVeselnessFilter.h>
#include <nbfGaussianFilter.h>

void main( int argv, char ** argc )
{
	if ( argv != 4 ){
		cout << "Parameters:" << endl;
		cout << "\t 1 : input file" << endl;
		cout << "\t 2 : laplacian iterations (=20)" << endl;
		cout << "\t 3 : output file" << endl;
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
	cout << "V = " << V.shape() << endl;
	cout << "V = [" << min(V) << "," << max(V) << "]" << endl;

	BordStrategyMirrorSimple< float,3 > bsForV( V, 1 );
	bsForV.refresh();

	t.start();

	Array< float, 3 > W( V.shape() );

	// gaussian smoothing
	for ( int i = 0; i < atoi( argc[2] ); i++ ){
		W = Laplacian3D(V);
		V = V + .1 * W;
		bsForV.refresh();
	}

	t.stop();
	cout << argc[2] << " steps of laplacian diffusion "  << t.elapsedSeconds() << " seconds." << endl;

	nbfVeselnessFilter< float, 3 > vesel( V );
	t.start();
	vesel.execute(W);
	t.stop();
	cout << "vesel measure "  << t.elapsedSeconds() << " seconds." << endl;

	// binarize
	float th = nbfArrayFilter<float,3>::graythresh(W);

	W = where( W > th, th, W ) / th;
	W = 1 - W;

	writer.setFileName( argc[3] );
	writer.write(W);
}