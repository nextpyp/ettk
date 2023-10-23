#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <nbf3DReader.h>
#include <nbfBordStrategyMirror.h>
#include <nbfMatlabReader.h>
#include <nbfMatlabWriter.h>
#include <nbfMorphologyFilters.h>

void main( int argv, char ** argc )
{
	if ( argv != 7 ){
		cout << "Parameters:" << endl;
		cout << "\t 1 : input file" << endl;
		cout << "\t 2 : gradient threshold (=2000)" << endl;
		cout << "\t 3 : erosion steps to remove small holes (=1)" << endl;
		cout << "\t 4 : dilation steps to extend holes (=2)" << endl;
		cout << "\t 5 : factor of \sigma to estimate gaussian width (=6)" << endl;
		cout << "\t 6 : output file" << endl;
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

	// 2. remove markers

	t.start();

	Array< bool, 3 > Vth( V.shape() );

	BordStrategyMirrorSimple< bool,3 > bsForVth( Vth, 1 );
	bsForVth.refresh();

	Vth = where( sqrt( pow2( central12n(V,firstDim) ) +
		pow2( central12n(V,secondDim)) +
		pow2( central12n(V,thirdDim) ) ) < atoi(argc[2]), false, true );

	Array< bool, 3 > R( V.shape() );

	// erode small holes : default = 1
	for ( int i = 0; i < atoi( argc[3] ); i++ ){
        R = erode3D(Vth);
		Vth = R;
		bsForVth.refresh();
	}

	// dilate remaining holes : default = 2
	for ( int i = 0; i < atoi( argc[4] ); i++ ){
		R = dilate3D(Vth);
		Vth = R;
		bsForVth.refresh();
	}

	// compute statistics to estimate the gaussian parameters
	double meanV = mean( V );
	double sigma = sqrt( sum( ( V - meanV ) * ( V - meanV ) / V.size() ) );

	// take a fixed fraction of the gaussian aperture : default = 6
	double factor = atof( argc[5] );
	double thV = factor * sigma;

	cout << "mean = " << meanV << ", sigma = " << sigma << endl;
	cout  << "[" << meanV - thV << "," << meanV + thV << "]" << endl;

	// threshold below and above (hopefully getting rid of outliers)

	V = where( V < meanV - thV, meanV - thV, V );
	V = where( V > meanV + thV, meanV + thV, V );
	V = where( Vth == true, meanV, V );

	cout << "V = [" << min(V) << "," << max(V) << "]" << endl;

	//V = V - min(V);
	//V = V / max(V);
	//bsForV.refresh();

	t.stop();
	cout << "time elapsed "  << t.elapsedSeconds() << " seconds." << endl;

	writer.setFileName( argc[6] );
	writer.write(V);
}