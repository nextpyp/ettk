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
#include <nbfPolarDomain3.h>
#include <nbfMinimalSurface.h>

void main( int argv, char ** argc )
{
	if ( argv != 9 ){
		cout << "Parameters:" << endl;
		cout << "\t 1 : input file" << endl;
		cout << "\t 2 : centerX" << endl;
		cout << "\t 3 : centerY" << endl;
		cout << "\t 4 : centerZ" << endl;
		cout << "\t 5 : rho" << endl;
		cout << "\t 6 : laplacian iterations (=20)" << endl;
		cout << "\t 7 : segmentation result" << endl;
		cout << "\t 8 : output file (.vtk)" << endl;
		exit(0);
	}

	nbfMatlabReader reader;

	// multi purpose writer
	nbfMatlabWriter writer;
	writer.setFileName( "ipath" );

	Timer t;

	// 1. V - input image
	Array< float, 3 > V;

	reader.setFileName( argc[1] );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	// set to [0-255]
	V = V - min(V);
	V = V / max(V) * 255;

	// estimate gaussian parameters
	float v = sum( pow2( V - mean(V) ) ) / V.size() / 2.0;
	float lower = mean(V) - v;
	float upper = mean(V) + v;

	// restric to useful range
    V = where( V > lower, V, lower );
	V = where( V < upper, V, upper );
	
	// set to [0-1]
	V = V - min(V);
	V = V / max(V);

	V = V + 1e0;

	//nbfImageWriter vwriter;
	//vwriter.setFileName("V.vtk");
	//vwriter.write(V);

	//writer.write(V);

	//BordStrategyMirrorSimple< float,3 > bsForV( V, 1 );
	//bsForV.refresh();

	//Array< float, 3 > W( V.shape() );
	//for ( int i = 0; i < 5; i++ ){
	//	//W = Laplacian3D(V);
	//	W = curvature3D(V);
	//	V = V + .1 * W;
	//	bsForV.refresh();
	//}
	//nbfVeselnessFilter< float, 3 > vesel( V );
	//vesel.execute(W,thirdDim);
	//V = 1 - W + 1e0;

	//writer.write(V);
	//cout << min(W) << "," << max(W) << endl;

	// retrieve center
	float centerX = atof(argc[2]);
	float centerY = atof(argc[3]);
	float centerZ = atof(argc[4]);

	TinyVector< float, 3 > center( centerX, centerY, centerZ );

	float rho = atof(argc[5]);
	//float rho = 80;

	Array< float, 3 > P;
	Array< bool, 3 > B;

	nbfPolarDomain3< float > polar3;
	polar3.setCenter( TinyVector<float,3>( V.ubound(firstDim) / 2, V.ubound(secondDim) / 2, centerZ ) );	
	polar3.setMaxRho( rho );
	polar3.setResRho( 60 );
	polar3.setResTheta( 180 );
	polar3.setZScale(1.0);
	polar3.cartesian2polar(V,P,B);
	//polar3.polar2cartesian(P,V,B);
	//writer.write(P);
	//writer.write(V);
	P = where( B == true, P, mean(V) );
	//P = where( B == true, P, numeric_limits<float>::max() );
	P( Range(0,5), Range::all(), Range::all() ) = numeric_limits<float>::max();
	
	//P = where( P < numeric_limits<float>::max(), P, min(P) );
	//writer.write(P);

	nbfMinimalSurface< float > ms(P);

	Array< float, 3 > I;
	t.start();
	ms.search(P,I);
	t.stop();
	cout << "gauss blur "  << t.elapsedSeconds() << " seconds." << endl;

	P.resize( V.shape() );

	// get the raw image back
	P = V - 1.0;
	BordStrategyMirrorDouble< float, 3 > bsForP( P, 1 );
	bsForP.refresh();

	// gaussian smoothing to compute gradient
	for ( int i = 0; i < 5; i++ ){
		V = Laplacian3D(P);
		P = P + .1 * V;
		bsForP.refresh();
	}

	// compute image gradient
	Array< float, 3 > Vx( P.shape() ), Vy( P.shape() ), Vz( P.shape() );
	Vx = - central12n(P,firstDim);
	Vy = - central12n(P,secondDim);
	Vz = - central12n(P,thirdDim);

	polar3.polar2cartesian(I,V,B);	

	BordStrategyMirrorDouble< float, 3 > bsForV( V, 1 );
	bsForV.refresh();

	//Array< float, 3 > S( V.shape() );

	// parameters
	float timeStep = .25;	// evolution time step
	int iterations = 2;		// evolution iterations
	for ( int i = 0; i < iterations; i++ ){
		P = 1.0 * curvature3D(V) + 0.25 * advection3D(V,Vx,Vy,Vz);
		V = V + timeStep * P;
		bsForV.refresh();

	//	// reinitialization
	//	if ( i % 10 == 0 ){
			//int iterMorel = 2;

	//		// compute distances
			//S = sussmanSign( V );
			//do {
			//	P = morelSussman2D(V,S);
			//	V = V + 0.1 * P;
			//	bsForW.refresh();
			//	iterMorel--;
			//} while( iterMorel > 0 );
	//	}
	}

	//vwriter.setFileName("Q.vtk");
	//vwriter.write(V);

	writer.setFileName( argc[7] );
	writer.write(V);

	P = where( V < 0, 1, 0 );
	int inside = sum( P );

	// reload input image
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	cout << "Inside voxels = " << inside << " (of " << V.rows()*V.cols()*V.depth() << ")" << endl;

	Array< float, 1 > Ip( inside );
	Array< float, 1 > :: iterator iIp = Ip.begin();
	Array< float, 3 > :: iterator iterPhi = P.begin(),
		iterV = V.begin();
	while( iterPhi != P.end() ){
		if ( (*iterPhi) == 1 ){
            *iIp = *iterV;
			++iIp;
		}
		++iterPhi;
		++iterV;
	}
	cout << "Mean gray value = " << sum( Ip ) / inside << ", [" << min(V) << "," << max(V) << "]" << endl;

	writer.setFileName( argc[8] );
	writer.write(Ip);

}