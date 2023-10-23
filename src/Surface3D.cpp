#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <nbfDifferentials.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfImageWriter.h>
#include <fm/nbfFastMarching.h>
#include <nbfGeodesicPath.h>
#include <nbfLinearInterpolator.h>
#include <nbfPolarDomain.h>
#include <vtkMath.h>
#include <nbfVeselnessFilter.h>
#include <nbfGaussianFilter.h>
#include <nbfMinimalSurface.h>
#include <bs/nbfBordStrategyMirror.h>

void main( int argv, char ** argc )
{
	if ( argv != 12 ){
		cout << "Parameters:" << endl;
		cout << "\t 1 : input file" << endl;
		cout << "\t 2 : center[x]" << endl;
		cout << "\t 3 : center[y]" << endl;
		cout << "\t 4 : center[z]" << endl;
		cout << "\t 5 : rho" << endl;
		cout << "\t 6 : P file (polar weights)" << endl;
		cout << "\t 7 : I file (polar implicit)" << endl;
		cout << "\t 8 : D file (polar distances)" << endl;
		cout << "\t 9 : Q file (implicit)" << endl;
		cout << "\t 10: V file" << endl;
		cout << "\t 11: output data file" << endl;
		exit(0);
	}

	Array< float, 3 > V, W, P, Q, D, I;

	nbfMatlabReader reader;
	nbfImageWriter writer;

	Timer t;
	double total = 0;

	reader.setFileName( argc[1] );
	if ( reader.read( W ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	writer.setFileName("tmp.vtk");
	writer.write(W);

	nbfPolarDomain< float, 3 > polar;

	// B2
	int centerX = atoi(argc[2]);
	int centerY = atoi(argc[3]);
	int centerZ = atoi(argc[4]);

	float rho = atof(argc[5]);

//	cout << "center = [" << centerX << "," << centerY << "," << centerZ << "]" << endl;
//	cout << "rho = " << rho << endl;

	polar.setCenter( TinyVector<float,3>( W.ubound(firstDim) / 2, W.ubound(secondDim) / 2, centerZ ) );

//	cout << "relative center = [" << W.ubound(firstDim) / 2 << "," << W.ubound(secondDim) / 2 << "," << centerZ << "]" << endl;

	polar.setMaxRho( rho );
	polar.setMinRho( 10.0 );
	polar.setResRho( 1.5 * rho );
	polar.setResTheta( 2 * vtkMath::Pi() / atanf( 1.0 / rho ) );
	polar.setResTheta( vtkMath::Pi() / atanf( 1.0 / rho ) );

	t.start();

#if 0
	reader.setFileName( argc[6] );
	if ( reader.read( P ) ){
		cout << "Error reading file.\n";
		exit(0);
	}
#else

	Array< bool, 3 > A;
	polar.cartesian2polar(W,P,A);
	
	//// set out of range point to fixed weight (!= \infty)
	//P = where( P < numeric_limits<float>::max(), P, 1 );

	// avoid small circles (mask weights close to rho = 0)
	int nullZoneX = P.rows() / 10;
	//int nullZoneX = P.rows() / 5;
	//P( Range( fromStart, nullZoneX - 1), Range::all(), Range::all() ) = 1;
	P( Range( fromStart, nullZoneX - 1), Range::all(), Range::all() ) = numeric_limits<float>::max();

	writer.setFileName( argc[6] );
	writer.write(P);

	// avoid angular distortion (mask weights close to: \theta = 0 and \theta = pi)
	int nullZoneZ = P.depth() / 8;
	Array< float, 3 > lslice( P( Range::all(), Range::all(), Range( fromStart, nullZoneZ - 1) ) );
	lslice = where( lslice < numeric_limits<float>::max(), 2, lslice );
	Array< float, 3 > uslice( P( Range::all(), Range::all(), Range( P.depth() - nullZoneZ, toEnd) ) );
	uslice = where( uslice < numeric_limits<float>::max(), 2, uslice );

	//polar.polar2cartesian(P,Q);	

	//writer.setFileName( argc[9] );
	//writer.write(Q);

	t.stop();
	cout << "\tto polar = "  << t.elapsedSeconds() << " seconds." << endl;
	total += t.elapsedSeconds();

	cout << "\tP = " << P.shape() << endl;

	P = where( P < numeric_limits<float>::max(), P + 1e-4, numeric_limits<float>::max() );

	writer.setFileName( argc[6] );
	writer.write(P);
#endif

	nbfMinimalSurface< float > ms(P);

	t.start();
	ms.executeNew(I,D);
	t.stop();
	cout << "\tminimal surface = "  << t.elapsedSeconds() << " seconds." << endl;
	//total += t.elapsedSeconds();

	writer.setFileName( argc[7] );
	writer.write(I);

	t.start();
	polar.polar2cartesian(I,Q,A);
	t.stop();
	cout << "\tback to cartesian = "  << t.elapsedSeconds() << " seconds." << endl;
	//total += t.elapsedSeconds();

	// fix missing interpolated values
	Q = where( Q < numeric_limits<float>::max(), Q, - numeric_limits<float>::max() );
	float mQ = max(Q);
	Q = where( Q > - numeric_limits<float>::max(), Q, mQ ); 

	writer.setFileName( argc[9] );
	writer.write(Q);

	// REGULARIZATION

	// driving distance function
	Array< float, 3 > G( Q.shape() );
	polar.polar2cartesian(D,G,A);

	writer.setFileName("ipath");
	writer.write(G);

	t.start();

	// initial phi function
	Array< float, 3 > phi( Q.shape() );
	phi = Q;

	// tmp arrays for computations
	Array< float, 3 > R( phi.shape() );
	Array< float, 3 > S( phi.shape() );

	BordStrategyMirrorDouble< float, 3 > BSforPhi( phi, 2 );
	BSforPhi.refresh();

	BordStrategyMirrorDouble< float, 3 > BSforQ( Q, 2 );
	BSforQ.refresh();

	Array< float, 3 > Gx( phi.shape() ), Gy( phi.shape() ), Gz( phi.shape() );

	//BordStrategyMirrorDouble< float, 3 > BSforG( G, 1 );
	//BSforG.refresh();

	Gx = - central12n(Q,firstDim);
	Gy = - central12n(Q,secondDim);
	Gz = - central12n(Q,thirdDim);

	BordStrategyMirrorDouble< float, 3 > BSforGx( Gx, 1 );
	BSforGx.refresh();
	BordStrategyMirrorDouble< float, 3 > BSforGy( Gy, 1 );
	BSforGy.refresh();
	BordStrategyMirrorDouble< float, 3 > BSforGz( Gz, 1 );
	BSforGz.refresh();

	phi = Gx * Gx + Gy * Gy + Gz * Gz;

	for ( int i = 0; i < 5; i++ ){
		R = Gx + .15 * ( Laplacian3D(Gx) - G * ( Gx - central12n(Q,firstDim) ) * phi );
		Gx = R;
		BSforGx.refresh();
		R = Gy + .15 * ( Laplacian3D(Gy) - G * ( Gy - central12n(Q,secondDim) ) * phi );
		Gy = R;
		BSforGy.refresh();
		R = Gz + .15 * ( Laplacian3D(Gz) - G * ( Gz - central12n(Q,thirdDim) ) * phi );
		Gz = R;
		BSforGz.refresh();
	}

	// re-normalize
	phi = sqrt( Gx * Gx + Gy * Gy + Gz * Gz );
	Gx = where( phi > 0, Gx / phi, 0 );
	Gy = where( phi > 0, Gy / phi, 0 );
	Gz = where( phi > 0, Gz / phi, 0 );

	// parameters
	float timeStep = .25;	// evolution time step
	int iterations = 20;	// 20 evolution iterations

	double totalReinit;
	totalReinit = 0;
	Timer tReinit;

	//W = where( W < 1, 1-W, 0 );
	//cout << "W = [" << min(W) << "," << max(W) << "]" << endl;

	phi = Q;
	BSforPhi.refresh();

	for ( int i = 0; i < iterations; i++ ){

		R = 1.0 * curvature3D(phi) - 0.25 * G * advection3D(phi,Gx,Gy,Gz);

		phi = phi + timeStep * R;

		BSforPhi.refresh();

		// reinitialization
		if ( i % 10 == 0 ){
			int iterMorel = 10;
			tReinit.start();

			// compute distances
			S = sussmanSign( phi );
			double stop = 1;
			do {
				R = morelSussman2D(phi,S);
				stop = abs( max(R) );
				phi = phi + 0.1 * R;
				BSforPhi.refresh();
				iterMorel--;
			} while( ( iterMorel > 0 ) & ( stop > 1 ) );

			tReinit.stop();
			totalReinit += tReinit.elapsedSeconds();
		}
		cout << "iter = " << i << endl;
	}

	Q = phi;
	//I = phi;

	t.stop();
	cout << "\tpde regularization = " << t.elapsedSeconds() << "(secs)" << endl;
	cout << "\tredistancing = " << tReinit.elapsedSeconds() << "(secs)" << endl;
	//total += t.elapsedSeconds();

	writer.setFileName( argc[9] );
	writer.write(Q);

	//cout << "Total Algorithm Time = " << total << "(secs)" << endl;

	phi = where( phi < 0, 1, 0 );
	int inside = sum( phi );

	//cout << "Inside voxels = " << inside << " (of " << phi.rows()*phi.cols()*phi.depth() << ")" << endl;

	//cout << argc[10] << endl;
	reader.setFileName( argc[10] );
	reader.read( V );
	//cout << V.shape() << endl;

	// compute volume statistics
	Array< float, 2 > Ip( inside, 1 ); // store only points inside the volume
	Array< float, 2 > :: iterator iIp = Ip.begin();
	Array< float, 3 > :: iterator iterPhi = phi.begin(),
		iterV = V.begin();
	while( iterPhi != phi.end() ){
		if ( (*iterPhi) == 1 ){ // if inside volume
            *iIp = *iterV;
			++iIp;
		}
		++iterPhi;
		++iterV;
	}
	//cout << "Mean gray value = " << sum( Ip ) / inside << ", [" << min(V) << "," << max(V) << "]" << endl;

	// save to disk
	writer.setFileName( argc[11] );
	writer.write( Ip );
	
	writer.setFileName( argc[9] );
	writer.write(Q);
}