#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#define PIXEL float

#define NBZ_DEBUG

#include <nbfTimer.h>

#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <vtkMath.h>
#include <cut/nbfCutMinimalSurface.h>
#include <nbfDifferentials.h>
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

	Array< PIXEL, 3 > V, W, P, Q, D, I;

	nbfMatlabReader reader;
	nbfMatlabWriter writer;

	Timer t;
	double total = 0;

	reader.setFileName( argc[1] );
	if ( reader.read( W ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	// B2
	int centerX = atoi(argc[2]);
	int centerY = atoi(argc[3]);
	int centerZ = atoi(argc[4]);
	TinyVector< int, 3 > center( W.ubound(firstDim) / 2, W.ubound(secondDim) / 2, centerZ );

	 nbfCutMinimalSurface< PIXEL > ms(W,center);

#if 1 // recompute D or read from file
	reader.setFileName( argc[8] );
	if ( reader.read( D ) ){
		cout << "Error reading file.\n";
		exit(0);
	}
	cout << "D = " << D.shape() << ", [" << min(D) << "," << max(D) << "]" << endl;
#else
	t.start();
	//ms.skeleton(D);
	ms.getDistance(D);
	t.stop();
	cout << "\t3d distances = "  << t.elapsedSeconds() << " seconds." << endl;
	total += t.elapsedSeconds();
	writer.setFileName( argc[8] );
	writer.write(D);
#endif

#if 1
	t.start();
	ms.execute(I,D);
	t.stop();
	cout << "\tminimal surface = "  << t.elapsedSeconds() << " seconds." << endl;
	total += t.elapsedSeconds();

	writer.setFileName( argc[7] );
	writer.write(I);
#else
	reader.setFileName( argc[7] );
	if ( reader.read( I ) ){
		cout << "Error reading file.\n";
		exit(0);
	}
#endif

	// REGULARIZATION
#if 1
	// run pde

	// tmp arrays for computations
	Array< PIXEL, 3 > phi( I.shape() );
	Array< PIXEL, 3 > G( phi.shape() );
	Array< PIXEL, 3 > R( phi.shape() );
	Array< PIXEL, 3 > S( phi.shape() );

	reader.setFileName( argc[1] );
	if ( reader.read( W ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	BordStrategyMirrorDouble< PIXEL, 3 > BSforW( W, 1 );
	BSforW.refresh();

	t.start();

	// build vector field as in gvf
	Array< PIXEL, 3 > Gx( phi.shape() ), Gy( phi.shape() ), Gz( phi.shape() );
	Gx = central12n(W,firstDim);
	Gy = central12n(W,secondDim);
	Gz = central12n(W,thirdDim);

	BordStrategyMirrorDouble< PIXEL, 3 > BSforGx( Gx, 1 );
	BSforGx.refresh();
	BordStrategyMirrorDouble< PIXEL, 3 > BSforGy( Gy, 1 );
	BSforGy.refresh();
	BordStrategyMirrorDouble< PIXEL, 3 > BSforGz( Gz, 1 );
	BSforGz.refresh();

	phi = Gx * Gx + Gy * Gy + Gz * Gz;

	for ( int i = 0; i < 5; i++ ){
		R = Gx + .15 * ( Laplacian3D(Gx) - ( 1 - W ) * ( Gx - central12n(W,firstDim) ) * phi );
		Gx = R;
		BSforGx.refresh();
		R = Gy + .15 * ( Laplacian3D(Gy) - ( 1 - W ) * ( Gy - central12n(W,secondDim) ) * phi );
		Gy = R;
		BSforGy.refresh();
		R = Gz + .15 * ( Laplacian3D(Gz) - ( 1 - W ) * ( Gz - central12n(W,thirdDim) ) * phi );
		Gz = R;
		BSforGz.refresh();
	}

	// re-normalize
	phi = sqrt( Gx * Gx + Gy * Gy + Gz * Gz );
	Gx = where( phi > 0, Gx / phi, 0 );
	Gy = where( phi > 0, Gy / phi, 0 );
	Gz = where( phi > 0, Gz / phi, 0 );

	t.stop();
	cout << "\tGVF field = "  << t.elapsedSeconds() << " seconds." << endl;
	total += t.elapsedSeconds();

	t.start();

	phi = I;

	BordStrategyMirrorDouble< PIXEL, 3 > BSforPhi( phi, 2 );
	BSforPhi.refresh();

	// parameters
	float timeStep = .15;	// evolution time step
	int iterations = 15;	// 20 evolution iterations

	double totalReinit;
	totalReinit = 0;
	Timer tReinit;

	for ( int i = 0; i < iterations; i++ ){

		// reinitialization
		if ( (i % 10 == 0) ){
			int iterMorel = 10; // 10
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

		//		R = 1.0 * curvature3D(phi) - 0.25 * W * advection3D(phi,Gx,Gy,Gz);
		R = 1.0 * curvature3D(phi) + 0.25 * advection3D(phi,Gx,Gy,Gz);

		phi = phi + timeStep * R;

		BSforPhi.refresh();
	}

	t.stop();
	cout << "\tpde regularization = " << t.elapsedSeconds() << "(secs)" << endl;
	cout << "\tredistancing = " << tReinit.elapsedSeconds() << "(secs)" << endl;
	total += t.elapsedSeconds();

	cout << "Total Algorithm Time = " << total << "(secs)" << endl;

	writer.setFileName( argc[9] );
	writer.write(phi);

	phi = where( phi < 0, 1, 0 );
	int inside = sum( phi );

	//cout << "Inside voxels = " << inside << " (of " << phi.rows()*phi.cols()*phi.depth() << ")" << endl;

	reader.setFileName( argc[10] );
	reader.read( V );

	// compute volume statistics
	Array< PIXEL, 2 > Ip( inside, 1 ); // store only points inside the volume
	Array< PIXEL, 2 > :: iterator iIp = Ip.begin();
	Array< PIXEL, 3 > :: iterator iterPhi = phi.begin(),
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
	
#endif
}