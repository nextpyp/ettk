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
#include <nbfMinimalSurface.h>

#define PIXEL float
#define DIM 3

void main( int argv, char ** argc )
{
	if ( argv != 7 ){
		cout << "Usage: Distancer.exe inputFile.array outputFile.array initialPoint" << endl;
		cout << "\t inputFile.array : input file name (from matlab)" << endl;
		cout << "\t outputFile.array : output file name (to be read in matlab)" << endl;
		cout << "\t initial point : select one of the four corners 1-4" << endl;
		exit(0);
	}

	Array< PIXEL, 3 > V, W, P, Q, D, I;

	nbfMatlabReader reader;
	nbfMatlabWriter writer;

	Timer t;
	double total;
	total = 0;

#if 0 // recompute weights W or read from file (this is done once for each 3d set)
	reader.setFileName( argc[3] );
	if ( reader.read( W ) ){
		cout << "Error reading file.\n";
		exit(0);
	}
	W.transposeSelf(firstDim,thirdDim,secondDim);

	cout << "W = " << W.shape() << endl;
#else
	reader.setFileName( argc[1] );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}
	V.transposeSelf(firstDim,thirdDim,secondDim);
	cout << "V = " << V.shape() << endl;

	BordStrategyMirrorSimple< PIXEL,3 > bsForV( V, 1 );
	bsForV.refresh();

	W.resize( V.shape() );

	t.start();

	// gaussian smoothing
#if 0
	//for ( int i = 1; i < 30; i++ ){
	for ( int i = 1; i < 10; i++ ){
		W = Laplacian3D(V);
		V = V + .1 * W;
		bsForV.refresh();
	}
#else
	nbfGaussianFilter< PIXEL, 3 > gauss( V );
	for ( int i = 0; i < 10; i++ ){
        gauss.execute(W);
		V = W;
		bsForV.refresh();
	}
#endif
	t.stop();
	cout << "gauss blur "  << t.elapsedSeconds() << " seconds." << endl;
	total += t.elapsedSeconds();

	nbfVeselnessFilter< PIXEL, 3 > vesel( V );
	t.start();
	vesel.execute(W);
	t.stop();
	cout << "vesel measure "  << t.elapsedSeconds() << " seconds." << endl;

	// don't need V anymore (free from memory)
	//V.free();

	cout << "W = " << W.shape() << endl;

	// binarize
//	PIXEL th = nbfArrayFilter<PIXEL,3>::graythresh(W);
//	W = where( W > th, th, W ) / th;
//	W = 1 - W;

	W = 1 - sqrt( sqrt( sqrt( sqrt(W) ) ) );

#if 1
	writer.setFileName( argc[3] );
	writer.write(W);
	//exit(0);
#endif

#endif

	nbfPolarDomain< PIXEL, 3 > polar;

	// ill
	polar.setCenter( TinyVector<PIXEL,3>(139,119,49) );
	//int centerX = 139;
	//int centerY = 119;
	//int centerZ = 49;

	// healthy
//	polar.setCenter( TinyVector<PIXEL,3>(69,89,44) );

	//polar.setCenter( TinyVector<PIXEL,3>(86,131,45) );
	//polar.setCenter( TinyVector<PIXEL,3>(81,39,45) );
	//polar.setCenter( TinyVector<PIXEL,3>(43,153,45) );

	// V2
	//polar.setCenter( TinyVector<PIXEL,3>(137,116,50) );

	// B2
	//int centerX = 275;
	//int centerY = 236;
	//int centerZ = 65;

	int rho = 60;
	//int lboundX = blitz::minmax::max( centerX - rho * 3 / 4, W.lbound(firstDim) );
	//int lboundY = blitz::minmax::max( centerY - rho * 3 / 4, W.lbound(secondDim) );
	//int uboundX = blitz::minmax::min( centerX + rho * 3 / 4, W.ubound(firstDim) );
	//int uboundY = blitz::minmax::min( centerY + rho * 3 / 4, W.ubound(secondDim) );

	//cout << "[" << lboundX << "," << uboundX << "]" << endl;
	//cout << "[" << lboundY << "," << uboundY << "]" << endl;

	//Array< PIXEL, 3 > w( uboundX - lboundX + 1, uboundY - lboundY + 1, W.depth() );
	//w = W( Range(lboundX,uboundX), Range(lboundY,uboundY), Range::all() );
	//W.free();

#if 0
	Array< PIXEL, 3 > v( uboundX - lboundX + 1, uboundY - lboundY + 1, V.depth() );
	v = V( Range(lboundX,uboundX), Range(lboundY,uboundY), Range::all() );
	V.free();
	writer.setFileName( argc[1] );
	writer.write(v);
	exit(0);
#endif

	//cout << "w = " << w.shape() << endl;

	//polar.setCenter( TinyVector<PIXEL,3>( w.ubound(firstDim) / 2, w.ubound(secondDim) / 2, centerZ ) );

	polar.setMaxRho( rho );
	polar.setResRho( 120 );
	polar.setResTheta( 2 * vtkMath::Pi() / atanf( 1.0 / rho ) );

	t.start();
	polar.cartesian2polar(W,P);
	
	// avoid small circles
	P( Range(fromStart,40), Range::all(), Range::all() ) = 1;

	t.stop();
	cout << "to polar "  << t.elapsedSeconds() << " seconds." << endl;
	total += t.elapsedSeconds();

	cout << "P = " << P.shape() << endl;

	P = P + 1e-1;

	writer.setFileName( argc[4] );
	writer.write(P);

	nbfMinimalSurface< PIXEL > ms(P);

#if 0 // recompute D or read from file
	reader.setFileName( argc[2] );
	if ( reader.read( D ) ){
		cout << "Error reading file.\n";
		exit(0);
	}
	cout << "D = " << D.shape() << endl;
#else
	t.start();
	ms.getDistance(D);
	t.stop();
	cout << "3d distances "  << t.elapsedSeconds() << " seconds." << endl;
	total += t.elapsedSeconds();
	writer.setFileName( argc[2] );
	writer.write(D);
#endif

	t.start();
	ms.execute(I,D);
	t.stop();
	cout << "minimal surface "  << t.elapsedSeconds() << " seconds." << endl;
	total += t.elapsedSeconds();

	t.start();
	polar.polar2cartesian(I,Q);
	t.stop();
	cout << "back to cartesian "  << t.elapsedSeconds() << " seconds." << endl;
	total += t.elapsedSeconds();

	Q = where( Q < numeric_limits<PIXEL>::max(), Q, - numeric_limits<PIXEL>::max() );
	PIXEL mQ = max(Q);
	Q = where( Q > - numeric_limits<PIXEL>::max(), Q, mQ ); 

	writer.setFileName( argc[5] );
	writer.write(I);

#if 1
	// driving distance function
	Array< float, 3 > G( Q.shape() );
	polar.polar2cartesian(D,G);

#endif

	// REGULARIZATION
#if 1
	// run pde
	t.start();

	// initial phi function
	Array< float, 3 > phi( Q.shape() );
	phi = Q;

	/* build narrow band
	typedef TinyVector<int,3> coord;
	std::vector<coord> Ip;
	Array< float, 3 > :: iterator iter = phi.begin();
	while( iter != phi.end() ){
		if ( fabs(*iter) < 10 ){
            Ip.push_back( iter.position() );
		}
		++iter;
	}

	cout << "narrow band = " << Ip.size() << " of " << phi.rows()*phi.cols()*phi.depth() <<endl;
	*/

	// tmp arrays for computations
	Array< float, 3 > R( phi.shape() );
	Array< float, 3 > S( phi.shape() );

	BordStrategyMirrorDouble< float, 3 > BSforPhi( phi, 2 );
	BSforPhi.refresh();

	Array< float, 3 > Gx( phi.shape() ), Gy( phi.shape() ), Gz( phi.shape() );

	BordStrategyMirrorDouble< float, 3 > BSforG( G, 1 );
	BSforG.refresh();

	Gx = - central12n(G,firstDim);
	Gy = - central12n(G,secondDim);
	Gz = - central12n(G,thirdDim);

	// normalization
	R = sqrt( pow2(Gx) + pow2(Gy) + pow2(Gz) );
	Gx = where( R > 0, Gx / R, Gx );
	Gy = where( R > 0, Gy / R, Gy );
	Gz = where( R > 0, Gz / R, Gz );

	// parameters
	float timeStep = .25;	// evolution time step
	int iterations = 25;	// evolution iterations

	double totalReinit;
	totalReinit = 0;
	Timer tReinit;

	W = where( W < 1, 1-W, 0 );

	for ( int i = 0; i < iterations; i++ ){

//		R = 1.0 * curvature3D(phi) - 0.25 * advection3D(phi,Gx,Gy,Gz);
		R = 1.0 * curvature3D(phi) - 0.25 * W * advection3D(phi,Gx,Gy,Gz);

		phi = phi + timeStep * R;

		BSforPhi.refresh();

		// reinitialization
		if ( i % 10 == 0 ){
			int iterMorel = 10;
			tReinit.start();

			// compute distances
			S = sussmanSign( phi );
			do {
				R = morelSussman2D(phi,S);
				phi = phi + 0.1 * R;
				BSforPhi.refresh();
				iterMorel--;
			} while( iterMorel > 0 );

			tReinit.stop();
			totalReinit += tReinit.elapsedSeconds();
		}

	}

	Q = phi;

	t.stop();
	cout << "pde regularization = " << t.elapsedSeconds() << "(secs)" << endl;
	cout << "redistancing = " << tReinit.elapsedSeconds() << "(secs)" << endl;
	total += t.elapsedSeconds();

	cout << "Total Algorithm Time = " << total << "(secs)" << endl;

	phi = where( phi < 0, 1, 0 );
	int inside = sum( phi );

	cout << "Inside voxels = " << inside << " (of " << phi.rows()*phi.cols()*phi.depth() << ")" << endl;

	Array< float, 1 > Ip( inside );
	Array< float, 1 > :: iterator iIp = Ip.begin();
	Array< float, 3 > :: iterator iterPhi = phi.begin(),
		iterV = V.begin();
	while( iterPhi != phi.end() ){
		if ( (*iterPhi) == 1 ){
            *iIp = *iterV;
			++iIp;
		}
		++iterPhi;
		++iterV;
	}
	cout << "Mean gray value = " << sum( Ip ) / inside << ", [" << min(V) << "," << max(V) << "]" << endl;


#endif
	writer.setFileName( argc[6] );
	writer.write(Q);

}