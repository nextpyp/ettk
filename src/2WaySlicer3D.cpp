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
#include <nbfSliceDomain.h>
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
	writer.setFileName( "ipath" );

	Timer t;

	// 1. V - input image
	Array< float, 3 > V;

	reader.setFileName( argc[1] );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	//reader.setFileName("ipath");
	//reader.read(V);
	//
	//nbfImageWriter vtkw;
	//vtkw.setFileName("cell.vtk");
	//vtkw.write(V);

	BordStrategyMirrorSimple< float,3 > bsForV( V, 1 );
	bsForV.refresh();

	Array< float, 3 > W( V.shape() );

	//// gaussian smoothing
	//for ( int i = 0; i < atoi( argc[6] ); i++ ){
	//	W = Laplacian3D(V);
	//	V = V + .1 * W;
	//	bsForV.refresh();
	//}

	//nbfVeselnessFilter< float, 3 > veselV( V );
	//veselV.execute(W,thirdDim);
	//W = 1 - sqrt( sqrt( sqrt( sqrt(W) ) ) ) + 1e-6;
	//writer.write(W);

	//Array< float, 3 > D( V.shape() );
	//V = 1;
	//nbfFastMarching3D< float > fm3(V);
	//vector< TinyVector< int, 3 > > points;
	//vector< float > distances;
	//Array< float, 3 > :: iterator iter = W.begin();
	//while( iter != W.end() ){
	//	if ( (*iter) < 1 ){
	//		points.push_back( iter.position() );
	//		distances.push_back(*iter);
	//	}
	//	++iter;
	//}
	//fm3.setAliveSet( points, distances );
	//fm3.execute(D);

	//W = D;
	//writer.write(W);

	// B2
	int centerX = atoi(argc[2]);
	int centerY = atoi(argc[3]);
	int centerZ = atoi(argc[4]);
	TinyVector< int, 3 > center( centerX, centerY, centerZ );

	// slice domain in direction perpendicular to z
	nbfSliceDomain< float, 3 > slicer;
	Array< float, 3 > P;
	Array< bool, 3 > B;

	slicer.setCenter( center );

	Array< float, 3 > I;

	V = V - min(V);
	//V = V / max(V) * 10;
	V = V + 1;

	//slicer.slice(V,P,thirdDim,B);
	//writer.write(P);

	slicer.slice(V,P,secondDim,B);
	P = where( B == true, P, numeric_limits<float>::max() );
	writer.write(P);

	TinyVector< int, 3 > center1( P.shape() / 2 );

	nbfMinimalSurface< float > ms1(P);
	//ms1.search(P,center1,I,2.0);
	ms1.search(P,center1,I,numeric_limits<float>::max());
	writer.write(I);

	//// compute confidence from slices
	//ms1.computeConfidence( I, secondDim, 2.0, B );
	////writer.write(B);

#if 1
	slicer.unslice(I,W,secondDim,B);
	slicer.slice(W,P,firstDim,B);
	P = where( fabs(P) < 1.0, fabs(P), 0.0 );
	P = 1 - P + 1e0;
	P = where( B == true, P, numeric_limits<float>::max() );
#else
	//I = where( B == true, 1, 0 );
	I = where( fabs(I) < 1.0, 1, 0 );
	writer.write(I);

	slicer.unslice(I,P,secondDim,B);
	writer.write(P);

	//W = where( ( P > 0 ), W / 10.0, W );
	W = 1 - P + 1e-0;
	//writer.write(W);

	slicer.slice(W,P,firstDim,B);
	P = where( B == true, P, numeric_limits<float>::max() );
	writer.write(P);
#endif

	//writer.write(P);

	TinyVector< int, 3 > center2( P.shape() / 2 );
	nbfMinimalSurface< float > ms2(P);
	//ms2.search(P,center2,I,numeric_limits<float>::max());
	ms2.search(P,center2,I,2.0);

	//ms2.computeConfidence( I, secondDim, 2.0, B );
	//writer.write(B);

	//writer.write(P);
	//writer.write(I);
	slicer.unslice(I,P,firstDim,B);
	//writer.write(P);
	//vtkw.write(P);

	slicer.slice(P,W,secondDim,B);
	W = where( fabs(W) < 1.0, fabs(W), 0.0 );
	W = 1 - W + 1e0;
	W = where( B == true, W, numeric_limits<float>::max() );
	nbfMinimalSurface< float > ms11(W);
	ms11.search(W,center1,I,numeric_limits<float>::max());
	writer.write(W);
	writer.write(I);
	slicer.unslice(I,P,secondDim,B);
	writer.write(P);
}