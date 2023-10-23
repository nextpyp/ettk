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

	// multi purpose writer
	nbfMatlabWriter writer;
	writer.setFileName( "ipath" );

	//Array< float, 2 > pistola(160,160);
	//pistola = 1;
	Timer t;

	// 1. V - input image
	Array< float, 3 > V;

	reader.setFileName( argc[1] );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	//BordStrategyMirrorSimple< float, 3 > bfForV(V,1);
	//bfForV.refresh();
	Array< float, 3 > G( V.shape() );
	//G = sqrt( pow2( central12(V,firstDim) ) + pow2( central12(V,secondDim) ) + pow2( central12(V,thirdDim) ) );

	//// set to [0-1]
	//G = G - min(G);
	//G = G / max(G);
	//G = 1 - G + 1e-5;
	////G = where( G > .99, 10 * G, G );

	cout << max(V) << "," << min(V) << endl;

	G = V - min(V);
	G = G / max(G);
	//writer.write(G);
	G = 1*(1 - G) + 1e-2;
	cout << max(G) << "," << min(G) << endl;

	//// 1. V - input image
	//Array< float, 2 > Im;

	//reader.setFileName( argc[1] );
	//if ( reader.read( Im ) ){
	//	cout << "Error reading file.\n";
	//	exit(0);
	//}
	//Im = 5 * ( 1 - Im ) + 1e-5;
	//cout << min(Im) << max(Im) << endl;

	// retrieve center
	float centerX = atof(argc[2]);
	float centerY = atof(argc[3]);
	float centerZ = atof(argc[4]);
	TinyVector< float, 3 > center( centerX, centerY, centerZ );

	Array< float, 2 > P;
	Array< float, 2 > W;
	Array< bool, 2 >  B;

	Array< float, 3 > implicit( V.shape() );

#if 1

	// slice domain in direction perpendicular to z
	nbfPolarDomain< float, 2 > polar;
	float rho = 300;
	polar.setCenter( TinyVector< float, 2 >( centerX, centerY ) );
	polar.setMaxRho( rho );
	polar.setResRho( rho );
	polar.setResTheta( 180 );

	// retrieve center
	float cX1 = atof(argc[2]);
	float cY1 = atof(argc[3]);
	float cZ1 = atof(argc[4]);
	float cX2 = atof(argc[5]);
	float cY2 = atof(argc[6]);
	float cZ2 = atof(argc[7]);

	for ( float k =  0; k <= cZ2; k++ ){
		cout << k << endl;

		//int k = 0;
		//nbfCutGeodesics< float > geodesicP( G( Range::all(), Range::all(), k ), TinyVector< int, 2 >(centerX,centerY) );
		vector< TinyVector< float, 2 > > path;
		//geodesicP.getFullCircularPath( G( Range::all(), Range::all(), k ), path, true, 0.0 );
		//nbfCutGeodesics< float > geodesicP( Im, TinyVector< int, 2 >(centerX,centerY) );
		//geodesicP.getFullCircularPath( Im, path, true, 0.0 );
		//Array< float, 2 > ipath( Im.shape() );
		Array< float, 2 > ipath( V.rows(), V.cols() );
		//geodesicP.getImplicitPath(path,ipath);
		//writer.write(G( Range::all(), Range::all(), k ));
		//writer.write(ipath);

		//Array< float, 2 > current( Im.shape() );
		//Array< float, 2 > current( V.rows(), V.cols() );
		//current = G( Range::all(), Range::all(), k );
		//polar.cartesian2polar( current, P );
		//polar.cartesian2polar( Im, 0 ), P );

		float centerX = ( 1 - k / cZ2 ) * cX1 + k / cZ2 * cX2;
		float centerY = ( 1 - k / cZ2 ) * cY1 + k / cZ2 * cY2;
		polar.setCenter( TinyVector< float, 2 >( centerX, centerY ) );

		polar.cartesian2polar( G( Range::all(), Range::all(), k ), P );
		P( Range(0,10), Range::all() ) = numeric_limits<float>::max();
		//writer.write(G( Range::all(), Range::all(), 0 ));
		//writer.write(P);

		nbfGeodesicPath< float > geodesic(P);
		float d = geodesic.getCircularPath(P,path);
		//cout << k << " - d = " << d << ", r = " << ( d - k * path.size() * .25) / path.size() << endl;

		Array< float, 2 > I( P.shape() );
		geodesic.getImplicitPath(path,I);
		//writer.write(I);

		ipath.resize( V.rows(), V.cols() );
		polar.polar2cartesian(I,ipath);
		//writer.write(ipath);
		implicit( Range::all(), Range::all(), k ) = ipath;
	}


	//nbfMinimalSurface< float > ms(implicit);
	//Array< bool, 3 > Bo( implicit.shape() );
	//ms.computeConfidence(implicit,thirdDim,2,Bo);
	//writer.write(implicit);

	//implicit = where( fabs(implicit) < 1.0, fabs(implicit), 1.0 ) + 1e0;
	//writer.setFileName("pre");
	writer.write(implicit);
	//writer.write(Bo);

#else

	for ( int i = 0; i < V.depth(); i++ ){
	cout << "i = " << i << " - ";
	//polar.cartesian2polar( V( Range::all(), Range::all(), i ), P );	

	P.reference( V( Range(fromStart,49), Range::all(), i ) );
	Array< float, 2 > distances( P.shape() );

	nbfFastMarching2D< float > fm( P );
	vector< TinyVector< int, 2 > > aliveP;
	vector< float > aliveD;

	float minDistance = numeric_limits<float>::max();
	TinyVector< int, 2 > minPoint;

	for ( int k = 0; k < P.rows(); k++ ){
		aliveP.clear();
		aliveD.clear();
		aliveP.push_back( TinyVector< int, 2 >(k,0) );
		aliveD.push_back(0);
		fm.setAliveSet( aliveP, aliveD );
		fm.setStopPoint( TinyVector< int, 2 >(k,P.ubound(secondDim)) );
		TinyVector< int, 2 > last = fm.execute(distances);
		if ( distances(last) < minDistance ){
			minDistance = distances(last);
			minPoint = last;
		}
	}

	// cout << minPoint << endl;

	aliveP.clear();
	aliveD.clear();
	
	// set alive (d = 0) all points on the left side of the image
	for ( int i = 0; i < P.rows(); i++ )
	{
		TinyVector< int, 2 > current( i, 0 );
		aliveP.push_back(current);
		aliveD.push_back(0);
	}			

	fm.setAliveSet( aliveP, aliveD );

	// compute distances
	TinyVector< int, 2 > first = fm.execute( distances );

	// new way of computing circular paths

	Array< float, 1 > rightLine( distances.rows() );
	rightLine = distances( Range::all(), distances.ubound(secondDim) );

	aliveP.clear();
	aliveD.clear();
	for ( int i = 0; i < P.rows(); i++ )
	{
		TinyVector< int, 2 > current( i, P.ubound(secondDim) );
		aliveP.push_back(current);
		aliveD.push_back(0);
	}			

	fm.setAliveSet( aliveP, aliveD );
	fm.execute(distances);

	rightLine += distances( Range::all(), 0 );
	TinyVector< int, 1 > point = minIndex( rightLine );

	cout << minPoint[0] << " ? " << point[0] << " ==> ";

	if ( minPoint[0] == point[0] ){
        cout << " success " << endl;
	}
	else{
		cout << " *** ERROR *** " << endl;
	}
	}

#endif

	////W.resize( P.shape() );
	////TinyVector< int, 2 > end = fm.execute(W);
	////writer.write(W);

	////nbfGeodesicPath< float > geodesic(W);
	////vector< TinyVector< float, 2 > > path;
	////geodesic.getSimplePath(W,end,path);

	////Array< float, 2 > I( W.shape() );
	////geodesic.getImplicitPath(path,I);
	////writer.write(I);

	//// second pass
	//W.resize( P.shape() );
	//points.clear();
	//distances.clear();
	//for ( int i = 0; i < P.rows(); i++ ){
	//	points.push_back( TinyVector< int, 2 >(i,W.ubound(secondDim)) );
	//	distances.push_back(0.0);
	//}
	//fm.setAliveSet(points,distances);
	//TinyVector< int, 2 > end = fm.execute(W);
	//TinyVector< int, 1 > point = minIndex( W( Range::all(), 0 ) );
	//end = TinyVector< int, 2 >( point[0], W.ubound(secondDim) );
	//cout << end << endl;
	//writer.write(W);
	//W.reverseSelf(secondDim);
	//writer.write(W);

	//nbfGeodesicPath< float > geodesic(W);
	//vector< TinyVector< float, 2 > > path;
	//geodesic.getSimplePath(W,end,path);
	//Array< float, 2 > I( W.shape() );
	//geodesic.getImplicitPath(path,I);
	//I.reverseSelf(secondDim);
	//writer.write(I);

	//int minI = 0;
	//float minD = numeric_limits<float>::max();
	//for ( int i = 0; i < W.rows(); i++ ){
	//	points.clear();
	//	distances.clear();
	//	points.push_back( TinyVector< int, 2 >(i,0) );
	//	distances.push_back(0);
	//	fm.setAliveSet(points,distances);
	//	TinyVector< int, 2 > end = fm.execute(W);
	//	float d = W(TinyVector< int, 2 >(i,W.ubound(secondDim)));
	//	if ( d < minD ){
	//		minD = d;
	//		minI = i;
	//	}
	//	cout << i << " - " <<  d << endl;
	//}
	//cout << minI << " - " << minD << endl;

	////

	//polar.polar2cartesian(I,W);
	//writer.write(W);

	//V = V3( Range::all(), Range::all(), 19 );
	//V = V - min(V);
	//V = 1e0 * V + 1e-2;

}