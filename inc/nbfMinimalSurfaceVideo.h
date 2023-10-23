#ifndef FILE_nbfMinimalSurfaceVideo
#define FILE_nbfMinimalSurfaceVideo

#include <nbfGeodesicPath.h>
#include <nbfPolarDomain.h>
#include <nbfEdgeFilter.h>

#include <vtkPriorityQueue.h>

using namespace blitz;

template< class Pixel >
class nbfMinimalSurfaceVideo
{
public:

	// constructor takes weight array as input
	nbfMinimalSurfaceVideo(){};

	~nbfMinimalSurfaceVideo(){};

	// convert sequence to cylindrical domain (assumes axis already specified)
	void toCylindrical( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< bool, 3 > & );

	// convert sequence from cylindrical domain to original sequence
	void toCartesian( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< bool, 3 > & );

	// minimal surface computation
	void search( Array< Pixel, 3 > &, Array< Pixel, 3 > & );

	// specify central axis as sequence of 3D points
	void addPointToAxis( TinyVector< Pixel, 3 > & );

	vector< TinyVector< Pixel, 3 > > centralAxis;

protected:

	nbfPolarDomain< Pixel, 2 > polar;

};


template< class Pixel >
void nbfMinimalSurfaceVideo< Pixel > :: addPointToAxis( TinyVector< Pixel, 3 > & point )
{
	this->centralAxis.push_back( point );
}

// ICIP-05 paper version


template< class Pixel >
void nbfMinimalSurfaceVideo< Pixel > :: toCylindrical( Array< Pixel, 3 > & I, Array< Pixel, 3 > & P, Array< bool, 3 > & B )
{
	// check if central axis already specified
	if ( this->centralAxis.size() > 1 ){
		
		float rho = 100;
		this->polar.setMaxRho( rho );
		this->polar.setResRho( rho );
		this->polar.setResTheta( 180 );

		P.resize( this->polar.getMaxRho(), this->polar.getResTheta(), I.depth() );
		B.resize( P.shape() );

		Array< Pixel, 2 > P2;
		Array< bool, 2 > B2;

		// for all line elements of medial axis
		for ( int i = 0; i < this->centralAxis.size() - 1; i++ ){

			// retrieve current start and end points
			TinyVector< Pixel, 3 > start( this->centralAxis[i] );
			TinyVector< Pixel, 3 >   end( this->centralAxis[i+1] );

			// for all slices in between
			for ( int k = start[thirdDim]; k <= end[thirdDim]; k++ ){

				// compute current center
				float d = ( k - start[thirdDim] ) / ( end[thirdDim] - start[thirdDim] );
				TinyVector< Pixel, 2 > center = 
					( 1 -  d ) * TinyVector< Pixel, 2 >( start[firstDim], start[secondDim] ) +
					(      d ) * TinyVector< Pixel, 2 >( end[firstDim], end[secondDim] );

				// set center of polar domain
				this->polar.setCenter( center );
				Array< Pixel, 2 > input( I( Range::all(), Range::all(), k ) );
				this->polar.cartesian2polar( input, P2, B2 );
				P( Range::all(), Range::all(), k ) = P2;		
				B( Range::all(), Range::all(), k ) = B2;		
			}
		}

	}
}

template< class Pixel >
void nbfMinimalSurfaceVideo< Pixel > :: toCartesian( Array< Pixel, 3 > & P, Array< Pixel, 3 > & I, Array< bool, 3 > & B )
{
	Array< Pixel, 2 > I2;
	Array< bool, 2 > B2;

	// for all line elements of medial axis
	for ( int i = 0; i < this->centralAxis.size() - 1; i++ ){

		// retrieve current start and end points
		TinyVector< Pixel, 3 > start( this->centralAxis[i] );
		TinyVector< Pixel, 3 >   end( this->centralAxis[i+1] );

		// for all slices in between
		for ( int k = start[thirdDim]; k <= end[thirdDim]; k++ ){

			// compute current center
			float d = ( k - start[thirdDim] ) / ( end[thirdDim] - start[thirdDim] );
			TinyVector< Pixel, 2 > center = 
				( 1 -  d ) * TinyVector< Pixel, 2 >( start[firstDim], start[secondDim] ) +
				(      d ) * TinyVector< Pixel, 2 >( end[firstDim], end[secondDim] );

			// set center of polar domain
			this->polar.setCenter( center );
			Array< Pixel, 2 > input( P( Range::all(), Range::all(), k ) );
			this->polar.polar2cartesian( input, I2, B2 );
			I( Range::all(), Range::all(), k ) = I2;		
			B( Range::all(), Range::all(), k ) = B2;		
		}
	}
}

template< class Pixel >
void nbfMinimalSurfaceVideo< Pixel > :: search( 
	Array< Pixel, 3 > & input,
	Array< Pixel, 3 > & output )
{
	nbfMatlabWriter writer;
	writer.setFileName("ipath");

	// store rectified 3D polar transformed domain
	Array< Pixel, 3 > P;
	Array< bool, 3 > B;

	// transform to cylindrical
	this->toCylindrical(output,P,B);
	output.resize( P.shape() );
	output = P;
	this->toCylindrical(input,P,B);

	//writer.write(output);
 	writer.write(P);

	BordStrategyMirrorSimple< float,3 > bsForP( P, 1 );
	bsForP.refresh();

	Array< float, 3 > Wy( P.shape() );
	Array< float, 3 > Wz( P.shape() );

	// mini -
	//			Laplacian iters  = 0;
	//			Laplacian lambda = DC;
	//			Wz = where( Wz > 0, 0, 1 ) + 1.0e-1;
	//			Wz( Range(fromStart,6), .... )

	// fedex -
	//			Laplacian iters  = 4;
	//			Laplacian lambda = .108;
	//			Wz = where( Wz > 0, 0, 2 ) + 1.1e0;
	//			Wz( Range(fromStart,9), .... )

	// liron -
	//			Laplacian iters  = 4;
	//			Laplacian lambda = .108;
	//			Wz = where( Wz > 0, 0, 2 ) + 1.1e0;
	//			Wz( Range(fromStart,9), .... )

	// gaussian smoothing
	for ( int i = 0; i < 4; i++ ){
		Wz = Laplacian2D(P);
		P = P + .108 * Wz;
		//P = P + .108 * Wz;
		bsForP.refresh();
	}

	nbfEdgeFilter< Pixel, 3 > edges( P );

	// compute periodic geodesics on individual movie frames

	edges.execute(Wz,thirdDim);
	//writer.write(Wz);

	// 1.0e-0; fedexStart
	// 5.0e-1; fedexEnd
	// 0.9e-1; miniStart
	// 0.9e-1; miniEnd
	//		 ; lironStart	
	//       ; lironEnd

	Pixel wzw = 1e-2;

	Pixel toR = 2;
	Wz = where( Wz > 0, 0, 2 ) + wzw;
	Wz( Range(fromStart,toR), Range::all(), Range::all() ) = numeric_limits<Pixel>::max();
	Wz = where( B, Wz, numeric_limits<Pixel>::max() );

	vector< TinyVector< Pixel, 2 > > path;

	vtkPriorityQueue * pqueueZ = vtkPriorityQueue::New();
	pqueueZ->Allocate( P.depth() );

	for ( int k = 1; k < P.ubound(thirdDim); k++ ){
	//for ( int k = 0; k < P.extent(thirdDim); k+=P.ubound(thirdDim) ){
		Array< Pixel, 2 > slice;
		slice.reference( Wz(  Range::all(), Range::all(), k ) );
		nbfGeodesicPath< Pixel > geodesic( slice );
		Pixel d = geodesic.getCircularPath( slice, path );
		if ( ( k == 0 ) || ( k == P.ubound(thirdDim) ) ){
			geodesic.getImplicitPath( path, output( Range::all(), Range::all(), k ) );
			writer.write( slice );
			writer.write( output( Range::all(), Range::all(), k ) );
		}
		else{
			pqueueZ->Insert( d / path.size(), k );
		}
	}

	//output( Range::all(), Range::all(), 1 ) = output( Range::all(), Range::all(), 0 );
	//output( Range::all(), Range::all(), output.ubound(thirdDim) - 1 ) = output( Range::all(), Range::all(), output.ubound(thirdDim) );
	//this->toCartesian(output,input,B);
	//writer.write(input);

	//Array< Pixel, 2 > smooth( input.rows(), input.cols() );
	//BordStrategyMirrorSimple< float, 2 > bsForS( smooth, 1 );
	//bsForS.refresh();
	//Array< Pixel, 2 > smoothTmp( input.rows(), input.cols() );
	//smooth = input( Range::all(), Range::all(), input.ubound(thirdDim) );
	//smoothTmp = smooth + .25 * curvature2D( smooth );
	//input( Range::all(), Range::all(), input.ubound(thirdDim) ) = smoothTmp;
	//writer.write(input);

	//this->toCartesian(Wz,input,B);
	//writer.write(input);
	//return;

	// compute geodesics across time dimension

	edges.execute(Wy,secondDim);

	Wy = where( Wy > 0, 0, 2 ) + wzw;
	Array< Pixel, 2 > fix1( Wy( Range::all(), Range::all(), Wy.ubound(thirdDim) ) );
	Array< Pixel, 2 > fix2( output( Range::all(), Range::all(), Wy.ubound(thirdDim) ) );
	fix1 = where( fabs( fix2 ) < 1, fix1, numeric_limits< Pixel > :: max() );
	Array< Pixel, 2 > fix3( Wy( Range::all(), Range::all(), Wy.lbound(thirdDim) ) );
	Array< Pixel, 2 > fix4( output( Range::all(), Range::all(), Wy.lbound(thirdDim) ) );
	fix3 = where( fabs( fix4 ) < 1, fix3, numeric_limits< Pixel > :: max() );

	Wy( Range(fromStart,toR), Range::all(), Range::all() ) = numeric_limits<Pixel>::max();
	Wy = where( B, Wy, numeric_limits<Pixel>::max() );
	//writer.write(Wy);

	vtkPriorityQueue * pqueueY = vtkPriorityQueue::New();
	pqueueY->Allocate( P.cols() );

	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;

	for ( int k = 0; k < P.extent(secondDim); k++ ){
		Array< Pixel, 2 > slice, distance;
		slice.reference( Wy(  Range::all(), k, Range::all() ) );
		distance.resize( slice.shape() );
		nbfFastMarching2D8< Pixel > fm( slice );
		aliveP.clear(); aliveD.clear();
		TinyVector< int, 1 > start = minIndex( fabs( output( Range::all(), k, output.lbound(thirdDim) ) ) );
		aliveP.push_back( TinyVector< int, 2 >( start[0], output.lbound(thirdDim) ) ); aliveD.push_back(0);
		fm.setAliveSet( aliveP, aliveD );
		start = minIndex( fabs( output( Range::all(), k, output.ubound(thirdDim) ) ) );
		TinyVector< int, 2 > endP( start[0], output.ubound(thirdDim) );
		fm.setStopPoint( endP );
		fm.execute(distance);
		//writer.write(slice);
		//writer.write(distance);
		nbfGeodesicPath< Pixel > geodesic(distance);
		geodesic.getSimplePath( distance, endP, path );
		pqueueY->Insert( distance(endP) / path.size(), k );
	}

	//writer.write( output( Range::all(), Range::all(), output.lbound(thirdDim) ) );
	//writer.write( output( Range::all(), Range::all(), output.ubound(thirdDim) ) );

	Pixel lowBound = min(P);

	//writer.write(output);
	//return;

	// process all slices in increasing priority order

	double distanceY = numeric_limits<Pixel>::max();
	double distanceZ = numeric_limits<Pixel>::max();
	while( ( pqueueY->GetNumberOfItems() > 0 ) | ( pqueueZ->GetNumberOfItems() > 0 ) ){

		int kY, kZ;

		if ( ( pqueueY->GetNumberOfItems() > 0 ) & ( distanceY == numeric_limits<Pixel>::max() ) ){
			kY = pqueueY->Pop(0,distanceY);
		}

		if ( ( pqueueZ->GetNumberOfItems() > 0 ) & ( distanceZ == numeric_limits<Pixel>::max() ) ){
			kZ = pqueueZ->Pop(0,distanceZ);
		}

		// compute geodesic across time frames
		if ( distanceY < distanceZ ){
			cout << "y = " << kY << endl;
			Array< Pixel, 2 > slice, distance, ipath;
			slice.reference( Wy( Range::all(), kY, Range::all() ) );
			ipath.reference( output( Range::all(), kY, Range::all() ) );
			distance.resize( slice.shape() );
			nbfFastMarching2D8< Pixel > fm( slice );
			TinyVector< int, 1 > start = minIndex( fabs( output( Range::all(), kY, output.lbound(thirdDim) ) ) );
			aliveP.clear(); aliveD.clear();
			aliveP.push_back( TinyVector< int, 2 >( start[0], output.lbound(thirdDim) ) ); aliveD.push_back(0);
			fm.setAliveSet( aliveP, aliveD );
			start = minIndex( fabs( output( Range::all(), kY, output.ubound(thirdDim) ) ) );
			TinyVector< int, 2 > endP( start[0], output.ubound(thirdDim) );
			fm.setStopPoint( endP );
			fm.execute(distance);
			nbfGeodesicPath< Pixel > geodesic(distance);
			geodesic.getSimplePath( distance, endP, path );
			geodesic.getImplicitPath( path, ipath, 1.0 );
			Range I = Range::all();
			Range J( fromStart + 1, toEnd - 1 );
			slice(I,J) = where( fabs(ipath(I,J)) < 1, fabs(ipath(I,J)), 1 ) + 1 * lowBound;
			
			Wz( I, kY, J ) = slice;
			distanceY = numeric_limits<Pixel>::max();
		}
		else{ // compute geodesic on given time frame

			cout << "z = " << kZ << endl;
			Array< Pixel, 2 > slice, ipath;
			slice.reference( Wz(  Range::all(), Range::all(), kZ ) );
			ipath.reference( output( Range::all(), Range::all(), kZ ) );
			nbfGeodesicPath< Pixel > geodesic( slice );
			geodesic.getCircularPath( slice, path );
			geodesic.getImplicitPath( path, ipath );
			slice = where( fabs(ipath) < 1, fabs(ipath), 1 ) + 1 * lowBound;
			Wy( Range::all(), Range::all(), kZ ) = slice;

			//writer.write(slice);
			distanceZ = numeric_limits<Pixel>::max();
		}
	}

	//writer.write(output);
	//this->toCartesian(output,input,B);
	//writer.write(input);
	//return;

	// second pass to guarantee smoothness

	for ( int k = 0; k < P.extent(secondDim); k++ ){
		Array< Pixel, 2 > slice, distance, ipath;
		slice.reference( Wy(  Range::all(), k, Range::all() ) );
		ipath.reference( output(  Range::all(), k, Range::all() ) );
		distance.resize( slice.shape() );
		nbfFastMarching2D8< Pixel > fm( slice );
		TinyVector< int, 1 > start = minIndex( fabs( output( Range::all(), k, output.lbound(thirdDim) ) ) );
		aliveP.clear(); aliveD.clear();
		aliveP.push_back( TinyVector< int, 2 >( start[0], output.lbound(thirdDim) ) ); aliveD.push_back(0);
		fm.setAliveSet( aliveP, aliveD );
		start = minIndex( fabs( output( Range::all(), k, output.ubound(thirdDim) ) ) );
		TinyVector< int, 2 > endP( start[0], output.ubound(thirdDim) );
		fm.setStopPoint( endP );
		fm.execute(distance);
		nbfGeodesicPath< Pixel > geodesic(distance);
		geodesic.getSimplePath( distance, endP, path );
		geodesic.getImplicitPath( path, ipath );
		slice = lowBound + fabs(ipath);
		Wz( Range::all(), k, Range::all() ) = slice;
	}

	for ( int k = 0; k < P.extent(thirdDim); k++ ){
		Array< Pixel, 2 > slice, ipath;
		slice.reference( Wz(  Range::all(), Range::all(), k ) );
		ipath.reference( output( Range::all(), Range::all(), k ) );
		nbfGeodesicPath< Pixel > geodesic( slice );
		geodesic.getCircularPath( slice, path );
		geodesic.getImplicitPath( path, ipath );
	}

	//writer.write(output);

	//return;

	B.resize( input.shape() );
	this->toCartesian(output,input,B);
	//output.resize( input.shape() );
	//output = input;
	//writer.write(input);

	return;
}

#endif /* FILE_nbfMinimalSurfaceVideo */