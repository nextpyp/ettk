#ifndef FILE_nbfMinimalSurface
#define FILE_nbfMinimalSurface

#include <fm/nbfFastMarching3D.h>
#include <fm/nbfFastFastMarching2D.h>
#include <cut/nbfCutGeodesics.h>

#include <map>

#include <vtkPriorityQueue.h>

using namespace blitz;

template< class Pixel >
class nbfMinimalSurface
{
public:

	// constructor takes weight array as input
	nbfMinimalSurface( Array< Pixel, 3 > & );

	~nbfMinimalSurface(){};

	void getDistance( Array< Pixel, 3 > & );
	void execute( Array< Pixel, 3 > &,  Array< Pixel, 3 > & );
	void executeOld( Array< Pixel, 3 > &,  Array< Pixel, 3 > & );
	void executeNew( Array< Pixel, 3 > &,  Array< Pixel, 3 > & );

	void executeUltimate( Array< Pixel, 3 > & );
	
	void fillSlice( Array< Pixel, 2 > &, Array< Pixel, 2 > &, Array< Pixel, 2 > & );

	void search( Array< Pixel, 3 > &, Array< Pixel, 3 > & );
	void searchGUI( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Pixel );
	void search( Array< Pixel, 3 > &, TinyVector< int, 3 > &, Array< Pixel, 3 > &, Pixel );
	void search( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< Pixel, 3 > &,	TinyVector< int, 3 > &,	Array< Pixel, 3 > & );

	void computeGeodesics( Array< Pixel, 3 > &,	int, int, TinyVector< int, 3 > &, Array< Pixel, 3 > & );
	void computeGeodesics( Array< Pixel, 3 > &, int, TinyVector< int, 2 > &, Array< Pixel, 3 > & );

	void computeConfidence( Array< Pixel, 3 > &, int, int, Array< bool, 3 > & );
	TinyVector< int, 2 > updateBaricenter( TinyVector< int, 2 > &, vector< TinyVector< Pixel, 2 > > & );

	void computeGeodesic( Array< Pixel, 2 > &, TinyVector< int, 2 > &, TinyVector< int, 2 > &, Array< Pixel, 2 > &, Pixel = numeric_limits<Pixel>::max() );

	void setGivenDistances(){ this->givenDistances = 1; };
	void unSetGivenDistances(){ this->givenDistances = 0; };

protected:

	bool givenDistances;

	Array< Pixel, 3 > distance3d;
	Array< Pixel, 3 > weights3d;
	Array< Pixel, 3 > implicit3d;
};


template< class Pixel >
nbfMinimalSurface< Pixel > :: nbfMinimalSurface( Array< Pixel, 3 > & weights )
{
	// set defaults
	this->weights3d.reference( weights );

	this->givenDistances = 0;
}

template< class Pixel >
void nbfMinimalSurface< Pixel > :: getDistance( Array< Pixel, 3 > & distance3d )
{
	// make sure array shapes are compatible
	distance3d.resize( this->weights3d.shape() );

#if 0 // compute full 3d distance (an not slice-wise)

	// 3d fast marching - compute distances from back face to rest of volume
	nbfFastMarching3D< Pixel > fastMarching3d( weights3d );
	vector< TinyVector< int, 3 > > positions;
	vector< Pixel > distances;

	// set to Alive all points in the back face with distance < \infty.
	for ( int i = 0; i < this->weights3d.rows(); i++ ){
		for ( int k = 0; k < this->weights3d.depth(); k++ ){
			if ( this->weights3d( i, this->weights3d.lbound(secondDim), k ) < numeric_limits<Pixel>::max() ){
//			if ( (*this->weights3d)( i, this->weights3d->lbound(secondDim), k ) < .5 ){
				TinyVector< int, 3 > last( i, this->weights3d.lbound(secondDim), k );
				positions.push_back( last );
				distances.push_back( 0 );
			}
		}
	}
	fastMarching3d.setAliveSet( positions, distances );

	fastMarching3d.execute( distance3d );
#if 0 // recompute 3 distances from a point IN the surface
	TinyVector< int, 2 > firstPoint;
	firstPoint = minIndex( distance3d( Range::all(), distance3d.ubound(secondDim), Range::all() ) );

	cout << "Start = " << firstPoint << endl;

	positions.clear();
	distances.clear();
	TinyVector< int, 3 > tmp( firstPoint(firstDim), distance3d.lbound(secondDim), firstPoint(secondDim) );
	positions.push_back(tmp);
	distances.push_back(0);
	fastMarching3d.setAliveSet( positions, distances );

	fastMarching3d.execute( distance3d );
#endif

#else 
	// compute distances on each slice (faster but biased).
	// Computing N NxN fast marchings is much faster than 1 NxNxN fast marching
	// but distances across slices may not be comparable and scaling is nedeed.

	nbfFastMarching2D< Pixel > fastMarching2d( weights3d( Range::all(), Range::all(), 0 ) );
	vector< TinyVector< int, 2 > > positions;
	vector< Pixel > distances;

	for ( int k = 0; k < weights3d.depth(); k++ ){
		positions.clear();
		distances.clear();
		for ( int j = 0; j < weights3d.cols(); j++ ){
			if ( weights3d( weights3d.lbound(firstDim), j, k ) < numeric_limits<Pixel>::max() ){
				TinyVector< int, 2 > last( weights3d.lbound(firstDim), j );
				positions.push_back( last );
				distances.push_back( 0 );
			}
		}
		fastMarching2d.setAliveSet( positions, distances );

		Array< Pixel, 2 > weights2d( weights3d( Range::all(), Range::all(), k ) );
		fastMarching2d.setWeights( weights2d );

		Array< Pixel, 2 > distance2d( distance3d( Range::all(), Range::all(), k ) );
		fastMarching2d.execute( distance2d );

		// sample scaling (this should be done on all slices)
		Array< Pixel, 2 > A( distance2d.shape() );
		A = where( distance2d < 1e9, distance2d, 0 );
		distance2d = distance2d / max(A);
	}

#endif

}

template< class Pixel >
void nbfMinimalSurface< Pixel > :: execute( Array< Pixel, 3 > & implicit3d,
											  Array< Pixel, 3 > & distance3d )
{
	// make sure array shapes are compatible
	implicit3d.resize( this->weights3d.shape() );

	// 2. Improved minimal surface determination
	//	- compute 3D distances from front slice
	//	- look for first arriving point to the back slice
	//	- find the circular geodesic (going through that point) from the
	//	  back to the front slice.
	//	- start computing geodesics (from back to front slice) to the right
	//	  restricted to a band around the path in the slice immediately to the left.
	//	- compute geodesics (from back to front slice) to the left
	//	  restricted to a band around the path in the slice immediately to the right.

	// look for first arriving point in the other side
	TinyVector< int, 2 > firstPoint;
	firstPoint = minIndex( distance3d( Range::all(), distance3d.ubound(secondDim), Range::all() ) );

	cout << "Start = " << firstPoint << endl;

	// size of the band (radius)
	int narrowBandMax = 3;
	int narrowBandMin = 1;

	Pixel bound = max( firstPoint(secondDim), distance3d.depth() - firstPoint(secondDim) );

	// get initial minimal path joining both sides and then go
	// to the right restricting each new path to a band around
	// the path in the previous slice.
	for ( int k = firstPoint(secondDim); k < distance3d.depth(); k++ ){

		Pixel x = k - firstPoint(secondDim);
		int narrowBand =  ceil( ( narrowBandMin - narrowBandMax ) / bound * x + narrowBandMax );
		// narrowBand = 5;

		// use accumulated 3D distances as image metric
		Array< Pixel, 2 > weights2d( distance3d( Range::all(), Range::all(), k ) );

		// if not in the first slice, enforce path to lay in the band
		if ( k > firstPoint(secondDim) ){
			weights2d = where( fabs( implicit3d( Range::all(), Range::all(), k - 1 ) ) <= narrowBand, weights2d, numeric_limits<Pixel>::max() );
			//for ( int h = 0; h < narrowBand; h++ ){
			//	Array< Pixel, 1 > sliceL( weights2d( Range::all(), h ) );
			//	Array< Pixel, 1 > implicitL( implicit3d( Range::all(), h, k - 1 ) );
			//	sliceL = where( fabs( implicitL ) <= h + 1, sliceL, numeric_limits<Pixel>::max() );

			//	Array< Pixel, 1 > sliceH( weights2d( Range::all(), weights2d.cols() - 1 - h ) );
			//	Array< Pixel, 1 > implicitH( implicit3d( Range::all(), weights2d.cols() - 1 - h, k - 1 ) );
			//	sliceH = where( fabs( implicitH ) <= h + 1, sliceH, numeric_limits<Pixel>::max() );
			//}

			//nbfMatlabWriter writer;
			//writer.setFileName( "data/distance.array" );
			//writer.write( weights2d );
			//exit(0);

		}


		// get circular path in the second dimension
		//nbfGeodesicPath< Pixel > geodesic2d( weights2d, secondDim );
		nbfGeodesicPath< Pixel > geodesic2d( weights2d );
		vector< TinyVector< int, 2 > > path;

		// store implicit geodesic representation in a 3D volume
		Array< Pixel,2 > ipath( implicit3d( Range::all(), Range::all(), k ) );

		// get implicit geodesic representation
		//geodesic2d.getForwardCircularPath( path, ipath );
	
		//cout << "k = " << k << endl;

		if ( path.size() == 1001 ){
			cout << "cannot find geodesic: = " << k << endl;
		}
	}

	// now go left of the initial path restricting the search to a band
	// around the minimal path in the slice inmediately right.
	for ( int k = firstPoint(secondDim) - 1; k >= 0; k-- ){

		Pixel x = firstPoint(secondDim) - k;
		int narrowBand =  ceil( ( narrowBandMin - narrowBandMax ) / bound * x + narrowBandMax );
		// narrowBand = 5;

		// use accumulated 3D distances as image metric
		//Array< Pixel, 2 > weights2d( (*weights3d)( Range::all(), Range::all(), k ) );
		Array< Pixel, 2 > weights2d( distance3d( Range::all(), Range::all(), k ) );

		// enforce path to lay in the band
		weights2d = where( fabs( implicit3d( Range::all(), Range::all(), k + 1 ) ) <= narrowBand, weights2d, numeric_limits<Pixel>::max() );

		//for ( int h = 0; h < narrowBand; h++ ){
		//	//Array< Pixel, 1 > sliceL( weights2d( Range::all(), h ) );
		//	//Array< Pixel, 1 > implicitL( implicit3d( Range::all(), h, k + 1 ) );
		//	//sliceL = where( fabs( implicitL ) <= h + 1, sliceL, numeric_limits<Pixel>::max() );

		//	//Array< Pixel, 1 > sliceH( weights2d( Range::all(), weights2d.cols() - 1 - h ) );
		//	//Array< Pixel, 1 > implicitH( implicit3d( Range::all(), weights2d.cols() - 1 - h, k + 1 ) );
		//	//sliceH = where( fabs( implicitH ) <= h + 1, sliceH, numeric_limits<Pixel>::max() );
		//}

		// get circular path in the first dimesion
		//nbfGeodesicPath< Pixel > geodesic2d( weights2d, secondDim );
		nbfGeodesicPath< Pixel > geodesic2d( weights2d );
		vector< TinyVector< int, 2 > > path;

		// store implicit geodesic representation in a 3D volume
		Array< Pixel,2 > ipath( implicit3d( Range::all(), Range::all(), k ) );

		// get implicit geodesic representation
		//geodesic2d.getForwardCircularPath( path, ipath );

		//cout << "k = " << k << endl;

		if ( path.size() == 1001 ){

			//nbfMatlabWriter writer;
			//writer.setFileName( "data/distance.array" );
			//writer.write( weights2d );
			cout << "cannot find geodesic: = " << k << endl;
		}
	}
}

template< class Pixel >
void nbfMinimalSurface< Pixel > :: executeOld( Array< Pixel, 3 > & implicit3d,
											  Array< Pixel, 3 > & distance3d )
{
	// make sure array shapes are compatible
	implicit3d.resize( this->weights3d->shape() );

	// use distances as new g ( g > 0 )
	distance3d = distance3d + 1e-5;

	// 2. Improved minimal surface determination
	//	- compute 3D distances from front slice
	//	- look for first arriving point to the back slice
	//	- find the circular geodesic (going through that point) from the
	//	  back to the front slice.
	//	- start computing geodesics (from back to front slice) to the right
	//	  restricted to a band around the path in the slice immediately to the left.
	//	- compute geodesics (from back to front slice) to the left
	//	  restricted to a band around the path in the slice immediately to the right.

	// look for first arriving point in the other side
	TinyVector< int, 2 > firstPoint;
	firstPoint = minIndex( distance3d( Range::all(), distance3d.ubound(secondDim), Range::all() ) );

	cout << "Start = " << firstPoint << endl;

	// size of the band (radius)
	int narrowBandMax = 3;
	int narrowBandMin = 1;

	Pixel bound = max( firstPoint(secondDim), distance3d.depth() - firstPoint(secondDim) );

	// get initial minimal path joining both sides and then go
	// to the right restricting each new path to a band around
	// the path in the previous slice.
	for ( int k = firstPoint(secondDim); k < distance3d.depth(); k++ ){

		Pixel x = k - firstPoint(secondDim);
		int narrowBand =  ceil( ( narrowBandMin - narrowBandMax ) / bound * x + narrowBandMax );
		//narrowBand = 3;

		// use accumulated 3D distances as image metric
		//Array< Pixel, 2 > weights2d( (*weights3d)( Range::all(), Range::all(), k ) );
		Array< Pixel, 2 > weights2d( distance3d( Range::all(), Range::all(), k ) );

		// if not in the first slice, enforce path to lay in the band
		if ( k > firstPoint(secondDim) ){
			weights2d = where( fabs( implicit3d( Range::all(), Range::all(), k - 1 ) ) < narrowBand, weights2d, numeric_limits<Pixel>::max() );
		}

		// get circular path in the first dimension
		nbfGeodesicPath< Pixel > geodesic2d( weights2d, secondDim );
		vector< TinyVector< int, 2 > > path;

		// store implicit geodesic representation in a 3D volume
		Array< Pixel,2 > ipath( implicit3d( Range::all(), Range::all(), k ) );

		// get implicit geodesic representation
		geodesic2d.getCircularPath( path, ipath );
		if ( path.size() == 1001 ){
			cout << "tancose = " << k << endl;
			//nbfMatlabWriter writer;
			//writer.setFileName( "data/distance.array" );
			//writer.write( weights2d );
			//exit(0);
		}
	}

	// now go left of the initial path restricting the search to a band
	// around the minimal path in the slice inmediately right.
	for ( int k = firstPoint(secondDim) - 1; k >= 0; k-- ){

		Pixel x = firstPoint(secondDim) - k;
		int narrowBand =  ceil( ( narrowBandMin - narrowBandMax ) / bound * x + narrowBandMax );
		//narrowBand = 3;

		// use accumulated 3D distances as image metric
		//Array< Pixel, 2 > weights2d( (*weights3d)( Range::all(), Range::all(), k ) );
		Array< Pixel, 2 > weights2d( distance3d( Range::all(), Range::all(), k ) );

		// enforce path to lay in the band
		weights2d = where( fabs( implicit3d( Range::all(), Range::all(), k + 1 ) ) < narrowBand, weights2d, numeric_limits<Pixel>::max() );

		// get circular path in the first dimesion
		nbfGeodesicPath< Pixel > geodesic2d( weights2d, secondDim );
		vector< TinyVector< int, 2 > > path;

		// store implicit geodesic representation in a 3D volume
		Array< Pixel,2 > ipath( implicit3d( Range::all(), Range::all(), k ) );

		// get implicit geodesic representation
		geodesic2d.getCircularPath( path, ipath );
		if ( path.size() == 1001 ){

			//nbfMatlabWriter writer;
			//writer.setFileName( "data/distance.array" );
			//writer.write( weights2d );
			cout << "tancose = " << k << endl;
		}
	}
}


template< class Pixel >
void nbfMinimalSurface< Pixel > :: executeUltimate( Array< Pixel, 3 > & implicit3d )
{
	nbfMatlabWriter writer;
	writer.setFileName("ipath");

	// make sure array shapes are compatible
	implicit3d.resize( this->weights3d.shape() );

	Array< Pixel, 3 > distance3d( implicit3d.shape() );

	// tmp array for computations
	Array< Pixel, 3 > temporal( implicit3d.shape() );

	// tolerance to consider areas where at least three geodesics are close apart
	Pixel band = 2;

	// for accurate back propagation
	// Pixel band = 3;
	distance3d = weights3d + 2;

	vector< TinyVector< int, 2 > > path;
	vector< TinyVector< int, 2 > > :: iterator iter;
	Pixel accumulator = 0;
		
	// compute theta circular geodesics - only for useful slices
	for ( int k = 0; k < weights3d.depth(); k++ ){

		Array< Pixel, 2 > weights2d( distance3d( Range::all(), Range::all(), k ) );

		Array< Pixel,2 > ipath( temporal( Range::all(), Range::all(), k ) );

		nbfGeodesicPath< Pixel > geodesic2d( weights2d );
		vector< TinyVector< Pixel, 2 > > path;
		geodesic2d.getCircularPath( weights2d, path );
		//geodesic2d.getImplicitOpenPath( path, ipath );
		geodesic2d.getImplicitPath( path, ipath, band + 1 );
	
		cout << "k = " << k << endl;
	}

	// 2.
	// get robust area (three consecutive geodesics are not far apart)

	writer.write(temporal);

	distance3d = 0;

	Array< Pixel, 2 > prevSlice( implicit3d.rows(), implicit3d.cols() );
	Array< Pixel, 2 > currSlice( prevSlice.shape() );
	Array< Pixel, 2 > nextSlice( prevSlice.shape() );

	for ( int z = 1; z < weights3d.depth() - 1; z++ ){
		prevSlice = fabs( temporal( Range::all(), Range::all(), z - 1 ) );
		currSlice = fabs( temporal( Range::all(), Range::all(), z ) );
		nextSlice = fabs( temporal( Range::all(), Range::all(), z + 1 ) );
		distance3d( Range::all(), Range::all(), z ) = where( ( prevSlice < band ) & ( currSlice < band ) & ( nextSlice < band ), 1, 0 );
	}

	// do first slice
	currSlice = fabs( temporal( Range::all(), Range::all(), 0 ) );
	nextSlice = fabs( temporal( Range::all(), Range::all(), 1 ) );
	distance3d( Range::all(), Range::all(), 0 ) = where( ( currSlice < band ) & ( nextSlice < band ), 1, 0 );

	// do last slice
	currSlice = fabs( temporal( Range::all(), Range::all(), weights3d.ubound(thirdDim) ) );
	prevSlice = fabs( temporal( Range::all(), Range::all(), weights3d.ubound(thirdDim) - 1 ) );
	distance3d( Range::all(), Range::all(), weights3d.ubound(thirdDim) ) = where( ( currSlice < band ) & ( prevSlice < band ), 1, 0 );

	// NOW 'distance3d' has the y-confidence area

	// compute geodesics on the phi direction

	// build new weights
	temporal = where( ( distance3d > 0 ) & ( weights3d < 1 ), weights3d / 10.0, weights3d );

	temporal += 2;

	// build common alive set
	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;
	for ( int i = 0; i < weights3d.rows(); i++ )
	{
		TinyVector< int, 2 > current( i, 0 );
		aliveP.push_back(current);
		aliveD.push_back(0);
	}			

	Array< Pixel, 1 > startSet( weights3d.cols() );
	Array< Pixel, 1 > endSet( weights3d.cols() );

	for ( int j = 0; j < distance3d.cols(); j++ ){

		Array< Pixel, 2 > weights2d( temporal( Range::all(), j, Range::all() ) );
		Array< Pixel, 2 > ipath( implicit3d( Range::all(), j, Range::all() ) );

		nbfGeodesicPath< Pixel > geodesic2d( weights2d );
		vector< TinyVector< Pixel, 2 > > path;
		
		nbfFastMarching2D8< Pixel > fm( weights2d );
		fm.setAliveSet( aliveP, aliveD );
		fm.setStopBorder( secondDim );

		Array< Pixel, 2 > distance2d( weights2d.shape() );

		TinyVector< int, 2 > firstArrival = fm.execute( distance2d );
		//geodesic2d.getForwardLinearPath( distance2d, firstArrival, path );
		geodesic2d.getSimplePath( distance2d, firstArrival, path );
		geodesic2d.getImplicitPath( path, ipath, band + 1 );
	
		cout << "j = " << j << endl;
	}

	writer.write(implicit3d);

	// get robust area into temporal
	temporal = 0;
	prevSlice.resize( implicit3d.rows(), implicit3d.cols() );
	currSlice.resize( implicit3d.rows(), implicit3d.cols() );
	nextSlice.resize( implicit3d.rows(), implicit3d.cols() );
	for ( int j = 1; j < implicit3d.cols() - 1; j++ ){
		prevSlice = fabs( implicit3d( Range::all(), j - 1, Range::all() ) );
		currSlice = fabs( implicit3d( Range::all(), j, Range::all() ) );
		nextSlice = fabs( implicit3d( Range::all(), j + 1, Range::all() ) );
		temporal( Range::all(), j, Range::all() ) = where( ( prevSlice < band ) & ( currSlice < band ) & ( nextSlice < band ), 1, 0 );
	}

	currSlice = fabs( implicit3d( Range::all(), implicit3d.ubound(secondDim), Range::all() ) );
	temporal( Range::all(), implicit3d.ubound(secondDim), Range::all() ) = where( currSlice < band, 1, 0 );
	
	currSlice = fabs( implicit3d( Range::all(), implicit3d.lbound(secondDim), Range::all() ) );
	temporal( Range::all(), implicit3d.lbound(secondDim), Range::all() ) = where( currSlice < band, 1, 0 );

	// set to intersection of robust areas
	distance3d = where( ( distance3d > 0 ) & ( temporal > 0 ), 1, 0 );

	writer.write(distance3d);

	implicit3d = distance3d;
	return;

	writer.write(implicit3d);
}


template< class Pixel >
void nbfMinimalSurface< Pixel > :: executeNew( Array< Pixel, 3 > & implicit3d,
											   Array< Pixel, 3 > & distance3d )
{
	nbfMatlabWriter writer;
	writer.setFileName("ipath");

	// make sure array shapes are compatible
	implicit3d.resize( this->weights3d.shape() );
	distance3d.resize( this->weights3d.shape() );

	// tmp array for computations
	Array< Pixel, 3 > temporal( implicit3d.shape() );

	// compute first useful slice
	int firstThetaSlice = 0;
	while ( true ){
		Array< Pixel, 1 > line( weights3d( Range::all(), weights3d.lbound(secondDim), firstThetaSlice ) );
		if ( min( line ) < 2 ){
			break;
		}
		firstThetaSlice += 1;
	}

	// compute last useful slice
	int lastThetaSlice = weights3d.ubound(thirdDim);
	while ( true ){
		Array< Pixel, 1 > line( weights3d( Range::all(), weights3d.lbound(secondDim), lastThetaSlice ) );
		if ( min( line ) < 2 ){
			break;
		}
		lastThetaSlice -= 1;
	}

	// tolerance to consider areas where at least three geodesics are close apart
	Pixel band = 3;

	// for accurate back propagation
	// Pixel band = 3;

	vector< TinyVector< int, 2 > > path;
	vector< TinyVector< int, 2 > > :: iterator iter;
	Pixel accumulator = 0;
		
	//// compute r_0 (average of geodesic coordinates in the first valid slice)
	//Array< Pixel, 2 > weightsr_0( weights3d( Range::all(), Range::all(), firstThetaSlice ) );
	//nbfGeodesicPath< Pixel > geodesic2d( weightsr_0 );
	//geodesic2d.getCircularPath( weightsr_0, path );
	//for( iter = path.begin(); iter != path.end(); ++iter ){
	//	accumulator += (*iter)[firstDim];
	//}
	//int r_0 = floor( accumulator / path.size() );

	//// compute r_1 (average of geodesic coordinates in the last valid slice)
	//accumulator = 0;
	//Array< Pixel, 2 > weightsr_1( weights3d( Range::all(), Range::all(), lastThetaSlice ) );
	//geodesic2d.getCircularPath( weightsr_1, path );
	//for( iter = path.begin(); iter != path.end(); ++iter ){
	//	accumulator += (*iter)[firstDim];
	//}
	//int r_1 = floor( accumulator / path.size() );

	//TinyVector< int, 2 > startPoint1( r_0, weights3d.lbound(thirdDim) );
	//TinyVector< int, 2 > endPoint1( r_1, weights3d.ubound(thirdDim) );

	//cout << "start = " << startPoint1 << endl;
	//cout << "end = " << endPoint1 << endl;

	//// compute boundary curve C0
	//Array< Pixel, 2 > weightsC0( weights3d.rows(), weights3d.depth() );
	//weightsC0 = weights3d( Range::all(), weights3d.ubound(secondDim), Range::all() );

	//vector< int > C0;
	//vector< Pixel > C0_d;
	//for ( int k = firstThetaSlice; k <= lastThetaSlice; k++ ){
	//	Array< Pixel, 2 > weights2d( weights3d( Range::all(), Range::all(), k ) );
	//	nbfGeodesicPath< Pixel > geodesic2d( weights2d );
	//	TinyVector< int, 2 > first = geodesic2d.getFirstPointOnRightSide( weights2d );
	//	C0.push_back( first[firstDim] );
	//	C0_d.push_back( weightsC0( first[firstDim], k ) );
	//}

	//for ( int i = 1; i < C0.size() - 1; i++ ){
	//	if ( ( C0_d[i-1] < 1 ) & ( C0_d[i] < 1 ) & ( C0_d[i+1] < 1 ) ){
	//		if ( ( fabs( C0[i] - C0[i-1] + 0.0 ) < 3 ) & ( fabs( C0[i] - C0[i+1] + 0.0 ) < 3 ) ){
	//			weightsC0( Range( C0[i] - 1, C0[i] + 1 ), i + firstThetaSlice ) /= 10.0;
	//		}
	//	}
	//}
	//writer.write(weightsC0);

	//nbfFastMarching2D< Pixel > fmC0( weightsC0 );
	//fmC0.setAliveSet( startPoint, 0 );
	//fmC0.setStopPoint( endPoint );
	//Array< Pixel, 2 > distanceC0( weightsC0.shape() );
	//fmC0.execute( distanceC0 );
	////writer.write(distanceC0);
	//nbfGeodesicPath< Pixel > geodesicC0( distanceC0 );
	//geodesicC0.getSimplePath( distanceC0, endPoint, path );
	//Array< Pixel, 2 > ipathC0( weightsC0.shape() );
	//geodesicC0.getImplicitPath( path, ipathC0 );
	//writer.write(ipathC0);

	// compute theta circular geodesics - only for useful slices
	 for ( int k = firstThetaSlice; k <= lastThetaSlice; k++ ){

		Array< Pixel, 2 > weights2d( weights3d( Range::all(), Range::all(), k ) );

		Array< Pixel,2 > ipath( temporal( Range::all(), Range::all(), k ) );

		nbfGeodesicPath< Pixel > geodesic2d( weights2d );
		vector< TinyVector< Pixel, 2 > > path;
		geodesic2d.getCircularPath( weights2d, path );
		//geodesic2d.getImplicitOpenPath( path, ipath );
		geodesic2d.getImplicitPath( path, ipath );
	
		writer.write(weights2d);
		writer.write(ipath);

		cout << "k = " << k << endl;
	}

	// 2.
	// get robust area (three consecutive geodesics are not far apart)

	writer.write(temporal);

	distance3d = 0;

	Array< Pixel, 2 > prevSlice( implicit3d.rows(), implicit3d.cols() );
	Array< Pixel, 2 > currSlice( prevSlice.shape() );
	Array< Pixel, 2 > nextSlice( prevSlice.shape() );

	for ( int z = firstThetaSlice + 1; z < lastThetaSlice; z++ ){
		prevSlice = fabs( temporal( Range::all(), Range::all(), z - 1 ) );
		currSlice = fabs( temporal( Range::all(), Range::all(), z ) );
		nextSlice = fabs( temporal( Range::all(), Range::all(), z + 1 ) );
		distance3d( Range::all(), Range::all(), z ) = where( ( prevSlice < band ) & ( currSlice < band ) & ( nextSlice < band ), 1, 0 );
	}

	// do first slice
	currSlice = fabs( temporal( Range::all(), Range::all(), firstThetaSlice ) );
	nextSlice = fabs( temporal( Range::all(), Range::all(), firstThetaSlice + 1 ) );
	distance3d( Range::all(), Range::all(), firstThetaSlice ) = where( ( currSlice < band ) & ( nextSlice < band ), 1, 0 );

	// do last slice
	currSlice = fabs( temporal( Range::all(), Range::all(), lastThetaSlice ) );
	prevSlice = fabs( temporal( Range::all(), Range::all(), lastThetaSlice - 1 ) );
	distance3d( Range::all(), Range::all(), lastThetaSlice ) = where( ( currSlice < band ) & ( prevSlice < band ), 1, 0 );

	// NOW 'distance3d' has the y-confidence area

	// compute geodesics on the phi direction

	// build new weights
	temporal = where( ( distance3d > 0 ) & ( weights3d < 1 ), weights3d / 10.0, weights3d );

	temporal += 1;

	// find r_0 and r_1
	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;
	for ( int i = 0; i < weights3d.rows(); i++ )
	{
		TinyVector< int, 2 > current( i, 0 );
		aliveP.push_back(current);
		aliveD.push_back(0);
	}			

	Array< Pixel, 1 > startSet( weights3d.cols() );
	Array< Pixel, 1 > endSet( weights3d.cols() );

	for ( int j = 0; j < distance3d.cols(); j++ ){

		Array< Pixel, 2 > weights2d( temporal( Range::all(), j, Range( firstThetaSlice, lastThetaSlice ) ) );
		Array< Pixel, 2 > ipath( implicit3d( Range::all(), j, Range( firstThetaSlice, lastThetaSlice ) ) );

		nbfGeodesicPath< Pixel > geodesic2d( weights2d );
		vector< TinyVector< Pixel, 2 > > path;
		
		nbfFastMarching2D8< Pixel > fm( weights2d );
		fm.setAliveSet( aliveP, aliveD );
		fm.setStopBorder( secondDim );

		Array< Pixel, 2 > distance2d( weights2d.shape() );

		TinyVector< int, 2 > firstArrival = fm.execute( distance2d );
		//geodesic2d.getForwardLinearPath( distance2d, firstArrival, path );
		geodesic2d.getSimplePath( distance2d, firstArrival, path );
		geodesic2d.getImplicitPath( path, ipath );
	
		cout << "j = " << j << endl;
	}

	writer.write(implicit3d);

	TinyVector< int, 2 > startPoint( sum(startSet)/startSet.size(), weights3d.lbound(thirdDim) );
	TinyVector< int, 2 > endPoint( sum(endSet)/endSet.size(), weights3d.ubound(thirdDim) );

	cout << "start = " << startPoint << endl;
	cout << "end = " << endPoint << endl;

	//// find geodesics in the phi direction (fix start and ending point)
	//for ( int j = 0; j < distance3d.cols(); j++ ){

	//	Array< Pixel, 2 > weights2d( temporal( Range::all(), j, Range( firstThetaSlice, lastThetaSlice ) ) );
	//	Array< Pixel,2 > ipath( implicit3d( Range::all(), j, Range( firstThetaSlice, lastThetaSlice ) ) );

	//	nbfGeodesicPath< Pixel > geodesic2d( weights2d );
	//	vector< TinyVector< Pixel, 2 > > path;
	//	
	//	nbfFastMarching2D8< Pixel > fm( weights2d );
	//	fm.setAliveSet( startPoint, 0 );
	//	fm.setStopPoint( endPoint );
	//	fm.execute( distance2d );
	//	geodesic2d.getForwardLinearPath( distance2d, endPoint, path );
	//	//geodesic2d.getSimplePath( distance2d, endPoint, path );
	//	geodesic2d.getImplicitOpenPath( path, ipath );
	//
	//	cout << "j = " << j << endl;
	//}

	// get robust area into temporal
	temporal = 0;
	prevSlice.resize( implicit3d.rows(), implicit3d.cols() );
	currSlice.resize( implicit3d.rows(), implicit3d.cols() );
	nextSlice.resize( implicit3d.rows(), implicit3d.cols() );
	for ( int j = 1; j < implicit3d.cols() - 1; j++ ){
		prevSlice = fabs( implicit3d( Range::all(), j - 1, Range::all() ) );
		currSlice = fabs( implicit3d( Range::all(), j, Range::all() ) );
		nextSlice = fabs( implicit3d( Range::all(), j + 1, Range::all() ) );
		temporal( Range::all(), j, Range::all() ) = where( ( prevSlice < band ) & ( currSlice < band ) & ( nextSlice < band ), 1, 0 );
	}

	currSlice = fabs( implicit3d( Range::all(), implicit3d.ubound(secondDim), Range::all() ) );
	temporal( Range::all(), implicit3d.ubound(secondDim), Range::all() ) = where( currSlice < band, 1, 0 );
	
	currSlice = fabs( implicit3d( Range::all(), implicit3d.lbound(secondDim), Range::all() ) );
	temporal( Range::all(), implicit3d.lbound(secondDim), Range::all() ) = where( currSlice < band, 1, 0 );

	// set to intersection of robust areas
	distance3d = where( ( distance3d > 0 ) & ( temporal > 0 ), 1, 0 );

	// set side lines to "good"
	distance3d( startPoint[firstDim], Range::all(), startPoint[secondDim] ) = 1;
	distance3d( endPoint[firstDim], Range::all(), endPoint[secondDim] ) = 1;

	writer.write(distance3d);

	implicit3d = distance3d;
	return;

	// extend mask to full x-extent
	Array< Pixel, 2 > splash( distance3d.cols(), distance3d.depth() );
	Array< Pixel, 2 > :: iterator iter2 = splash.begin();
	while ( iter2 != splash.end() ){
		Array< Pixel, 1 > line( distance3d( Range::all(), iter2.position()[firstDim], iter2.position()[secondDim] ) );
		Pixel suma = sum(line);
		if ( suma > 0 ){
			line = 1;
		}
		iter2++;
	}

	//for ( int i = 0; i < 100; i++ ){
	//	applyStencil( surfaceSmooth3D(), temporal, implicit3d );
	//	implicit3d = where( distance3d < 1, temporal, implicit3d );
	//	// fix borders
	//	implicit3d( Range::all(), implicit3d.lbound(secondDim), Range::all() ) = implicit3d( Range::all(), implicit3d.lbound(secondDim) + 1, Range::all() );
	//	implicit3d( Range::all(), implicit3d.ubound(secondDim), Range::all() ) = implicit3d( Range::all(), implicit3d.lbound(secondDim), Range::all() );
	//}

	Range I(0,implicit3d.ubound(firstDim));
	Range J(1,implicit3d.ubound(secondDim)-1);
	Range K(1,implicit3d.ubound(thirdDim)-1);

	for ( int i = 0; i < 100; i++ ){
		temporal(I,J,K) = 
			( implicit3d(I,J,K-1) + implicit3d(I,J,K) + implicit3d(I,J,K+1) + 
			implicit3d(I,J-1,K-1) + implicit3d(I,J-1,K) + implicit3d(I,J-1,K+1) + 
			implicit3d(I,J+1,K-1) + implicit3d(I,J+1,K) + implicit3d(I,J+1,K+1) ) / 9.0; 
		temporal(I,0,K) = 
			( implicit3d(I,0,K-1) + implicit3d(I,0,K) + implicit3d(I,0,K+1) + 
			implicit3d(I,implicit3d.ubound(secondDim),K-1) + implicit3d(I,implicit3d.ubound(secondDim),K) + implicit3d(I,implicit3d.ubound(secondDim),K+1) + 
			implicit3d(I,1,K-1) + implicit3d(I,1,K) + implicit3d(I,1,K+1) ) / 9.0; 
		temporal(I,temporal.ubound(secondDim),K) = 
			( implicit3d(I,temporal.ubound(secondDim),K-1) + implicit3d(I,temporal.ubound(secondDim),K) + implicit3d(I,temporal.ubound(secondDim),K+1) + 
			implicit3d(I,temporal.ubound(secondDim)-1,K-1) + implicit3d(I,temporal.ubound(secondDim)-1,K) + implicit3d(I,temporal.ubound(secondDim)-1,K+1) + 
			implicit3d(I,temporal.lbound(secondDim),K-1) + implicit3d(I,temporal.lbound(secondDim),K) + implicit3d(I,temporal.lbound(secondDim),K+1) ) / 9.0; 
		implicit3d = where( distance3d < 1, temporal, implicit3d );
	}

	// force to a plane
	//implicit3d( Range(fromStart,implicit3d.ubound(firstDim)/4),  Range::all(), Range::all() ) = -1;
	//implicit3d( Range(implicit3d.ubound(firstDim)/4,toEnd), Range::all(), Range::all() )= 1;

	//for ( int i = 0; i < 5; i++ ){
	//	applyStencil( surfaceSmooth3D(), temporal, implicit3d );
	//	implicit3d = temporal;
	//	bsForV.refresh();
	//}

	writer.write(implicit3d);

#if 0
	// initial surface inteprolation
	distance3d = implicit3d;
	for ( int j = 0; j < phi3d.cols(); j++ ){
		this->fillSlice( phi3d( Range::all(), j, Range::all() ), implicitX( Range::all(), j, Range::all() ),distanceX( Range::all(), j, Range::all() ) );
	}

	// set start and ending lines
	phi3d( Range::all(), phi3d.lbound(secondDim), Range::all() ) = where( fabs( implicitX( Range::all(), phi3d.lbound(secondDim), Range::all() ) ) < 2, 1, phi3d( Range::all(), phi3d.lbound(secondDim), Range::all() ) );
	phi3d( Range::all(), phi3d.ubound(secondDim), Range::all() ) = where( fabs( implicitX( Range::all(), phi3d.ubound(secondDim), Range::all() ) ) < 2, 1, phi3d( Range::all(), phi3d.ubound(secondDim), Range::all() ) );

	writer.write(phi3d);

	implicitY = implicitX;
	for ( int k = 0; k < phi3d.depth(); k++ ){
		this->fillSlice( phi3d( Range::all(), Range::all(), k ), implicitY( Range::all(), Range::all(), k ), distanceY( Range::all(), Range::all(), k ) );
	}


	for ( int j = 2; j < implicit3d.cols() - 2; j++ ){
		prevPhi = implicitX( Range::all(), j - 1, Range::all() );
		currPhi = implicitX( Range::all(), j, Range::all() );
		nextPhi = implicitX( Range::all(), j + 1, Range::all() );
		distance3d( Range::all(), j, Range::all() ) = prevPhi + currPhi + nextPhi;
		prevPhi = implicitX( Range::all(), j - 2, Range::all() );
		nextPhi = implicitX( Range::all(), j + 2, Range::all() );
		distance3d( Range::all(), j, Range::all() ) = ( distance3d( Range::all(), j, Range::all() ) + prevPhi + nextPhi ) / 5.0;
	}
	implicitX = distance3d;

	for ( int z = 2; z < implicit3d.depth() - 2; z++ ){
		prevSlice = implicitY( Range::all(), Range::all(), z - 1 );
		currSlice = implicitY( Range::all(), Range::all(), z );
		nextSlice = implicitY( Range::all(), Range::all(), z + 1 );
		distance3d( Range::all(), Range::all(), z ) = prevSlice + currSlice + nextSlice;
		prevSlice = implicitY( Range::all(), Range::all(), z - 2 );
		nextSlice = implicitY( Range::all(), Range::all(), z + 2 );
		distance3d( Range::all(), Range::all(), z ) = ( distance3d( Range::all(), Range::all(), z ) + prevSlice + nextSlice ) / 5.0;
	}
	implicitY( Range::all(), Range::all(), Range(2,implicit3d.depth()-3) ) = distance3d( Range::all(), Range::all(), Range(2,implicit3d.depth()-3) );

	//writer.write(distance3d);

//	implicit3d = where( distanceX + distanceY > 0, (distanceY) / ( (distanceX) + (distanceY) ) * implicitX + (distanceX) / ( (distanceX) + (distanceY) ) * implicitY, ( implicitX + implicitY ) / 2.0 );
	implicit3d = where( distanceX + distanceY > 0, pow2(distanceY) / ( pow2(distanceX) + pow2(distanceY) ) * implicitX + pow2(distanceX) / ( pow2(distanceX) + pow2(distanceY) ) * implicitY, ( implicitX + implicitY ) / 2.0 );

#endif


#if 0
	// extend area one pixel below and above in the Z axis to force 4-connectivity
	Array< Pixel, 3 > middle( phi3d( Range(1,phi3d.rows()-2), Range::all(), Range::all() ) );
	Array< Pixel, 3 > lower( phi3d( Range(0,phi3d.rows()-3), Range::all(), Range::all() ) );
	Array< Pixel, 3 > upper( phi3d( Range(2,phi3d.rows()-1), Range::all(), Range::all() ) );
	
	Array< Pixel, 3 > target( distance3d( Range(1,phi3d.rows()-2), Range::all(), Range::all() ) );
	
	target = where( ( middle > 0 ) & ( upper < 1 ) & ( lower < 1 ), 1, 0 );
	upper = where( target > 0, 1, upper );
	lower = where( target > 0, 1, lower );

	writer.write(phi3d);

	// set to infinity all weights that are not robust
	Array< Pixel, 2 > splash( distance3d.cols(), distance3d.depth() );
	Array< Pixel, 2 > :: iterator iter2 = splash.begin();
	while ( iter2 != splash.end() ){
		Array< Pixel, 1 > line( phi3d( Range::all(), iter2.position()[firstDim], iter2.position()[secondDim] ) );
		Pixel suma = sum(line);
		if ( suma > 0 ){
			line = where( line < 1, numeric_limits<Pixel>::max(), line );
		}
		else{
			line = 1;
		}
		iter2++;
	}
	phi3d = where( ( weights3d > 3 ) & ( weights3d < numeric_limits<Pixel>::max() ), 1, phi3d );

	//writer.write(phi3d);

	//// theta direction (periodic geodesics)
	//cout << sidePath.size() << " vs. "<< distance3d.depth() << endl;
	//int index = sidePath.size() - 1;
	//distance2d.resize( distance3d.rows(), distance3d.cols() );
	//for ( int k = 0; k < distance3d.depth(); k++ ){

	//	Array< Pixel, 2 > weights2d( phi3d( Range::all(), Range::all(), k ) );
	//	Array< Pixel, 2 > ipath( zone3d( Range::all(), Range::all(), k ) );

	//	nbfGeodesicPath< Pixel > geodesic2d( weights2d );
	//	vector< TinyVector< int, 2 > > path;

	//	writer.write(weights2d);
	//	
	//	geodesic2d.getCircularPath( weights2d, path );

	//	geodesic2d.getImplicitOpenPath( path, ipath );
	//	writer.write(ipath);

	//	cout << "k = " << k << endl;
	//}

	writer.write(zone3d);

	// phi direction (fix start and ending point)
	for ( int j = 0; j < distance3d.cols(); j++ ){

		Array< Pixel, 2 > weights2d( phi3d( Range::all(), j, Range::all() ) );
		Array< Pixel,2 > ipath( implicit3d( Range::all(), j, Range::all() ) );

		nbfGeodesicPath< Pixel > geodesic2d( weights2d );
		vector< TinyVector< int, 2 > > path;
		
		nbfFastMarching2D8< Pixel > fm( weights2d );
		fm.setAliveSet( startPoint, 0 );
		fm.setStopPoint( endPoint );
		fm.execute( distance2d );
		
		writer.write(distance2d);

		geodesic2d.getForwardLinearPath( distance2d, endPoint, path );
		geodesic2d.getImplicitOpenPath( path, ipath );
	
		cout << "j = " << j << endl;
	}

	writer.write(implicit3d);
	
	implicit3d = ( implicit3d + zone3d ) / 2.0;
	writer.write(implicit3d);
#endif
}


template< class Pixel >
void nbfMinimalSurface< Pixel > :: fillSlice( Array< Pixel, 2 > & regions,
											  Array< Pixel, 2 > & implicit,
											  Array< Pixel, 2 > & distance )
{
	nbfMatlabWriter writer;
	writer.setFileName("ipath");

	distance = 0;
	//writer.write(regions);

	Array< Pixel, 1 > binary( regions.cols() );
	//firstIndex indexI;
	//secondIndex indexJ;
	for ( int i = 0; i < regions.cols(); i++ ){
		binary(i) = sum( regions( Range::all(), i ) );
	}
	writer.write(binary);
	Pixel colStart = 0;
	// find first
	while ( binary( colStart ) > 0 ){
		colStart++;
	}
	colStart--;
	Pixel colEnd = colStart + 1;
	while ( colStart < ( regions.cols() - 1 ) ){
		// find end
		while ( binary( colEnd ) < 1 ){
			colEnd++;
		}
		
		//cout << colStart << ", " << colEnd << endl;

		// interpolate
		for ( int k = colStart + 1; k < colEnd; k++ ){
			implicit( Range::all(), k ) = implicit( Range::all(), colEnd ) * ( k - colStart ) / ( colEnd - colStart ) + implicit( Range::all(), colStart ) * ( 1.0 - ( k - colStart ) / ( colEnd - colStart ) );
			distance( Range::all(), k ) = min( k - colStart, colEnd - k );
		}

		// go to next
		colStart = colEnd;
		while ( ( binary( colStart ) > 0 ) & ( colStart < regions.cols() ) ){
			colStart++;
		}
		colStart--;
		colEnd = colStart + 1;
	}

	//writer.write(implicit);

}

BZ_DECLARE_STENCIL2(surfaceSmooth3D,A,B)
  A = ( B(0,0,0) + B(0,0,1) + B(0,0,-1) + B(0,1,0) + B(0,1,1) + B(0,1,-1) + B(0,-1,0) + B(0,-1,1) + B(0,-1,-1) ) / 9.0;
BZ_END_STENCIL


template< class Pixel >
void nbfMinimalSurface< Pixel > :: search( 
	Array< Pixel, 3 > & inputX,
	Array< Pixel, 3 > & inputY,
	Array< Pixel, 3 > & inputZ,
	TinyVector< int, 3 > & center,
	Array< Pixel, 3 > & output )
{
	// store temporal confidence values
	Array< bool, 3 > confidenceX( inputX.shape() );
	Array< bool, 3 > confidenceY( inputY.shape() );
	Array< bool, 3 > confidenceZ( inputZ.shape() );

	implicit3d.resize( this->weights3d.shape() );

	int band = 2;

	// compute circular geodesics in Z slices
	this->computeGeodesics( inputZ, thirdDim, band, center, implicit3d );

	nbfMatlabWriter writer;
	writer.setFileName("ipath");
	writer.write(implicit3d);

	// compute confidence from Z slices
	this->computeConfidence( implicit3d, thirdDim, band, confidenceZ );

	writer.write(confidenceZ);

	// update X and Y weights
	inputX = where( confidenceZ & ( inputX < 1 ), inputX / 10.0, inputX );
	inputY = where( confidenceZ & ( inputY < 1 ), inputY / 10.0, inputY );

	writer.write(inputX);

	// compute circular geodesics in X slices
	this->computeGeodesics( inputX, firstDim, band, center, implicit3d );

	writer.write(implicit3d);

	// compute confidence from X slices
	this->computeConfidence( implicit3d, firstDim, band, confidenceX );

	writer.write(confidenceX);

	// compute circular geodesics in Y slices
	this->computeGeodesics( inputY, secondDim, band, center, implicit3d );

	writer.write(implicit3d);

	// compute confidence from Y slices
	this->computeConfidence( implicit3d, secondDim, band, confidenceY );

	writer.write(confidenceY);

	output = where( confidenceZ & ( confidenceX | confidenceY ), 1, 0 );

	output = where( confidenceZ & confidenceX & confidenceY, 1, 0 );
}

// IP-SI paper version
// vuelta a theta/phi/rho coordinates with the restricted domain
template< class Pixel >
void nbfMinimalSurface< Pixel > :: search( 
	Array< Pixel, 3 > & input,
	Array< Pixel, 3 > & output )
{
	nbfImageWriter writer;
	writer.setFileName("ipath.vtk");

	vector< Pixel > positions;
	vector< Pixel > midPositions;

	vector< TinyVector< Pixel, 2 > > path;

	Array< Pixel, 2 > slice( input.rows(), input.depth() * 2 );
	Array< Pixel, 2 > part( input.rows(), input.depth() );
	
	// First, robustly compute boundary condition
	for ( int k = 0; k < 3; k++ ){
		// make image periodic in the \phi direction to compute closed geodesics in \theta_i
		part = input( Range::all(), k, Range::all() );
		part.reverseSelf(secondDim);
		slice( Range::all(), Range(fromStart,input.depth()-1) ) = part;
		slice( Range::all(), Range(input.depth(),toEnd) )		= input( Range::all(), k + input.cols() / 2, Range::all() );
		nbfGeodesicPath< Pixel > geodesic( slice );
		geodesic.getCircularPath( slice, path );
		
		//writer.write(slice);
		//Array< float, 2 > ipath( slice.shape() );
		//geodesic.getImplicitPath(path,ipath);
		//writer.write(ipath);

		positions.push_back( path[0][firstDim] );
		vector< TinyVector< Pixel, 2 > > :: iterator iter = path.begin();
		while ( (*iter)[secondDim] > input.depth() ){
			++iter;
		}
		midPositions.push_back( (*iter)[firstDim] );
	}

	for ( int k = input.cols() / 2 - 3; k < input.cols() / 2; k++ ){
		//cout << k << endl;
		part = input( Range::all(), k, Range::all() );
		part.reverseSelf(secondDim);
		slice( Range::all(), Range(fromStart,input.depth()-1) ) = part;
		slice( Range::all(), Range(input.depth(),toEnd) )		= input( Range::all(), k + input.cols() / 2, Range::all() );
		nbfGeodesicPath< Pixel > geodesic( slice );
		geodesic.getCircularPath( slice, path );
		positions.push_back( path[0][firstDim] );
		vector< TinyVector< Pixel, 2 > > :: iterator iter = path.begin();
		while ( (*iter)[secondDim] > input.depth() ){
			++iter;
		}
		midPositions.push_back( (*iter)[firstDim] );
	}

	vector< Pixel > :: iterator first = positions.begin();       
	vector< Pixel > :: iterator last  = positions.end();       
	sort( first, last );
	first = midPositions.begin();
	last = midPositions.end();
	sort( first, last );

	//cout << positions[ floor( positions.size() / 2.0 ) ] << endl;
	//cout << midPositions[ floor( midPositions.size() / 2.0 ) ] << endl;
	
	Pixel m = positions[ floor( positions.size() / 2.0 ) ];
	if ( m - floor(m) < .5 ){
		m = floor(m);
	}
	else{
		m = ceil(m);
	}
	TinyVector< int, 2 > endP( floor(m), input.ubound(thirdDim) );
	//cout << endP << endl;

	m = midPositions[ floor( midPositions.size() / 2.0 ) ];
	if ( m - floor(m) < .5 ){
		m = floor(m);
	}
	else{
		m = ceil(m);
	}
	TinyVector< int, 2 > startP( floor(m), input.lbound(thirdDim) );
	//cout << startP << endl;

	output.resize( input.shape() );

	Pixel lowBound = min(input);
	Pixel  upBound = max(input);

	//cout << lowBound << endl;
	//cout <<  upBound << endl;

	TinyVector< int, 2 > start, end;

#if 0
	for ( int k = 0; k < 45; k++ ){

		cout << k << endl;

		Array< Pixel, 2 > slice, ipath;
	
		// do all four geodesics in \theta

		start = startP;
		end	  = endP;

		// from \theta = 0
		slice.reference( input(  Range::all(), k, Range::all() ) );
		ipath.reference( output( Range::all(), k, Range::all() ) );
		this->computeGeodesic( slice, start, end, ipath );
		//writer.write(slice);
		//writer.write(ipath);
		slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), upBound );
		
		// from \theta = \pi forward
		slice.reference( input(  Range::all(), 90 + k, Range::all() ) );
		ipath.reference( output( Range::all(), 90 + k, Range::all() ) );
		this->computeGeodesic( slice, start, end, ipath );
		//writer.write(slice);
		//writer.write(ipath);
		if ( k == 0 ){
			slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), numeric_limits<Pixel>::max() );
		}
		else{
			slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), upBound );
		}

		// from \theta = \pi backward
		slice.reference( input(  Range::all(), 89 - k, Range::all() ) );
		ipath.reference( output( Range::all(), 89 - k, Range::all() ) );
		this->computeGeodesic( slice, start, end, ipath );
		//writer.write(slice);
		//writer.write(ipath);
		slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), upBound );

		// from \theta = 2 \pi
		slice.reference( input(  Range::all(), 179 - k, Range::all() ) );
		ipath.reference( output( Range::all(), 179 - k, Range::all() ) );
		this->computeGeodesic( slice, start, end, ipath );
		//writer.write(slice);
		//writer.write(ipath);
		slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), upBound );

		// do all four geodesics in \phi

		if ( fmod(k,2.0) == 0 ){

			TinyVector< int, 1 > p;

			// from \phi = \phi / 4 forward
			slice.reference( input(  Range::all(), Range::all(), 22 + k / 2 ) );
			ipath.reference( output( Range::all(), Range::all(), 22 + k / 2 ) );

			//writer.write( output( Range::all(), 0, 22 + k / 2 ) );

			p = minIndex( fabs( output( Range::all(), 0, 22 + k / 2 ) ) );
			start = TinyVector< int, 2 >(     p[firstDim], input.lbound(secondDim) );
			end	  = TinyVector< int, 2 >( start[firstDim], input.ubound(secondDim) );

			//cout << start << endl;
			//cout << end << endl;
			//writer.write(slice);

			this->computeGeodesic( slice, start, end, ipath );
		//writer.write(slice);
		//writer.write(ipath);
			//writer.write(ipath);
			slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), upBound );
			//writer.write(slice);

			// from \phi = \phi / 4 backward
			slice.reference( input(  Range::all(), Range::all(), 22 - k / 2 ) );
			ipath.reference( output( Range::all(), Range::all(), 22 - k / 2 ) );

			p = minIndex( fabs( output( Range::all(), 0, 22 - k / 2 ) ) );
			start = TinyVector< int, 2 >(     p[firstDim], input.lbound(secondDim) );
			end	  = TinyVector< int, 2 >( start[firstDim], input.ubound(secondDim) );

			this->computeGeodesic( slice, start, end, ipath );
		//writer.write(slice);
		//writer.write(ipath);
			slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), upBound );

			// from \phi = 3 \phi / 4 forward
			slice.reference( input(  Range::all(), Range::all(), 67 + k / 2 ) );
			ipath.reference( output( Range::all(), Range::all(), 67 + k / 2 ) );

			p = minIndex( fabs( output( Range::all(), 0, 67 + k / 2 ) ) );
			start = TinyVector< int, 2 >(     p[firstDim], input.lbound(secondDim) );
			end	  = TinyVector< int, 2 >( start[firstDim], input.ubound(secondDim) );

			this->computeGeodesic( slice, start, end, ipath );
		//writer.write(slice);
		//writer.write(ipath);
			slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), upBound );

			// from \phi = 3 \phi / 4 backward
			slice.reference( input(  Range::all(), Range::all(), 67 - k / 2 ) );
			ipath.reference( output( Range::all(), Range::all(), 67 - k / 2 ) );

			p = minIndex( fabs( output( Range::all(), 0, 67 - k / 2 ) ) );
			start = TinyVector< int, 2 >(     p[firstDim], input.lbound(secondDim) );
			end	  = TinyVector< int, 2 >( start[firstDim], input.ubound(secondDim) );

			this->computeGeodesic( slice, start, end, ipath );
		//writer.write(slice);
		//writer.write(ipath);
			slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), upBound );
		}
	}

	writer.write(output);

	return;

	input( Range::all(), 0, Range::all() ) = output( Range::all(), 0, Range::all() );
	input( Range::all(), input.ubound(secondDim), Range::all() ) = output( Range::all(), input.ubound(secondDim), Range::all() );

	for ( int k = 1; k < input.extent(secondDim) - 1; k++ ){
		input( Range::all(), k, Range::all() ) = ( output( Range::all(), k - 1, Range::all() ) + output( Range::all(), k, Range::all() ) + output( Range::all(), k + 1, Range::all() ) ) / 3.0;
	}

	output( Range::all(), Range::all(), 0 ) = input( Range::all(), Range::all(), 0 );
	output( Range::all(), Range::all(), input.ubound(thirdDim) ) = input( Range::all(), Range::all(), input.ubound(thirdDim) );
	for ( int k = 1; k < input.extent(thirdDim) - 1; k++ ){
		output( Range::all(), Range::all(), k ) = ( input( Range::all(), Range::all(), k - 1 ) + input( Range::all(), Range::all(), k ) + input( Range::all(), Range::all(), k + 1 ) ) / 3.0;
	}

	//for ( int k = 0; k < input.extent(secondDim); k++ ){
	//	cout << k << endl;
	//	Array< Pixel, 2 > slice, ipath;
	//	slice.reference( input(  Range::all(), k, Range::all() ) );
	//	ipath.reference( output( Range::all(), k, Range::all() ) );
	//	this->computeGeodesic( slice, startP, endP, ipath );
	//}

	writer.write( input );
	writer.write( output );

	return;

#elif 1

	TinyVector< int, 1 > p;

	start = startP;
	end	  = endP;

	// compute boundary conditions

	Array< Pixel, 2 > ipath;

	slice.reference( input(  Range::all(), input.lbound(secondDim), Range::all() ) );
	ipath.reference( output( Range::all(), input.lbound(secondDim), Range::all() ) );
	this->computeGeodesic( slice, startP, endP, ipath, 1.0 );
	//slice = where( fabs(ipath) < 1, slice, numeric_limits<Pixel>::max() );

	Array< Pixel, 2 > sum( slice.shape() );
	sum = ( input(  Range::all(), input.ubound(secondDim) / 2, Range::all() ) + 
		    input(  Range::all(), input.ubound(secondDim) / 2 + 1, Range::all() ) ) / 2.0;
	ipath.reference( output( Range::all(), input.ubound(secondDim) / 2, Range::all() ) );
	this->computeGeodesic( sum, startP, endP, ipath, 1.0 );
	//slice.reference( input(  Range::all(), input.ubound(secondDim) / 2, Range::all() ) );
	//slice = where( fabs(ipath) < 1, slice, numeric_limits<Pixel>::max() );
	//slice.reference( input(  Range::all(), input.ubound(secondDim) / 2 + 1, Range::all() ) );
	//slice = where( fabs(ipath) < 1, slice, numeric_limits<Pixel>::max() );
	output(  Range::all(), input.ubound(secondDim) / 2 + 1, Range::all() ) = ipath;

	slice.reference( input(  Range::all(), input.ubound(secondDim), Range::all() ) );
	ipath.reference( output( Range::all(), input.ubound(secondDim), Range::all() ) );
	this->computeGeodesic( slice, startP, endP, ipath, 1.0 );
	//slice = where( fabs(ipath) < 1, slice, numeric_limits<Pixel>::max() );

	//writer.write(input);

	for ( int portion = 0; portion < 2; portion++ ){

		Array< Pixel, 3 > halfInput;
		Array< Pixel, 3 > halfOutput;

		Range Jl(0,input.ubound(secondDim)/2);
		Range Jh(input.ubound(secondDim)/2+1,input.ubound(secondDim));

		if ( portion == 0 ){
			halfInput.reference( input( Range::all(), Jl, Range::all() ) );
			halfOutput.reference( output( Range::all(), Jl, Range::all() ) );
		}
		else{
			halfInput.reference( input( Range::all(), Jh, Range::all() ) );
			halfOutput.reference( output( Range::all(), Jh, Range::all() ) );
		}

		//cout << halfInput.shape() << endl;

		vtkPriorityQueue * pqueueY = vtkPriorityQueue::New();
		pqueueY->Allocate( halfInput.cols() );

		for ( int k = 1; k < halfInput.extent(secondDim) - 1; k++ ){
			Array< Pixel, 2 > slice, distance;
			//slice.reference( halfInput(  Range::all(), k, Range::all() ) );
			slice.resize( halfInput.rows(), halfInput.depth() );
			slice = halfInput(  Range::all(), k, Range::all() );
			distance.resize( slice.shape() );
			//nbfFastMarching2D8< Pixel > fm( slice );
			nbfFastFastMarching2D< Pixel > fm( slice );
			vector< TinyVector< int, 2 > > aliveP;
			vector< Pixel > aliveD;
			aliveP.push_back( startP );
			aliveD.push_back(0);
			fm.setAliveSet( aliveP, aliveD );
			fm.setStopPoint( endP );
			fm.execute(distance);
			nbfGeodesicPath< Pixel > geodesic(distance);
			vector< TinyVector< Pixel, 2 > > path;
			geodesic.getSimplePath( distance, endP, path );
			//cout << distance(endP) / path.size() << endl;
			pqueueY->Insert( distance(endP) / path.size(), k );
		}

		vtkPriorityQueue * pqueueZ = vtkPriorityQueue::New();
		pqueueZ->Allocate( halfInput.depth() );

		int boundary = floor( 200.0 / ( startP[0] + endP[0] ) );

		//cout << "b = " << boundary << endl;

		for ( int k = 0; k < halfInput.extent(thirdDim); k++ ){
			if ( ( k < boundary ) | ( k > halfInput.extent(thirdDim) - boundary ) ){
				Array< Pixel, 2 > slice, ipath;
				slice.reference( halfInput(  Range::all(), Range::all(), k ) );
				ipath.reference( halfOutput( Range::all(), Range::all(), k ) );
				p = minIndex( fabs( halfOutput( Range::all(), 0, k ) ) );
				start = TinyVector< int, 2 >( p[firstDim], halfInput.lbound(secondDim) );
				p = minIndex( fabs( halfOutput( Range::all(), halfOutput.ubound(secondDim), k ) ) );
				end	  = TinyVector< int, 2 >( p[firstDim], halfInput.ubound(secondDim) );
				slice = where( slice < numeric_limits< Pixel > :: max(), 1.0, slice );
				this->computeGeodesic( slice, start, end, ipath, 1.0 );
				//slice = lowBound + fabs(ipath);
				//slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), upBound );
				slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), 2 * lowBound );
				//slice = where( fabs(ipath) < 1, slice, meanV );
			}
			else{
				Array< Pixel, 2 > slice, distance;
				//slice.reference( halfInput(  Range::all(), Range::all(), k ) );
				slice.resize( halfInput.rows(), halfInput.cols() );
				slice = halfInput( Range::all(), Range::all(), k );
				distance.resize( slice.shape() );
				//nbfFastMarching2D8< Pixel > fm( slice );
				nbfFastFastMarching2D< Pixel > fm( slice );

				p = minIndex( fabs( halfOutput( Range::all(), 0, k ) ) );
				start = TinyVector< int, 2 >( p[firstDim], halfInput.lbound(secondDim) );

				p = minIndex( fabs( halfOutput( Range::all(), halfOutput.ubound(secondDim), k ) ) );
				end	  = TinyVector< int, 2 >( p[firstDim], halfInput.ubound(secondDim) );

				vector< TinyVector< int, 2 > > aliveP;
				vector< Pixel > aliveD;
				aliveP.push_back( start );
				aliveD.push_back(0);
				fm.setAliveSet( aliveP, aliveD );
				fm.setStopPoint( end );
				fm.execute(distance);
				nbfGeodesicPath< Pixel > geodesic(distance);
				vector< TinyVector< Pixel, 2 > > path;
				geodesic.getSimplePath( distance, end, path );
				//cout << distance(end) / path.size() << endl;
				pqueueZ->Insert( distance(end) / path.size(), k );
			}
		}

		//writer.write(halfInput);
		//writer.write(halfOutput);

		//while ( pqueueZ->GetNumberOfItems() > 0 ){
		//	float vdistance;
		//	int id = pqueueZ->Pop(0,vdistance);
		//	cout << id << " - " << vdistance << endl;
		//}

		// traverse all slices
		double distanceY = numeric_limits<Pixel>::max();
		double distanceZ = numeric_limits<Pixel>::max();
		while( ( pqueueY->GetNumberOfItems() > 0 ) | ( pqueueZ->GetNumberOfItems() > 0 ) ){

			int kY = 0;
			int kZ = 0;

			if ( ( pqueueY->GetNumberOfItems() > 0 ) & ( distanceY == numeric_limits<Pixel>::max() ) ){
				kY = pqueueY->Pop(0,distanceY);
			}

			if ( ( pqueueZ->GetNumberOfItems() > 0 ) & ( distanceZ == numeric_limits<Pixel>::max() ) ){
				kZ = pqueueZ->Pop(0,distanceZ);
			}

			if ( distanceY < distanceZ ){
				Array< Pixel, 2 > slice, ipath;
				slice.reference( halfInput(  Range::all(), kY, Range::all() ) );
				ipath.reference( halfOutput( Range::all(), kY, Range::all() ) );
				this->computeGeodesic( slice, startP, endP, ipath, 1.0 );
				//slice = lowBound + fabs(ipath);
				//slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), upBound );
				//slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), slice );
				slice = where( fabs(ipath) < 1, fabs(ipath), 1 ) + 1 * lowBound;
				//slice = where( fabs(ipath) < 1, slice, meanV );
				distanceY = numeric_limits<Pixel>::max();
				//cout << "y = " << kY << endl;
			}
			else{

				Array< Pixel, 2 > slice, ipath;
				slice.reference( halfInput(  Range::all(), Range::all(), kZ ) );
				ipath.reference( halfOutput( Range::all(), Range::all(), kZ ) );

				p = minIndex( fabs( halfOutput( Range::all(), 0, kZ ) ) );
				start = TinyVector< int, 2 >( p[firstDim], halfInput.lbound(secondDim) );

				p = minIndex( fabs( halfOutput( Range::all(), halfOutput.ubound(secondDim), kZ ) ) );
				end	  = TinyVector< int, 2 >( p[firstDim], halfOutput.ubound(secondDim) );

				this->computeGeodesic( slice, start, end, ipath, 1.0 );
				//slice = lowBound + fabs(ipath);
				//slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), upBound );
				//slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), slice );
				slice = where( fabs(ipath) < 1, fabs(ipath), 1 ) + 1 * lowBound;
				//slice = where( fabs(ipath) < 1, slice, meanV );
				distanceZ = numeric_limits<Pixel>::max();
				//cout << "z = " << kZ << endl;
			}
			//Array< Pixel, 3 > vtk( halfOutput.shape() );
			//vtk = where( halfInput < numeric_limits<Pixel>::max(), halfInput, min(halfInput) );
			//writer.write(vtk);
			//writer.write(halfOutput);
		}
	}

	//writer.write(input);
	//writer.write(output);

	//return;

	for ( int k = 1; k < input.extent(secondDim) - 1; k++ ){
		Array< Pixel, 2 > slice, ipath;
		slice.reference( input(  Range::all(), k, Range::all() ) );
		ipath.reference( output( Range::all(), k, Range::all() ) );
		this->computeGeodesic( slice, startP, endP, ipath );
		// slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), 2 * lowBound );
		slice = lowBound + fabs(ipath);
	}

	for ( int k = 0; k < input.extent(thirdDim); k++ ){
		Array< Pixel, 2 > slice, ipath;
		slice.reference( input(  Range::all(), Range::all(), k ) );
		ipath.reference( output( Range::all(), Range::all(), k ) );
		p = minIndex( fabs( output( Range::all(), 0, k ) ) );
		start = TinyVector< int, 2 >( p[firstDim], input.lbound(secondDim) );
		end	  = TinyVector< int, 2 >( p[firstDim], output.ubound(secondDim) );
		this->computeGeodesic( slice, start, end, ipath );
	}

	//writer.write(output);
	return;

	//input( Range::all(), 0, Range::all() ) = where( fabs( output( Range::all(), 0, Range::all() ) ) < 1, input( Range::all(), 0, Range::all() ), numeric_limits<Pixel>::max() );

	//input( Range::all(), input.ubound(secondDim), Range::all() ) = where( fabs( output( Range::all(), input.ubound(secondDim), Range::all() ) ) < 1, input( Range::all(), input.ubound(secondDim), Range::all() ), numeric_limits<Pixel>::max() );

	//input( Range::all(), input.ubound(secondDim) / 2, Range::all() ) = where( fabs( output( Range::all(), input.ubound(secondDim) / 2, Range::all() ) ) < 1, input( Range::all(), input.ubound(secondDim) / 2, Range::all() ), numeric_limits<Pixel>::max() );
	//	
#else

	Array< Pixel, 2 > ipath;

	slice.reference( input(  Range::all(), input.lbound(secondDim), Range::all() ) );
	ipath.reference( output( Range::all(), input.lbound(secondDim), Range::all() ) );
	this->computeGeodesic( slice, startP, endP, ipath );
	slice = where( fabs(ipath) < 1, slice, numeric_limits<Pixel>::max() );

	slice.reference( input(  Range::all(), input.ubound(secondDim) / 2, Range::all() ) );
	ipath.reference( output( Range::all(), input.ubound(secondDim) / 2, Range::all() ) );
	this->computeGeodesic( slice, startP, endP, ipath );
	slice = where( fabs(ipath) < 1, slice, numeric_limits<Pixel>::max() );

	slice.reference( input(  Range::all(), input.ubound(secondDim), Range::all() ) );
	ipath.reference( output( Range::all(), input.ubound(secondDim), Range::all() ) );
	this->computeGeodesic( slice, startP, endP, ipath );
	slice = where( fabs(ipath) < 1, slice, numeric_limits<Pixel>::max() );

	Array< Pixel, 3 > implicit( input.shape() );
	Array< bool, 3 > confidence( input.shape() );
	//this->computeConfidence( output, thirdDim, 1.0, confidence );
	//writer.write(confidence);

	//input = where( fabs(output) < 1.0, input / 1.1, input ) + 1;
	//writer.write(input);

	for ( int k = 0; k < input.extent(secondDim); k++ ){
		//cout << k << endl;
		Array< Pixel, 2 > slice, ipath;
		slice.reference( input(  Range::all(), k, Range::all() ) );
		ipath.reference( implicit( Range::all(), k, Range::all() ) );
		this->computeGeodesic( slice, startP, endP, ipath );
		cout << k << endl;
	}

	writer.write(implicit);

	for ( int k = 0; k < input.extent(thirdDim); k++ ){
		//cout << k << endl;
		Array< Pixel, 2 > slice, ipath;
		slice.reference( input(  Range::all(), Range::all(), k ) );
		ipath.reference( implicit( Range::all(), Range::all(), k ) );

		//Array< Pixel, 2 > averaged( slice.shape() );
		//averaged = slice + fabs( k - 44.0 ) / 10.0;
		//if ( k > 0 ){
		//	averaged += input(Range::all(), Range::all(), k - 1);
		//	// averaged += input(Range::all(), Range::all(), k - 2);
		//}
		//if ( k < input.ubound(thirdDim) ){
		//	averaged += input(Range::all(), Range::all(), k + 1);
		//	//averaged += input(Range::all(), Range::all(), k + 2);
		//}
		//averaged = averaged / 3;
		//slice.reference( averaged );

		TinyVector< int, 1 > p = minIndex( fabs( output( Range::all(), 0, k ) ) );
		start = TinyVector< int, 2 >(     p[firstDim], input.lbound(secondDim) );
		end	  = TinyVector< int, 2 >( start[firstDim], input.ubound(secondDim) );
		this->computeGeodesic( slice, start, end, ipath );
		cout << k << endl;

		//nbfGeodesicPath< Pixel > geodesic( slice );
		//geodesic.getCircularPath( slice, path );
		//geodesic.getImplicitPath( path, ipath, 2.0 );
	}

	writer.write(implicit);
	Array< bool, 3 > confidenceII( confidence.shape() );
	this->computeConfidence( implicit, thirdDim, 1.0, confidenceII );
	confidence = where( ( confidence == true ) & ( confidenceII == true ), true, false );
	writer.write(confidence);
	
#endif

	//
	//implicit3d.resize( this->weights3d.shape() );

	//int band = 2;

	//// compute phi geodesics in slices
	//// store temporal paths
	//vector< TinyVector< Pixel, 2 > > path;

	//vector< Pixel > startPoints, endPoints;
	//TinyVector< Pixel, 2 > startPoint, endPoint;

	//for ( int k = 0; k < input.extent(secondDim); k++ ){
	//	cout << k << endl;
	//	Array< Pixel, 2 > slice, ipath;
	//	slice.reference( input(  Range::all(), k, Range::all() ) );
	//	ipath.reference( output( Range::all(), k, Range::all() ) );
	//	nbfGeodesicPath< Pixel > geodesic( slice );
	//	geodesic.getFirstPointOnRightSide( slice, startPoint, endPoint );
	//	startPoints.push_back( floor( startPoint[firstDim] ) );
	//	endPoints.push_back( floor( endPoint[firstDim] ) );
	//}

	//vector< Pixel > :: iterator first = startPoints.begin();       
	//vector< Pixel > :: iterator last  = startPoints.end();       
	//sort( first, last );
	//startPoint = TinyVector< Pixel, 2 >( startPoints[ floor( startPoints.size() / 2.0 ) ], input.lbound(thirdDim) );
	//first = endPoints.begin();
	//last  = endPoints.end();
	//sort( first, last );
	//endPoint = TinyVector< Pixel, 2 >( endPoints[ floor( endPoints.size() / 2.0 ) ], input.ubound(thirdDim) );
	//
	//cout << startPoint << endl;
	//cout << endPoint << endl;

	//// restrict domain
	//for ( int i = 0; i < 10; i++ ){
	//	input( Range(fromStart,startPoint[firstDim]-i-1), Range::all(), i ) = numeric_limits<Pixel>::max();
	//	input( Range(startPoint[firstDim]+i+1,toEnd), Range::all(), i ) = numeric_limits<Pixel>::max();
	//	input( Range(fromStart,endPoint[firstDim]-i-1), Range::all(), input.ubound(thirdDim)-i ) = numeric_limits<Pixel>::max();
	//	input( Range(endPoint[firstDim]+i+1,toEnd), Range::all(), input.ubound(thirdDim)-i ) = numeric_limits<Pixel>::max();
	//}

	//nbfMatlabWriter writer;
	//writer.setFileName("ipath");
	//writer.write(input);

	//vector< TinyVector< int, 2 > > aliveP;
	//vector< Pixel > aliveD;
	//aliveP.push_back(startPoint);
	//aliveD.push_back(0);

	//for ( int k = 0; k < input.extent(secondDim); k++ ){
	//	cout << k << endl;
	//	Array< Pixel, 2 > slice, ipath;
	//	slice.reference( input(  Range::all(), k, Range::all() ) );
	//	ipath.reference( output( Range::all(), k, Range::all() ) );
	//	Array< Pixel, 2 > distances( slice.shape() );
	//	nbfFastMarching2D8< Pixel > fm( slice );
	//	fm.setAliveSet(aliveP,aliveD);
	//	TinyVector< int, 2 > stopFM = floor( endPoint );
	//	fm.setStopPoint( stopFM );
	//	fm.execute(distances);
	//	//writer.write(distances);
	//	nbfGeodesicPath< Pixel > geodesic( slice );
	//	geodesic.getSimplePath( distances, stopFM, path);
	//	geodesic.getImplicitPath( path, ipath, band + 1 );
	//}

	//// writer.write(input);
	//writer.write(output);

	//// compute confidence from Z slices
 //	this->computeConfidence( output, secondDim, band, confidence );

	//input = where( ( confidence == true ) & ( input < numeric_limits<Pixel>::max() ), input / 10.0, input );
	//writer.write(input);

	//int narrowBand = 3;

	//for ( int k = 0; k < input.extent(thirdDim); k++ ){
	//	cout << k << endl;
	//	Array< Pixel, 2 > slice, ipath;
	//	slice.reference( input(  Range::all(), Range::all(), k ) );
	//	ipath.reference( output( Range::all(), Range::all(), k ) );

	//	//if ( k > 0 ){
	//	//	slice = where( fabs( output( Range::all(), Range::all(), k - 1 ) ) <= narrowBand, slice, numeric_limits<Pixel>::max() );
	//	//}
	//

	//	//writer.write(slice);
	//	nbfGeodesicPath< Pixel > geodesic( slice );
	//	geodesic.getCircularPath( slice, path );
	//	geodesic.getImplicitPath( path, ipath );
	//	//writer.write(ipath);

	//}

	//writer.write(output);

	//// output = where( confidence == true, 1, 0 );

}


template< class Pixel >
void nbfMinimalSurface< Pixel > :: searchGUI( 
	Array< Pixel, 3 > & input,
	Array< Pixel, 3 > & output,
	Pixel omega )
{
	nbfMatlabWriter writer;
	writer.setFileName("joder.array");

	vector< Pixel > positions;
	vector< Pixel > midPositions;

	vector< TinyVector< Pixel, 2 > > path;

	Array< Pixel, 2 > slice( input.rows(), input.depth() * 2 );
	Array< Pixel, 2 > part( input.rows(), input.depth() );
	
	// First, robustly compute boundary condition
	// make image periodic in the \phi direction to compute closed geodesics in \theta_i
	part = input( Range::all(), 0, Range::all() );
	part.reverseSelf(secondDim);
	slice( Range::all(), Range(fromStart,input.depth()-1) ) = part;
	slice( Range::all(), Range(input.depth(),toEnd) )		= input( Range::all(), floor( input.ubound(secondDim) / 2.0 ), Range::all() );
	nbfGeodesicPath< Pixel > geodesic( slice );
	geodesic.getCircularPath( slice, path );

	vector< TinyVector< Pixel, 2 > > :: iterator iter = path.begin();
	while ( (*iter)[secondDim] > input.depth() ){
		++iter;
	}

	TinyVector< int, 2 > endP( floor(path[0][firstDim]), input.ubound(thirdDim) );
	TinyVector< int, 2 > startP( floor((*iter)[firstDim]), input.lbound(thirdDim) );

	output.resize( input.shape() );

	// AVOID!
	Pixel lowBound = omega;		// min(input);
	Pixel  upBound = omega + 1; // max(input);

	TinyVector< int, 2 > start, end;

	TinyVector< int, 1 > p;

	start = startP;
	end	  = endP;

	// compute boundary conditions

	Array< Pixel, 2 > ipath;

	slice.reference( input(  Range::all(), input.lbound(secondDim), Range::all() ) );
	ipath.reference( output( Range::all(), input.lbound(secondDim), Range::all() ) );
	this->computeGeodesic( slice, startP, endP, ipath, 1.0 );

	//writer.write(slice);
	//writer.write(ipath);

	slice.reference( input(  Range::all(), floor( input.ubound(secondDim) / 2.0 ), Range::all() ) );
	//Array< Pixel, 2 > sum( slice.shape() );
	//sum = ( input(  Range::all(), floor( input.ubound(secondDim) / 2.0 ), Range::all() ) + 
	//	    input(  Range::all(), floor( input.ubound(secondDim) / 2.0 ) + 1, Range::all() ) ) / 2.0;
	ipath.reference( output( Range::all(), floor( input.ubound(secondDim) / 2.0 ), Range::all() ) );
	this->computeGeodesic( slice, startP, endP, ipath, 1.0 );
	output(  Range::all(), floor( input.ubound(secondDim) / 2.0 ) + 1, Range::all() ) = ipath;
	
	slice.reference( input(  Range::all(), input.ubound(secondDim), Range::all() ) );
	ipath.reference( output( Range::all(), input.ubound(secondDim), Range::all() ) );
	this->computeGeodesic( slice, startP, endP, ipath, 1.0 );


	//writer.write(output);
	//AfxMessageBox("Inner file written");

	vector< int > starts;
	vector< int > ends;

	//input = input + 1;
	//lowBound += 1;
	//upBound  += 1;

	for ( int portion = 0; portion < 2; portion++ ){

		Array< Pixel, 3 > halfInput;
		Array< Pixel, 3 > halfOutput;

		Range Jl( 0, floor( input.ubound(secondDim) / 2.0 ) );
		Range Jh( floor( input.ubound(secondDim) / 2.0 ) + 1, input.ubound(secondDim));

		if ( portion == 0 ){
			halfInput.reference( input( Range::all(), Jl, Range::all() ) );
			halfOutput.reference( output( Range::all(), Jl, Range::all() ) );
		}
		else{
			halfInput.reference( input( Range::all(), Jh, Range::all() ) );
			halfOutput.reference( output( Range::all(), Jh, Range::all() ) );
		}

		vtkPriorityQueue * pqueueY = vtkPriorityQueue::New();
		pqueueY->Allocate( halfInput.cols() );

		Array< Pixel, 2 > slice, distance;
		slice.resize( halfInput.rows(), halfInput.depth() );
		distance.resize( slice.shape() );

		//nbfFastMarching2D8< Pixel > fm( slice );
		nbfFastFastMarching2D< Pixel > fm( slice );
		
		vector< TinyVector< int, 2 > > aliveP;
		vector< Pixel > aliveD;
		aliveP.push_back( startP );
		aliveD.push_back(0);
		
		fm.setAliveSet( aliveP, aliveD );
		fm.setStopPoint( endP );
		
		vector< TinyVector< Pixel, 2 > > path;
		nbfGeodesicPath< Pixel > geodesic(distance);
		
		for ( int k = 1; k < halfInput.ubound(secondDim); k++ ){
			slice = halfInput(  Range::all(), k, Range::all() );
			fm.execute(distance);
			geodesic.getSimplePath( distance, endP, path );
			pqueueY->Insert( distance(endP) / path.size(), k );
		}

		vtkPriorityQueue * pqueueZ = vtkPriorityQueue::New();
		pqueueZ->Allocate( halfInput.depth() );

		int boundary = floor( 200.0 / ( startP[0] + endP[0] ) );
		//boundary = halfInput.ubound(secondDim) / 5;
		boundary = 1;

		slice.resize( halfInput.rows(), halfInput.cols() );
		distance.resize( slice.shape() );
		for ( int k = 0; k <= halfInput.ubound(thirdDim); k++ ){
			if ( ( k < boundary ) || ( k > halfInput.ubound(thirdDim) - boundary ) ){
				slice.reference( halfInput(  Range::all(), Range::all(), k ) );
				ipath.reference( halfOutput( Range::all(), Range::all(), k ) );

				// if first pass, compute and store positions
				if ( portion == 0 ){
					p = minIndex( fabs( halfOutput( Range::all(), 0, k ) ) );
					start = TinyVector< int, 2 >( p[firstDim], halfInput.lbound(secondDim) );
					starts.push_back( p[firstDim] );
					p = minIndex( fabs( halfOutput( Range::all(), halfOutput.ubound(secondDim), k ) ) );
					end	  = TinyVector< int, 2 >( p[firstDim], halfInput.ubound(secondDim) );
					ends.push_back( p[firstDim] );
				}
				else{
					start = TinyVector< int, 2 >( ends[k], halfInput.lbound(secondDim) );
					end = TinyVector< int, 2 >( starts[k], halfInput.ubound(secondDim) );
				}

				slice = where( slice < numeric_limits< Pixel > :: max(), 1.0, slice );
				this->computeGeodesic( slice, start, end, ipath, 1.0 );
				slice = where( fabs(ipath) < 1, lowBound + fabs(ipath), 2 * lowBound );
			}
			else{
				slice = halfInput( Range::all(), Range::all(), k );

				// store positions if first pass
				if ( portion == 0 ){
					p = minIndex( fabs( halfOutput( Range::all(), 0, k ) ) );
					start = TinyVector< int, 2 >( p[firstDim], halfInput.lbound(secondDim) );
					starts.push_back( p[firstDim] );

					p = minIndex( fabs( halfOutput( Range::all(), halfOutput.ubound(secondDim), k ) ) );
					end	  = TinyVector< int, 2 >( p[firstDim], halfInput.ubound(secondDim) );
					ends.push_back( p[firstDim] );
				}
				else{
					start = TinyVector< int, 2 >( ends[k], halfInput.lbound(secondDim) );
					end = TinyVector< int, 2 >( starts[k], halfInput.ubound(secondDim) );
				}

				aliveP.clear(); aliveD.clear();
				aliveP.push_back( start );
				aliveD.push_back(0);
				fm.setAliveSet( aliveP, aliveD );
				fm.setStopPoint( end );
				fm.execute(distance);
				geodesic.getSimplePath( distance, end, path );
				//cout << distance(end) / path.size() << endl;
				pqueueZ->Insert( distance(end) / path.size(), k );
			}
		}

		// compensate for distortion around poles of spherical transformation
		thirdIndex kk;
		//halfInput = halfInput * pow2( sin( 1.0f * kk / ( halfInput.depth() - 1.0f ) * vtkMath::Pi() ) ) + lowBound;
		halfInput = halfInput + upBound * ( 1 - pow2( sin( 1.0f * kk / ( halfInput.depth() - 1.0f ) * vtkMath::Pi() ) ) );
		//lowBound += 1;
		//upBound  += 1;
		//halfInput = 1.0f * ( 1.0f - ( sin( 1.0f * kk / halfInput.depth() * vtkMath::Pi() ) ) );

		//writer.write(halfInput);
		//return;

		// traverse all slices
		double distanceY = numeric_limits<Pixel>::max();
		double distanceZ = numeric_limits<Pixel>::max();
		while( ( pqueueY->GetNumberOfItems() > 0 ) || ( pqueueZ->GetNumberOfItems() > 0 ) ){

			int kY;
			int kZ;

			if ( ( pqueueY->GetNumberOfItems() > 0 ) && ( distanceY == numeric_limits<Pixel>::max() ) ){
				kY = pqueueY->Pop(0,distanceY);
			}

			if ( ( pqueueZ->GetNumberOfItems() > 0 ) && ( distanceZ == numeric_limits<Pixel>::max() ) ){
				kZ = pqueueZ->Pop(0,distanceZ);
			}

			if ( distanceY < distanceZ ){
				slice.reference( halfInput(  Range::all(), kY, Range::all() ) );
				ipath.reference( halfOutput( Range::all(), kY, Range::all() ) );
				this->computeGeodesic( slice, startP, endP, ipath, 1.0 );
				slice = where( ( fabs(ipath) < 1 ) || ( slice == numeric_limits<Pixel>::max() ), fabs(ipath), 1 ) + 1 * lowBound;
				distanceY = numeric_limits<Pixel>::max();
			}
			else{
				slice.reference( halfInput(  Range::all(), Range::all(), kZ ) );
				ipath.reference( halfOutput( Range::all(), Range::all(), kZ ) );

				if ( portion == 0 ){
					start = TinyVector< int, 2 >( starts[kZ], halfInput.lbound(secondDim) );
					end = TinyVector< int, 2 >( ends[kZ], halfInput.ubound(secondDim) );
				}
				else{
					start = TinyVector< int, 2 >( ends[kZ], halfInput.lbound(secondDim) );
					end = TinyVector< int, 2 >( starts[kZ], halfInput.ubound(secondDim) );
				}

				this->computeGeodesic( slice, start, end, ipath, 1.0 );
				slice = where( ( fabs(ipath) < 1 ) || ( slice == numeric_limits<Pixel>::max() ), fabs(ipath), 1 ) + 1 * lowBound;
				distanceZ = numeric_limits<Pixel>::max();
			}
		}
	}

	//writer.write(input);
	//AfxMessageBox("File written");

	// second sweep

	for ( int k = 1; k < input.ubound(secondDim); k++ ){
		slice.reference( input(  Range::all(), k, Range::all() ) );
		ipath.reference( output( Range::all(), k, Range::all() ) );
		this->computeGeodesic( slice, startP, endP, ipath );
		slice = lowBound + fabs(ipath);
	}

	for ( int k = 0; k < input.extent(thirdDim); k++ ){
		slice.reference( input(  Range::all(), Range::all(), k ) );
		ipath.reference( output( Range::all(), Range::all(), k ) );
		start = TinyVector< int, 2 >( starts[k], input.lbound(secondDim) );
		end	  = TinyVector< int, 2 >( starts[k], output.ubound(secondDim) );
		this->computeGeodesic( slice, start, end, ipath );
	}

	//for ( int j = 1; j < input.cols() - 1; j++ ){
	//	input( Range::all(), j, Range::all() ) = ( output( Range::all(), j - 1, Range::all() ) +
	//		output( Range::all(), j, Range::all() ) + output( Range::all(), j + 1, Range::all() ) ) / 3.0;
	//}
	//input( Range::all(), input.lbound(secondDim), Range::all() ) = output( Range::all(), input.lbound(secondDim), Range::all() );
	//input( Range::all(), input.ubound(secondDim), Range::all() ) = output( Range::all(), input.ubound(secondDim), Range::all() );

	//for ( int z = 1; z < input.depth() - 1; z++ ){
	//	input( Range::all(), Range::all(), z ) = ( output( Range::all(), Range::all(), z - 1 ) +
	//		output( Range::all(), Range::all(), z ) + output( Range::all(), Range::all(), z + 1 ) ) / 3.0;
	//}
	//input( Range::all(), Range::all(), input.lbound(thirdDim) ) = output( Range::all(), Range::all(), input.lbound(thirdDim) );
	//input( Range::all(), Range::all(), input.ubound(thirdDim) ) = output( Range::all(), Range::all(), input.ubound(thirdDim) );

	return;
}

template< class Pixel >
void nbfMinimalSurface< Pixel > :: search( 
	Array< Pixel, 3 > & input,
	TinyVector< int, 3 > & center,
	Array< Pixel, 3 > & implicit,
	Pixel maxDistance )
{
	implicit.resize( input.shape() );

	// compute circular geodesics in Z slices
	TinyVector< int, 2 > center2D( center[firstDim], center[thirdDim] );

	// store temporal paths
	vector< TinyVector< Pixel, 2 > > path;

	Array< Pixel, 2 > slice, ipath;

	for ( int k = 0; k < input.cols(); k++ ){
		cout << k << endl;
		slice.reference( input(  Range::all(), k, Range::all() ) );
		ipath.reference( implicit( Range::all(), k, Range::all() ) );
		nbfCutGeodesics< Pixel > geodesic( slice, center2D );
		geodesic.getFullCircularPath( slice, path, true );
		if ( maxDistance == numeric_limits<Pixel>::max() ){
			geodesic.getImplicitPath( path, ipath );
		}
		else{
			geodesic.getImplicitPath( path, ipath, maxDistance );
		}
	}

	vector< int > lower, upper;
	for ( int k = 0; k < input.cols(); k++ ){
		lower.push_back( minIndex( fabs( implicit( center[firstDim], k, Range( fromStart, center[thirdDim] ) ) ) )[0] );
		upper.push_back( minIndex( fabs( implicit( center[firstDim], k, Range( center[thirdDim], toEnd ) ) ) )[0] );
	}

	vector< int > :: iterator first = lower.begin();       
	vector< int > :: iterator last  = lower.end();       
	sort( first, last );
	int lowerLimit = lower[ floor( lower.size() / 2.0 ) ];
	first = upper.begin();
	last  = upper.end();
	sort( first, last );
	int upperLimit = upper[ floor( upper.size() / 2.0 ) ] + center[thirdDim];
	
	cout << lowerLimit << endl;
	cout << upperLimit << endl;

	input( center[firstDim], Range::all(), Range( fromStart, lowerLimit - 1 ) ) = numeric_limits<Pixel>::max();
	input( center[firstDim], Range::all(), Range( lowerLimit + 1, upperLimit - 1 ) ) = numeric_limits<Pixel>::max();
	input( center[firstDim], Range::all(), Range( upperLimit + 1, toEnd ) ) = numeric_limits<Pixel>::max();

	nbfMatlabWriter writer;
	writer.setFileName("ipath");
	writer.write(implicit);

	// second pass
	for ( int k = 0; k < input.cols(); k++ ){
		slice.reference( input(  Range::all(), k, Range::all() ) );
		ipath.reference( implicit( Range::all(), k, Range::all() ) );
		//writer.write(slice);
		//writer.write(ipath);
		// recompute only if needed
		if ( ( fabs( ipath( center[firstDim], lowerLimit ) ) > 1 ) | ( fabs( ipath( center[firstDim], upperLimit ) ) > 1 ) ){
			cout << "redo - " << k << endl;
			nbfCutGeodesics< Pixel > geodesic( slice, center2D );
			geodesic.getFullCircularPath( slice, path, true );
			if ( maxDistance == numeric_limits<Pixel>::max() ){
				geodesic.getImplicitPath( path, ipath );
			}
			else{
				geodesic.getImplicitPath( path, ipath, maxDistance );
			}
			//writer.write(ipath);
		}
	}

	writer.write(implicit);
}

template< class Pixel >
void nbfMinimalSurface< Pixel > :: computeGeodesics( 
	Array< Pixel, 3 > & input,
	int dim, int band,
	TinyVector< int, 3 > & center,
	Array< Pixel, 3 > & output )
{
	// store baricenter coordinates
	TinyVector< int, 2 > baricenter;

	// store temporal paths
	vector< TinyVector< Pixel, 2 > > path;

	switch (dim){
		case firstDim:
			baricenter = TinyVector< Pixel, 2 >( center[secondDim], center[thirdDim] );
			break;
		case secondDim:
			baricenter = TinyVector< Pixel, 2 >( center[firstDim], center[thirdDim] );
			break;
		case thirdDim:
			baricenter = TinyVector< Pixel, 2 >( center[firstDim], center[secondDim] );
			break;
	}

	Array< Pixel, 2 > slice, ipath;

	for ( int k = center[dim]; k > -1; k-- ){
		//cout << k << endl;
		switch (dim){
			case firstDim:
				slice.reference( input(  k, Range::all(), Range::all() ) );
				ipath.reference( output( k, Range::all(), Range::all() ) );
				break;
			case secondDim:
				slice.reference( input(  Range::all(), k, Range::all() ) );
				ipath.reference( output( Range::all(), k, Range::all() ) );
				break;
			case thirdDim:
				slice.reference( input(  Range::all(), Range::all(), k ) );
				ipath.reference( output( Range::all(), Range::all(), k ) );
				break;
		}
		nbfCutGeodesics< Pixel > geodesic( slice, baricenter );
		geodesic.getFullCircularPath( slice, path );
		// baricenter = updateBaricenter( baricenter, path );
		// cout << "b = " << baricenter << endl;
		geodesic.getImplicitPath( path, ipath, band + 1 );
		//nbfMatlabWriter writer;
		//writer.setFileName("ipath");
		//writer.write(slice);
		//writer.write(ipath);
	}

	// re-initialize baricenter

	switch (dim){
		case firstDim:
			baricenter = TinyVector< Pixel, 2 >( center[secondDim], center[thirdDim] );
			break;
		case secondDim:
			baricenter = TinyVector< Pixel, 2 >( center[firstDim], center[thirdDim] );
			break;
		case thirdDim:
			baricenter = TinyVector< Pixel, 2 >( center[firstDim], center[secondDim] );
			break;
	}

	for ( int k = center[dim] + 1; k < input.extent(dim); k++ ){
		//cout << k << endl;
		switch (dim){
			case firstDim:
				slice.reference( input(  k, Range::all(), Range::all() ) );
				ipath.reference( output( k, Range::all(), Range::all() ) );
				break;
			case secondDim:
				slice.reference( input(  Range::all(), k, Range::all() ) );
				ipath.reference( output( Range::all(), k, Range::all() ) );
				break;
			case thirdDim:
				slice.reference( input(  Range::all(), Range::all(), k ) );
				ipath.reference( output( Range::all(), Range::all(), k ) );
				break;
		}
		nbfCutGeodesics< Pixel > geodesic( slice, baricenter );
		geodesic.getFullCircularPath( slice, path );
		// baricenter = updateBaricenter( baricenter, path );
		geodesic.getImplicitPath( path, ipath, band + 1 );
	}
}

template< class Pixel >
void nbfMinimalSurface< Pixel > :: computeGeodesics( 
	Array< Pixel, 3 > & input, int band,
	TinyVector< int, 2 > & center,
	Array< Pixel, 3 > & output )
{

}

template< class Pixel >
void nbfMinimalSurface< Pixel > :: computeConfidence( 
	Array< Pixel, 3 > & input,
	int dim, int band,
	Array< bool, 3 > & output )
{
	Array< Pixel, 2 > prevSlice;
	switch (dim){
		case firstDim:
			prevSlice.resize( input.cols(), input.depth() );
			break;
		case secondDim:
			prevSlice.resize( input.rows(), input.depth() );
			break;
		case thirdDim:
			prevSlice.resize( input.rows(), input.cols() );
			break;
	}

	Array< Pixel, 2 > currSlice( prevSlice.shape() );
	Array< Pixel, 2 > nextSlice( prevSlice.shape() );
	Array< bool, 2 > outSlice;

	for ( int k = 1; k < input.ubound(dim); k++ ){
		switch (dim){
			case firstDim:
				prevSlice = fabs( input( k - 1, Range::all(), Range::all() ) );
				currSlice = fabs( input( k, Range::all(), Range::all() ) );
				nextSlice = fabs( input( k + 1, Range::all(), Range::all() ) );
				outSlice.reference( output( k, Range::all(), Range::all() ) );
				break;
			case secondDim:
				prevSlice = fabs( input( Range::all(), k - 1, Range::all() ) );
				currSlice = fabs( input( Range::all(), k, Range::all() ) );
				nextSlice = fabs( input( Range::all(), k + 1, Range::all() ) );
				outSlice.reference( output( Range::all(), k, Range::all() ) );
				break;
			case thirdDim:
				prevSlice = fabs( input( Range::all(), Range::all(), k - 1 ) );
				currSlice = fabs( input( Range::all(), Range::all(), k ) );
				nextSlice = fabs( input( Range::all(), Range::all(), k + 1 ) );
				outSlice.reference( output( Range::all(), Range::all(), k ) );
				break;
		}
		outSlice = where( ( prevSlice < band ) & ( currSlice < band ) & ( nextSlice < band ), true, false );
	}

	// do first slice separately
	switch (dim){
			case firstDim:
				prevSlice = fabs( input( input.lbound(dim)    , Range::all(), Range::all() ) );
				currSlice = fabs( input( input.lbound(dim) + 1, Range::all(), Range::all() ) );
				outSlice.reference( output( input.lbound(dim), Range::all(), Range::all() ) );
				break;
			case secondDim:
				prevSlice = fabs( input( Range::all(), input.lbound(dim)    , Range::all() ) );
				currSlice = fabs( input( Range::all(), input.lbound(dim) + 1, Range::all() ) );
				outSlice.reference( output( Range::all(), input.lbound(dim), Range::all() ) );
				break;
			case thirdDim:
				prevSlice = fabs( input( Range::all(), Range::all(), input.lbound(dim) ) );
				currSlice = fabs( input( Range::all(), Range::all(), input.lbound(dim) + 1 ) );
				outSlice.reference( output( Range::all(), Range::all(), input.lbound(dim) ) );
				break;
	}
	outSlice = where( ( currSlice < band ) & ( nextSlice < band ), true, false );

	// do last slice separately
	switch (dim){
			case firstDim:
				prevSlice = fabs( input( input.ubound(dim) - 1, Range::all(), Range::all() ) );
				currSlice = fabs( input( input.ubound(dim)    , Range::all(), Range::all() ) );
				outSlice.reference( output( input.lbound(dim), Range::all(), Range::all() ) );
				break;
			case secondDim:
				prevSlice = fabs( input( Range::all(), input.ubound(dim) - 1, Range::all() ) );
				currSlice = fabs( input( Range::all(), input.ubound(dim)   , Range::all() ) );
				outSlice.reference( output( Range::all(), input.ubound(dim), Range::all() ) );
				break;
			case thirdDim:
				prevSlice = fabs( input( Range::all(), Range::all(), input.ubound(dim) - 1 ) );
				currSlice = fabs( input( Range::all(), Range::all(), input.ubound(dim) ) );
				outSlice.reference( output( Range::all(), Range::all(), input.ubound(dim) ) );
				break;
	}
	outSlice = where( ( currSlice < band ) & ( nextSlice < band ), true, false );
	//cout << sum( where( outSlice == 1, 1, 0 ) ) << endl;
}


template< class Pixel >
TinyVector< int, 2 > nbfMinimalSurface< Pixel > :: updateBaricenter( 
	TinyVector< int, 2 > & current,
	vector< TinyVector< Pixel, 2 > > & path )
{
	Pixel factor = .99;
	TinyVector< Pixel, 2 > baricenter(0,0);
	vector< TinyVector< Pixel, 2 > > :: iterator iter = path.begin();
	while ( iter != path.end() ){
		baricenter += (*iter);
		++iter;
	}
	baricenter = factor * current + ( 1.0 - factor ) * baricenter / ( 0.0 + path.size() ); 
	return TinyVector< int, 2 >( floor(baricenter[firstDim]), floor(baricenter[secondDim]) );
}

template< class Pixel >
void nbfMinimalSurface< Pixel > :: computeGeodesic( 
	Array< Pixel, 2 > & input,
	TinyVector< int, 2 > & start,
	TinyVector< int, 2 > & end,
	Array< Pixel, 2 > & ipath,
	Pixel maxDistance )
{
	// make local copy because we need contiguously stored arrays
	//Array< Pixel, 2 > w( input.shape() );
	//w = input;

	//nbfFastMarching2D8< Pixel > fm( input );
	nbfFastFastMarching2D< Pixel > fm( input );
	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;
	aliveP.push_back( start );
	aliveD.push_back(0);
	fm.setAliveSet( aliveP, aliveD );
	fm.setStopPoint( end );
	Array< Pixel, 2 > distance( input.shape() );
	fm.execute(distance);

	//nbfMatlabWriter writer;
	//writer.setFileName("ipath");
	//writer.write(input);
	//writer.write(distance);

	nbfGeodesicPath< Pixel > geodesic( distance );
	vector< TinyVector< Pixel, 2 > > path;
	geodesic.getSimplePath( distance, end, path );
	geodesic.getImplicitPath( path, ipath, maxDistance );
	return;
}

#endif /* FILE_nbfMinimalSurface */