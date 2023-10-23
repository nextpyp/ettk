#ifndef FILE_nbfCutMinimalSurface
#define FILE_nbfCutMinimalSurface

// Class nbfCutMinimalSurface.
//
// Find minimal surfaces on cartesian domain. 
// The surface is restricted to be around a given point.
// 

#include <nbfArray.h>
#include <nbfGeodesicPath.h>
#include <cut/nbfCutGeodesics.h>
#include <fm/nbfFastMarching3D26.h>
#include <fm/nbfFastMarching3D.h>
#include <fm/nbfFastMarching2D.h>

using namespace blitz;

#define NBF_DEBUG

template< class Pixel >
class nbfCutMinimalSurface : public nbfArray< Pixel, 3 >
{
public:

	// constructor takes weight array as input and center point coordinates
	nbfCutMinimalSurface( Array< Pixel, 3 > &, TinyVector< int, 3 > & );

	~nbfCutMinimalSurface(){};

	// detach distance computation so we don't need to recompute each time
	void getDistance( Array< Pixel, 3 > & );

	// compute minimal surface. Write implicit surface on first argument and
	// distances on the second argument (neccessary to run the pde)
	void execute( Array< Pixel, 3 > &,  Array< Pixel, 3 > & );
	
protected:

	int geodesics( Array< Pixel, 3 > &, Array< Pixel, 3 > &, int &, TinyVector< int, 2 > &, Pixel &, Pixel &, int &, int &, Array< Pixel, 1 > &, vector< TinyVector< int, 2 > > &, bool = false, bool = false  );
												 
	// store polar distances
	Array< Pixel, 3 > distance3d;

	// store input fast marching weights
	Array< Pixel, 3 > weights3d;

	// store implicit surface representation
	Array< Pixel, 3 > implicit3d;
};


template< class Pixel >
nbfCutMinimalSurface< Pixel > :: nbfCutMinimalSurface( Array< Pixel, 3 > & weights,
												       TinyVector< int, 3 > & center )
													   : nbfArray< Pixel, 3 >( weights, center )
{
	// set weights
	this->weights3d.reference( weights );
}

template< class Pixel >
void nbfCutMinimalSurface< Pixel > :: getDistance( Array< Pixel, 3 > & distance3d )
{
	nbfMatlabWriter writer;
	writer.setFileName("ipath");
#if 1
	// slice by slice 2d fast marching
	distance3d.resize( this->weights3d.shape() );
	distance3d = numeric_limits<Pixel>::max();
	int band = 10;
	int ZuBound = min( center[thirdDim] + band, weights3d.ubound(thirdDim) );
	int ZlBound = max( center[thirdDim] - band, weights3d.lbound(thirdDim) );
	TinyVector< int, 2 > center2D( this->center[firstDim], this->center[secondDim] );
	Array< Pixel, 2 > weights2d( weights3d.rows(), weights3d.cols() );
	Array< Pixel, 2 > rho( weights2d.shape() );
	firstIndex indexI;
	secondIndex indexJ;
	rho = sqrt( pow2( center[0] - indexI + 0.0 ) + pow2( center[1] - indexJ + 0.0 ) );
	rho = where( rho < 3, 0, rho );
	rho = rho / max(rho);

	vector< TinyVector< int, 3 > > aliveP3D;
	vector< Pixel > aliveD3D;

	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;
	for ( int i = center2D[firstDim] + 1; i < weights2d.rows(); i++ ){
		aliveP.push_back( TinyVector<int,2>( i, center2D[secondDim] ) );
		aliveD.push_back(0);
	}
	for ( int z = ZlBound; z <= ZuBound; z++ ){
		weights2d = where( rho > 0, this->weights3d( Range::all(), Range::all(), z ) / rho + 1e-2, numeric_limits<Pixel>::max() );

		// set background points to a fixed big value
		//weights2d = where( this->weights3d( Range::all(), Range::all(), z ) > .99, 2, this->weights3d( Range::all(), Range::all(), z ) );
		//weights2d = where( this->weights3d( Range::all(), Range::all(), z ) < numeric_limits<Pixel>::max(), weights2d + 1e-10, this->weights3d( Range::all(), Range::all(), z ) );
		//weights2d = where( rho < 5, numeric_limits<Pixel>::max(), weights2d );

		Array< Pixel, 2 > distance2d( distance3d( Range::all(), Range::all(), z ) );
		nbfFastMarching2D8< Pixel > fm2( weights2d, center2D );
		fm2.setAliveSet( aliveP, aliveD );
		fm2.setStopBorder(firstDim);
		fm2.execute( distance2d );
		//cout << z << endl;
		//writer.write(weights2d);
		//writer.write(distance2d);

		TinyVector< int, 1 > point = minIndex( distance2d( Range(center2D[firstDim] + 1,toEnd), center2D[secondDim] + 1 ) );
		aliveP3D.push_back( TinyVector< int, 3 >( point[firstDim] + center2D[firstDim] + 1, center2D[secondDim], z ) );
		aliveD3D.push_back(0);

	}
//#elif 0
	// compute cut-transformed weights:

	Array< Pixel, 3 > rho3D( this->weights3d.shape() );
	this->getRho( rho3D );
	Pixel maxRho = max(rho3D);
	// overwrite to use as weights to compute distances
	rho3D = where( rho3D > 0, this->weights3d / rho3D * maxRho, numeric_limits<Pixel>::max() );

	// avoid center pixels
	//this->weights3d = where( phi < 2, numeric_limits<Pixel>::max(), this->weights3d );

	//Array< Pixel, 3 > phi( this->weights3d.shape() );
	//this->getPhi( phi );
	//this->weights3d = this->weights3d / pow2(rho) / sin(phi);
	// fix extremes
	//this->weights3d = where( phi < vtkMath::Pi() / 100, numeric_limits<Pixel>::max(), this->weights3d );
	//this->weights3d = where( phi > vtkMath::Pi() * 99 / 100, numeric_limits<Pixel>::max(), this->weights3d );

	// make sure array shapes are compatible
	distance3d.resize( this->weights3d.shape() );

	// 3d fast marching - compute distances from back face to rest of volume
	nbfFastMarching3D26< Pixel > fastMarching3d( rho3D, this->center );
	vector< TinyVector< int, 3 > > positions;
	vector< Pixel > distances;

	// set to Alive all points in the cut with distance < \infty.
	for ( int k = 0; k < this->weights3d.depth(); k++ ){
		for ( int i = center[firstDim] + 1; i < this->weights3d.rows(); i++ ){
			if ( this->weights3d( i, center[secondDim], k ) < numeric_limits<Pixel>::max() ){
				TinyVector< int, 3 > last( i, center[secondDim], k );
				positions.push_back( last );
				distances.push_back( 0 );
			}
		}
	}
//	fastMarching3d.setAliveSet( positions, distances );
	fastMarching3d.setAliveSet( aliveP3D, aliveD3D );

	// stop when final half-plane reached
//	fastMarching3d.setStopBorder(firstDim);

	fastMarching3d.execute( distance3d );
#else

	distance3d.resize( this->weights3d.shape() );

	//TinyVector< int, 2 > center2D( this->center[firstDim], this->center[secondDim] );
	//Array< Pixel, 2 > weights2d( weights3d.rows(), weights3d.cols() );
	//Array< Pixel, 2 > rho( weights2d.shape() );
	//firstIndex indexI;
	//secondIndex indexJ;
	//rho = sqrt( pow2( center[0] - indexI + 0.0 ) + pow2( center[1] - indexJ + 0.0 ) );
	//rho = where( rho < 3, 0, rho );
	//rho = rho / max(rho);

	//vector< TinyVector< int, 2 > > aliveP;
	//vector< Pixel > aliveD;
	//for ( int i = center2D[firstDim] + 1; i < weights2d.rows(); i++ ){
	//	aliveP.push_back( TinyVector<int,2>( i, center2D[secondDim] ) );
	//	aliveD.push_back(0);
	//}
	//weights2d = where( rho > 0, this->weights3d( Range::all(), Range::all(), center[thirdDim] ) / rho + 1e-2, numeric_limits<Pixel>::max() );

	//Array< Pixel, 2 > distance2d( distance3d( Range::all(), Range::all(), center[thirdDim] ) );
	//nbfFastMarching2D8< Pixel > fm2( weights2d, center2D );
	//fm2.setAliveSet( aliveP, aliveD );
	//fm2.setStopBorder(firstDim);
	//fm2.execute( distance2d );

	//Array< Pixel, 1 > lineA = distance2d( Range( center2D[firstDim], toEnd ), center2D[secondDim] + 1 );
	//TinyVector< int, 1 > mA = minIndex( lineA );
	//TinyVector< int, 2 > start2D( center2D[firstDim] + mA[0], center2D[secondDim] );

	//////////////////
	//
	//// compute cut-transformed weights:

	//Array< Pixel, 3 > rho3d( this->weights3d.shape() );
	//this->getRho3D( rho3d );
	//Pixel maxRho = max(rho3d);
	// overwrite to use as weights to compute distances
	//this->weights3d = where( this->weights3d > .99, 10, this->weights3d );
	//rho3d = where( rho3d > 0, this->weights3d / rho3d * maxRho + 1e-2, numeric_limits<Pixel>::max() );
	//rho3d = where( rho3d > 0, this->weights3d / rho3d * maxRho + 1e-10, numeric_limits<Pixel>::max() );
	//writer.write( rho3d );

	//Array< Pixel, 3 > phi( this->weights3d.shape() );
	//this->getPhi( phi );
	//writer.write(phi);
	//rho3d = where( ( phi > 0 ) & ( phi < vtkMath::Pi() ), rho3d / sin(phi) + 1e-10, rho3d );
	//rho3d = where( rho3d > 3, this->weights3d / pow2(rho3d/maxRho) / sin(phi) + 1e-10, numeric_limits<Pixel>::max() );

	// avoid center pixels
	//this->weights3d = where( phi < 2, numeric_limits<Pixel>::max(), this->weights3d );

	//Array< Pixel, 3 > phi( this->weights3d.shape() );
	//this->getPhi( phi );
	//this->weights3d = this->weights3d / pow2(rho) / sin(phi);
	// fix extremes
	//rho3d = where( phi < vtkMath::Pi() / 100, numeric_limits<Pixel>::max(), rho3d );
	//rho3d = where( phi > vtkMath::Pi() * 99 / 100, numeric_limits<Pixel>::max(), rho3d );

	// consider uniform weights
	distance3d.resize( this->weights3d.shape() );
	Array< Pixel, 3 > uniform( this->weights3d.shape() );
	uniform = 1;
	nbfFastMarching3D< Pixel > fastMarching3dWeights( uniform );
	vector< TinyVector< int, 3 > > positionsW;
	vector< Pixel > distancesW;

	// build alive set to all point with w < 1.
	Array< Pixel, 3 > :: iterator iter = this->weights3d.begin();
	while ( iter != this->weights3d.end() ){
		if ( *iter < 1 ){
			positionsW.push_back( iter.position() );
			distancesW.push_back( *iter );
		}
		++iter;
	}
	fastMarching3dWeights.setAliveSet( positionsW, distancesW );
	fastMarching3dWeights.execute( distance3d );

	// change weights !!!
	//this->weights3d = where( rho3d > 3, this->weights3d / rho3d * maxRho + 1e-10, numeric_limits<Pixel>::max() );
	//writer.write(weights3d);

//	// make sure array shapes are compatible
//	distance3d.resize( this->weights3d.shape() );
//
//	// 3d fast marching - compute distances from back face to rest of volume
//	//nbfFastMarching3D< Pixel > fastMarching3d( rho3d, this->center );
//	nbfFastMarching3D26< Pixel > fastMarching3d( rho3d );
//	vector< TinyVector< int, 3 > > positions;
//	vector< Pixel > distances;
//
//	//TinyVector< int, 3 > start3D( start2D(firstDim), start2D(secondDim), center[thirdDim] );
//
//	//positions.push_back( start3D );
//	//distances.push_back( 0 );
//
//	for ( int k = 0; k < center[thirdDim]; k++ ){
//		TinyVector< int, 3 > start3D( center(firstDim), center(secondDim), k );
//		positions.push_back(start3D);
//		distances.push_back(0);
//	}
//	fastMarching3d.setAliveSet( positions, distances );
//
//	// stop when final half-plane reached
////	fastMarching3d.setStopBorder(firstDim);
//
//	//fastMarching3d.execute( distance3d );
//
//	// SKELETON
//
//	Array< Pixel, 2 > w2d( weights3d.cols(), weights3d.depth() );
//	Array< Pixel, 2 > d2d( w2d.shape() );
//	w2d = numeric_limits< Pixel > :: max();
//	w2d( Range(fromStart,center[secondDim]), Range(center[thirdDim],toEnd) ) = weights3d( center[firstDim], Range(fromStart,center[secondDim]), Range(center[thirdDim],toEnd) );
//	writer.write(w2d);
//
//	nbfFastMarching2D8< Pixel > fm( w2d );
//	aliveP.clear();
//	aliveD.clear();
//	for ( int i = 0; i < center[secondDim]; i++ ){
//		aliveP.push_back( TinyVector< int, 2 >( i, center[thirdDim] ) );
//		aliveD.push_back(0);
//	}
//	fm.setAliveSet( aliveP, aliveD );
//	fm.execute( d2d );
//	writer.write(d2d);


#endif
}

template< class Pixel >
int nbfCutMinimalSurface< Pixel > :: geodesics( Array< Pixel, 3 > & implicit3d,
											     Array< Pixel, 3 > & weights3d,
												 int & firstZ,
												 TinyVector< int, 2 > & centerOnFirstSlice,
												 Pixel & narrowBandStrictIn,
												 Pixel & narrowBandStrictOut,
												 int & bottomZ,
												 int & topZ,
												 Array< Pixel, 1 > & position,
												 vector< TinyVector< int, 2 > > & centers,
												 bool circular = false,
												 bool computeBounds = false )
{
	int bestSlice = firstZ;
	int minPathLenght = numeric_limits<int>::max();

	nbfMatlabWriter writer;
	writer.setFileName("ipath");

	// store geodesic lenght for each slice
	Array< Pixel, 1 > energy( implicit3d.depth() );
	energy = 0;

	// stop when geodesic lenght exceeds threshold value
	Pixel thEnergy = .5;

	Array< Pixel, 2 > weights2d( implicit3d.rows(), implicit3d.cols() );
	nbfCutGeodesics< Pixel > geodesic2d( weights2d, centerOnFirstSlice );	

	centers[firstZ] = centerOnFirstSlice;

	for ( int z = firstZ; z <= topZ; z++ ){
		
		// if INTRA slice, skip slice
		if ( position(z) != 2 ){
			continue;
		}

		//cout << "z = " << z << endl;
		
		weights2d = weights3d( Range::all(), Range::all(), z );

		// restrict path to a band around geodesic on previous slice
		if ( ( 1 | circular ) & ( z > firstZ ) ){
			weights2d = where( - implicit3d( Range::all(), Range::all(), z - 1 ) <= narrowBandStrictIn, weights2d, numeric_limits<Pixel>::max() );
			weights2d = where( implicit3d( Range::all(), Range::all(), z - 1 ) <= narrowBandStrictOut, weights2d, numeric_limits<Pixel>::max() );
		}

		//writer.write(weights2d);

		// store implicit geodesic representation in a 3D volume
		Array< Pixel,2 > ipath( implicit3d( Range::all(), Range::all(), z ) );

		//writer.write( weights2d );

		// re-center according to center on previous slice
		Array< Pixel,2 > ipathPrevious( implicit3d( Range::all(), Range::all(), z - 1 ) );
		centers[z] = minIndex( ipathPrevious );

		// to avoid erratic changes in the position of the center, we do a weighted average
		// with the position in the precious slice and the one before 
		TinyVector< int, 2 > center;
		if ( z > firstZ ){
			center = ( .75 * centers[z-1] + .25 * centers[z] );
		}
		else{
			center = centerOnFirstSlice;
		}
		if ( computeBounds == true ){
			geodesic2d.setCenter( center );
		}
		centers[z] = center;

		vector< TinyVector< int, 2 > > tpath;
		geodesic2d.getFullCircularPath( weights2d, tpath, circular );

		// compute lenght of geodesic
		for ( int h = 0; h < tpath.size(); h++ ){
			energy(z) += this->weights3d( tpath[h](firstDim), tpath[h](secondDim), z );
		}
		energy(z) /= tpath.size();

		if ( ( fabs(firstZ-z+0.0) < 10 ) & ( energy(z) < .4 ) & ( tpath.size() < minPathLenght ) ){
			minPathLenght = tpath.size();
			bestSlice = z;
		}

		geodesic2d.getImplicitPath( tpath, ipath );

		// override
		geodesic2d.getCircularPath( weights2d, ipath );

		Array< Pixel, 2 > tmp( ipath.shape() );
		tmp = ipath;
		geodesic2d.curvature(tmp);
		ipath = tmp;

		// if two consecutive slices have lenght above the threshold => assume the surface has ended
		if ( computeBounds & ( energy(z) > thEnergy ) & ( energy(z-1) > thEnergy ) ){
			topZ = z;
			break;
		}
	}

	for ( int z = firstZ - 1; z >= bottomZ; z-- ){
	
		if ( position(z) != 2 ){
			continue;
		}

		//cout << "z = " << z << endl;

		weights2d = weights3d( Range::all(), Range::all(), z );
		
		// restrict path to a band around geodesic on previous slice
		if ( 1 | circular ){
			weights2d = where( - implicit3d( Range::all(), Range::all(), z + 1 ) <= narrowBandStrictIn, weights2d, numeric_limits<Pixel>::max() );
			weights2d = where( implicit3d( Range::all(), Range::all(), z + 1 ) <= narrowBandStrictOut, weights2d, numeric_limits<Pixel>::max() );
		}

		//writer.write(weights2d);

		// store implicit geodesic representation in a 3D volume
		Array< Pixel,2 > ipath( implicit3d( Range::all(), Range::all(), z ) );

		// re-compute center from previous slices
		Array< Pixel,2 > ipathPrevious( implicit3d( Range::all(), Range::all(), z + 1 ) );
		centers[z] = minIndex( ipathPrevious );
		TinyVector< int, 2 > center = ( .25 * centers[z] + .75 * centers[z+1] );
		if ( computeBounds == true ){
			geodesic2d.setCenter( center );
		}
		centers[z] = center;

		vector< TinyVector< int, 2 > > tpath;
		geodesic2d.getFullCircularPath( weights2d, tpath, circular );

		// compute lenght of geodesic		
		for ( int h = 0; h < tpath.size(); h++ ){
			energy(z) += this->weights3d( tpath[h](firstDim), tpath[h](secondDim), z );
		}
		energy(z) /= tpath.size();

		if ( ( fabs(firstZ-z+0.0) < 10 ) & ( energy(z) < .4 ) & ( tpath.size() < minPathLenght ) ){
			minPathLenght = tpath.size();
			bestSlice = z;
		}

		geodesic2d.getImplicitPath( tpath, ipath );

		// override
		geodesic2d.getCircularPath( weights2d, ipath );

		Array< Pixel, 2 > tmp( ipath.shape() );
		tmp = ipath;
		geodesic2d.curvature(tmp);
		ipath = tmp;

		if ( computeBounds & ( energy(z) > thEnergy ) & ( energy(z+1) > thEnergy ) ){
			bottomZ = z;
			//cout << "z = " << z << endl;
			break;
		}
	}
	writer.write(energy);
	return bestSlice;
}

template< class Pixel >
void nbfCutMinimalSurface< Pixel > :: execute( Array< Pixel, 3 > & implicit3d,
											   Array< Pixel, 3 > & distance3d )
{
	nbfMatlabWriter writer;
	writer.setFileName("ipath");

	// initialize array
	implicit3d.resize( this->weights3d.shape() );
	implicit3d = 1;

	// Even with uniform weights (g=1), the final arrival times will not be the same
	// because of fast marching errors.
	// This becomes a problem to detect the first arriving point.
	// We correct the arriving times re-scaling by values obtained for the uniform case.
#if 0
	// 2D uniform $g$
	Array< Pixel,2 > uniformG( distance3d.rows(), distance3d.cols() );
	uniformG = 1;
	
	// g_r = g / r
	firstIndex i;
	secondIndex j;
	Array< Pixel, 2 > rho( uniformG.shape() );
	rho = sqrt( pow2( center[0] - i + 0.0 ) + pow2( center[1] - j + 0.0 ) );
	Pixel maxRho = max(rho);
	uniformG = where( rho > 0, uniformG / rho * maxRho, numeric_limits<Pixel>::max() );

	// build alive set
	vector< TinyVector< int, 2 > > min1;
	vector< Pixel > ini;
	for ( int i = center[firstDim] + 1; i < uniformG.rows(); i++ ){
		min1.push_back( TinyVector< int, 2 >(i,center[secondDim]) );
		ini.push_back(0);
	}			
	
	// compute cut-distances (same center as in 3D)
	nbfFastMarching2D8< Pixel > fm2( uniformG, TinyVector< int, 2 >( center[0], center[1] ) );
	
	fm2.setAliveSet( min1, ini );
	Array< Pixel, 2 > uniformDistance( uniformG.shape() );
	fm2.execute( uniformDistance );

	// look at arriving distances (1D line)
	// (this is the correction term)
	Array< Pixel, 1 > line( uniformDistance( Range( center[firstDim], toEnd ), center[1] + 1 ) );

	for ( int k = 0; k < B.cols(); k++ ){
		B( Range::all(), k ) = B( Range::all(), k ) / line;
	}
#endif


	TinyVector< int, 2 > centerD( this->center[firstDim], this->center[secondDim] );	
	Array< Pixel, 2 > weightsD( distance3d( Range::all(), Range::all(), center[thirdDim] ) );
	nbfCutGeodesics< Pixel > geodesicD( weightsD, centerD );	
	Array< Pixel, 2 > implicitD( weightsD.shape() );
	geodesicD.getCircularPath( weightsD, implicitD );
	geodesicD.curvature(implicitD);
	writer.write(weightsD);
	writer.write(implicitD);

	// Minimal surface computation

	// Start slice-by-slice processing

	// big enough to allow fast surface termination, small enough to avoid zero minima
	Pixel narrowBandStrictIn = 200;
	// small to force regularity, big enough to follow surface
	Pixel narrowBandStrictOut = 200;

	// first point Z coordinate
	int firstZ = center[thirdDim];

	// tell where the surface starts and ends
	int topZ = weights3d.ubound(thirdDim);
	int bottomZ = weights3d.lbound(thirdDim);

	TinyVector< int, 2 > centerOnFirstSlice( this->center[firstDim], this->center[secondDim] );
	
	vector< TinyVector< int, 2 > > centers( weights3d.depth() );
	centers[firstZ] = centerOnFirstSlice;

	// this is to deal with interpolated slices
	Array< Pixel, 1 > axis( implicit3d.depth() );
	axis = 2;

	// find all circular geodesics in the Z dimension
	int bestZ = this->geodesics( implicit3d, distance3d, firstZ, centerOnFirstSlice, narrowBandStrictIn, narrowBandStrictOut, bottomZ, topZ, axis, centers, true, false );
	writer.write(implicit3d);

	// find robust segmented regions

	Array< Pixel, 3 > zone3d( implicit3d.shape() );
	zone3d = numeric_limits<Pixel>::max();
	Array< Pixel, 2 > prevSlice( implicit3d.rows(), implicit3d.cols() );
	Array< Pixel, 2 > currSlice( prevSlice.shape() );
	Array< Pixel, 2 > nextSlice( prevSlice.shape() );
	Pixel band = 2;
	for ( int z = bottomZ + 1; z < topZ; z++ ){
		prevSlice = fabs( implicit3d( Range::all(), Range::all(), z - 1 ) );
		currSlice = fabs( implicit3d( Range::all(), Range::all(), z ) );
		nextSlice = fabs( implicit3d( Range::all(), Range::all(), z + 1 ) );
		zone3d( Range::all(), Range::all(), z ) = where( ( prevSlice < band ) & ( currSlice < band ) & ( nextSlice < band ), 1, numeric_limits<Pixel>::max() );
	}

	currSlice = fabs( implicit3d( Range::all(), Range::all(), topZ ) );
	zone3d( Range::all(), Range::all(), topZ ) = where( currSlice < band, 1, numeric_limits<Pixel>::max() );
	currSlice = fabs( implicit3d( Range::all(), Range::all(), bottomZ ) );
	zone3d( Range::all(), Range::all(), bottomZ ) = where( currSlice < band, 1, numeric_limits<Pixel>::max() );

	// build new weights
	//zone3d = where( zone3d < 10, this->weights3d, 1 );
	zone3d = where( zone3d < 10, this->weights3d, 10 * this->weights3d );

	writer.write(zone3d);

	// XZ slice
	Array< Pixel, 2 > weightsXZ( zone3d( Range::all(), this->center[secondDim], Range::all() ) );

	TinyVector< int, 2 > centerXZ( this->center[firstDim], this->center[thirdDim] );	
	nbfCutGeodesics< Pixel > geodesicXZ( weightsXZ, centerXZ );	
	Array< Pixel, 2 > implicitXZ( weightsXZ.shape() );
	vector< TinyVector< int, 2 > > path;
	geodesicXZ.getFullCircularPath( weightsXZ, path, true );
	geodesicXZ.getImplicitPath( path, implicitXZ );

	geodesicXZ.curvature(implicitXZ);

	writer.write(weightsXZ);
	writer.write(implicitXZ);

	// YZ slice
	Array< Pixel, 2 > weightsYZ( zone3d( this->center[firstDim], Range::all(), Range::all() ) );

	TinyVector< int, 2 > centerYZ( this->center[secondDim], this->center[thirdDim] );	
	nbfCutGeodesics< Pixel > geodesicYZ( weightsYZ, centerYZ );	
	Array< Pixel, 2 > implicitYZ( weightsYZ.shape() );
	geodesicYZ.getFullCircularPath( weightsYZ, path, true, 1 );
	geodesicYZ.getImplicitPath( path, implicitYZ );

	geodesicXZ.curvature(implicitYZ);

	writer.write(weightsYZ);
	writer.write(implicitYZ);

	weights3d( Range::all(), this->center[secondDim], Range::all() ) = where( fabs( implicitXZ ) < 2, weights3d( Range::all(), this->center[secondDim], Range::all() ), numeric_limits<Pixel>::max() );

	weights3d( this->center[firstDim], Range::all(), Range::all() ) = where( fabs(implicitYZ) < 2, weights3d( this->center[firstDim], Range::all(), Range::all() ), numeric_limits<Pixel>::max() );
	
	writer.write(weights3d);

	//narrowBandStrictIn = 3;
	//narrowBandStrictOut = 3;
	implicit3d = 1;
	distance3d = where( weights3d < numeric_limits<Pixel>::max(), distance3d, numeric_limits<Pixel>::max() );
	distance3d = where( weights3d < numeric_limits<Pixel>::max(), 1, numeric_limits<Pixel>::max() );
	this->geodesics( implicit3d, distance3d, firstZ, centerOnFirstSlice, narrowBandStrictIn, narrowBandStrictOut, bottomZ, topZ, axis, centers, true, false );
	writer.write(implicit3d);

#if 0
	// SLICE INTERPOLATION 

	// \infty - slice out of surface
	// -1 - first intra slice
	//  0 - intermediate intra slice
	//  1 - last intra slice
	//  2 - interpolated slice

	// initialize all points to interpolated
	axis = 2;

	axis(bottomZ) = -1;
	axis(topZ) = 1;
	
	for ( int z = bottomZ + 1; z < topZ; z++ ){
		weights2d = zone3d( Range::all(), Range::all(), z );
		geodesic2d.setCenter( centers[z] );

		vector< TinyVector< int, 2 > > tpath;
		geodesic2d.getFullCircularPath( weights2d, tpath, true );

		if ( tpath.size() == 0 ){ // if not closed
			axis(z) = 2; // interpolated
		}
		else{
			axis(z) = 0; // last intra
		}
	}

	// interpolate between INTRA's and find corresponding geodesic

	// starting and ending points
	Pixel start = bottomZ;
	Pixel last = start;

	zone3d( Range::all(), Range::all(), bottomZ ) = implicit3d( Range::all(), Range::all(), bottomZ );
	zone3d( Range::all(), Range::all(), topZ ) = implicit3d( Range::all(), Range::all(), topZ );

	bool doneInterpolating = false;
	while ( !doneInterpolating ){
		start = last;
		last++;
		while ( axis(last) > 1 ){
			if ( last == axis.ubound(firstDim) ){
				doneInterpolating = true;
				break;
			}
			last++;
		}
		if ( axis(last) == 1 ){
			doneInterpolating = true;
		}
		if ( axis(last) < 2 ){
			for ( Pixel k = start + 1; k < last; k++ ){
				zone3d( Range::all(), Range::all(), k ) = implicit3d( Range::all(), Range::all(), last ) * (k-start)/(last-start) + implicit3d( Range::all(), Range::all(), start ) * ( 1.0 - (k-start)/(last-start) );
			}
		}
	}
#endif

	// build new weights
	for ( int k = bottomZ; k <= topZ; k++ ){
		Array< Pixel, 2 > weights2d( this->weights3d( Range::all(), Range::all(), k ) );
		Array< Pixel, 2 > zone2d( zone3d( Range::all(), Range::all(), k ) );
		weights2d = where( abs(zone2d) > 2, numeric_limits<Pixel>::max(), weights2d );
	}

	// re-compute geodesics on interpolated slices
	writer.write(weights3d);
	this->geodesics( implicit3d, distance3d, firstZ, centers[bottomZ], narrowBandStrictIn, narrowBandStrictOut, bottomZ, topZ, axis, centers, false, false );

	// initialize remaining slices
	for ( int i = bottomZ - 1; i >= 0; i-- ){
		implicit3d( Range::all(), Range::all(), i ) = abs( implicit3d( Range::all(), Range::all(), bottomZ ) ) + ( bottomZ - i );
	}
	for ( int i = topZ + 1; i < implicit3d.depth(); i++ ){
		implicit3d( Range::all(), Range::all(), i ) = abs( implicit3d( Range::all(), Range::all(), topZ ) ) + ( i - topZ );
	}

	// force closed surfaces
	implicit3d( Range::all(), Range::all(), implicit3d.ubound(thirdDim) ) = abs( implicit3d( Range::all(), Range::all(), implicit3d.ubound(thirdDim) - 1 ) ) + 1;
	implicit3d( Range::all(), Range::all(), implicit3d.lbound(thirdDim) ) = abs( implicit3d( Range::all(), Range::all(), implicit3d.lbound(thirdDim) + 1 ) ) + 1;

}

#if 0 // VERSION with cross sections
template< class Pixel >
void nbfCutMinimalSurface< Pixel > :: execute( Array< Pixel, 3 > & implicit3d,
											   Array< Pixel, 3 > & distance3d )
{
	nbfMatlabWriter writer;
	writer.setFileName("ipath");

	// initialize array
	implicit3d.resize( this->weights3d.shape() );
	implicit3d = 1;

	// Minimal surface computation

	// get slice with 3D arriving distances
	Array< Pixel, 2 > B( distance3d( Range::all(),
									 center[secondDim] + 1, 
									 Range::all() ) );
	writer.write(B);

	// work relative to the arriving times slice

	// compute geodesic on arriving slice (assume the object is entirely inside 3D volume)
	TinyVector< int, 2 > centerOnB( center[firstDim], center[thirdDim] );
	nbfCutGeodesics< Pixel > cutGeodesicsB( B, centerOnB );
	vector< TinyVector< int, 2 > > path;

	// B is a slice in the X-Z direction (Y constant)
	// path coordinates are (x,z) pairs, y = center[secondDim] + 1

	// get path surrounding center in arriving slice
	cutGeodesicsB.setCutDimension(secondDim);
	cutGeodesicsB.getFullCircularPath( B, path );

	Array< Pixel,2 > ipath( implicit3d( Range::all(), center[secondDim] + 1, Range::all() ) );
	cutGeodesicsB.getImplicitPath( path, ipath );
	writer.write(ipath);

	// look for crossings with z-axis (to restrict path in cross slice)
	int top = 0;
	int bottom = B.cols();
	for ( int i = 0; i < path.size(); i++ )
	{
		if ( path[i](firstDim) == center[firstDim] ){
			if ( path[i](secondDim) > top ){
				top = path[i](secondDim);
			}
			if ( path[i](secondDim) < bottom ){
				bottom = path[i](secondDim);
			}
		}
	}

	TinyVector< int, 2 > topP( center[secondDim], top );
	TinyVector< int, 2 > bottomP( center[secondDim], bottom );

	cout << topP << endl;
	cout << bottomP << endl;

	// C is a slice in the Y-Z direction (X constant)

	// get cross geodesic passing through top & bottom points
	Array< Pixel, 2 > C( distance3d( center[firstDim], Range::all(), Range::all() ) );
	//C = where( C == numeric_limits< Pixel > :: max(), 4, C );
	writer.write(C);

	TinyVector< int, 2 > centerOnC( center[secondDim], center[thirdDim] );
	nbfCutGeodesics< Pixel > cutGeodesicsOnC( C, centerOnC );
	cutGeodesicsOnC.setCutDimension( secondDim );
	vector< TinyVector< int, 2 > > crossPath;
	cutGeodesicsOnC.getConstrainedCircularPath( C, topP, bottomP, crossPath );

	writer.write(C);
	Array< Pixel, 2 > tmp( C.shape() );
	cutGeodesicsOnC.getImplicitPath( crossPath, tmp );
	writer.write(tmp);

	// here we store the top and bottom positions for each slice
	Array< int, 2 > fixedPoints( distance3d.cols(), 2 );
	// tops
	fixedPoints( Range::all(), 0 ) = 0;
	// bottoms
	fixedPoints( Range::all(), 1 ) = distance3d.depth();

	// B is a slice in the X-Z direction (Y constant)
	// path coordinates are (x,z) pairs, y = center[secondDim] + 1
	// C is a slice in the Y-Z direction (X constant)

	// extract top and bottom points on each slice
	for ( int i = 0; i < crossPath.size(); i++ )
	{
		if ( crossPath[i](secondDim) > fixedPoints(crossPath[i](firstDim),0) ){
			fixedPoints(crossPath[i](firstDim),0) = crossPath[i](secondDim);
		}
		if ( crossPath[i](secondDim) < fixedPoints(crossPath[i](firstDim),1) ){
			fixedPoints(crossPath[i](firstDim),1) = crossPath[i](secondDim);
		}
	}

	// restrict current geodesic to a band around geodesic in the previous slice (+/- narrowband)
	int narrowBand = 5;


	for ( int j = center[secondDim] + 2; j < distance3d.cols(); j++ )
	{
		// store implicit geodesic representation in a 3D volume
		Array< Pixel,2 > ipath( implicit3d( Range::all(), j, Range::all() ) );
		Array< Pixel,2 > weights2d( distance3d( Range::all(), j, Range::all() ) );

		// restrict to band around previous geodesic
		weights2d = where( fabs( implicit3d( Range::all(), j - 1, Range::all() ) ) <= narrowBand, weights2d, numeric_limits<Pixel>::max() );

		writer.write(weights2d);
		if ( ( fixedPoints(j,0) != 0 ) & ( fixedPoints(j,1) < distance3d.depth() ) )
		{
			TinyVector< int, 2 > topP( center[firstDim], fixedPoints(j,0) );
			TinyVector< int, 2 > bottomP( center[firstDim], fixedPoints(j,1) );
			cout << "top = " << topP << endl;
			cout << "bot = " << bottomP << endl;
			TinyVector< int, 2 > newCenter = ( topP + bottomP ) / 2;
			cutGeodesicsB.setCenter( newCenter );
			path.clear();
			cutGeodesicsB.getConstrainedCircularPath( weights2d, topP, bottomP, path );
			cout << "size = " << path.size() << endl;
			cutGeodesicsB.getImplicitPath( path, ipath );
			writer.write(ipath);

		}
		else
		{
			break;
		}
	}

	implicit3d( Range::all(), center[secondDim], Range::all() ) =
		implicit3d( Range::all(), center[secondDim] + 1, Range::all() );

	for ( int j = center[secondDim] - 1; j >= 0; j-- )
	{
		// store implicit geodesic representation in a 3D volume
		Array< Pixel,2 > ipath( implicit3d( Range::all(), j, Range::all() ) );
		Array< Pixel,2 > weights2d( distance3d( Range::all(), j, Range::all() ) );

		// restrict to band around previous geodesic
		weights2d = where( fabs( implicit3d( Range::all(), j + 1, Range::all() ) ) <= narrowBand, weights2d, numeric_limits<Pixel>::max() );

		writer.write(weights2d);
		if ( ( fixedPoints(j,0) != 0 ) & ( fixedPoints(j,1) < distance3d.depth() ) )
		{
			TinyVector< int, 2 > topP( center[firstDim], fixedPoints(j,0) );
			TinyVector< int, 2 > bottomP( center[firstDim], fixedPoints(j,1) );
			cout << "top = " << topP << endl;
			cout << "bot = " << bottomP << endl;
			TinyVector< int, 2 > newCenter = ( topP + bottomP ) / 2;
			cutGeodesicsB.setCenter( newCenter );
			path.clear();
			cutGeodesicsB.getConstrainedCircularPath( weights2d, topP, bottomP, path );
			cutGeodesicsB.getImplicitPath( path, ipath );
			writer.write(ipath);
		}
		else
		{
			break;
		}
	}

	implicit3d( implicit3d.lbound(firstDim), Range::all(), Range::all() ) = 1;
	implicit3d( implicit3d.ubound(firstDim), Range::all(), Range::all() ) = 1;
	implicit3d( Range::all(), Range::all(), implicit3d.lbound(thirdDim) ) = 1;
	implicit3d( Range::all(), Range::all(), implicit3d.ubound(thirdDim) ) = 1;
	implicit3d( Range::all(), implicit3d.lbound(secondDim), Range::all() ) = 1;
	implicit3d( Range::all(), implicit3d.ubound(secondDim), Range::all() ) = 1;

	writer.write(implicit3d);
}
#endif // execute version

#endif /* FILE_nbfCutMinimalSurface */