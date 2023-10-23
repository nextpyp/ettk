#ifndef FILE_nbfFastMarchingCurvature3D
#define FILE_nbfFastMarchingCurvature3D

// Class nbfFastMarchingCurvature.
//
// Implements 3D fast marching method with arbitrary weights.
// 

#include <fm/nbfFastMarching3D.h>
#include <bs/nbfBordStrategyMirror.h>
#include <nbf-stencil-et.h>

#include "io/nbfImageWriter.h"

template< class Pixel, class Scalar = Pixel >
class nbfFastMarchingCurvature3D : public nbfFastMarching3D< Pixel, Scalar >
{
public:

	// constructor takes weight array as input
	nbfFastMarchingCurvature3D( Array< Pixel, 3 > & );
	~nbfFastMarchingCurvature3D(){};

	// override from FastMarching 
	// (this is done only to return the 'state' as the output and not the 'distance' values)
	virtual TinyVector< int, 3 > execute( Array< Scalar, 3 > & );

	// set the smoothness term weight (defaults to .25)
	void setAlpha( Scalar a ){ this->alpha = a; }

	// set neighborhood size (defaults to 5)
	void setNeighborhoodSize( int s ){ this->nsize = s; }

	void getNeighbors( TinyVector< int, 3 > &, vector< TinyVector< int, 3 > > & );

protected:

	// override from FastMarching
	void updatePoint( TinyVector< int, 3 > &, TinyVector< int, 3 > & );

	// override from FastMarching
	virtual void initialize();

	// update global statistics
	void updateGlobalStatistics();

	// update mean and variance in the neighborhood around a point
	void updateLocalStatistics( TinyVector< int, 3 > & );

	// override from FastMarching
	void loopUpdate( TinyVector< int, 3 > & );

	// Statistical members

	// neighborhood size
	int nsize;

	// Mean: inside (piecewise smooth), outside (constant)
	Array< Pixel, 3 > mean;

	// Mean outside ( note that mean(outside) = meanOutside )
	Pixel meanOutside;

	// variance inside/outside
	Scalar sigmaInside, sigmaOutside;

	// current inside region size
	float foregroundSize;

	// current valid background (FA != 0)
	float backgroundSize;

	// size of region when last update of statistics was done
	// this is used to update the global statistics when the region has grow enough
	// (see the loopUpdate method)
	int previousSizeUpdate;

	// Local mean and variance estimates (independent of segmentation)
	Array< Pixel, 3 >  localMean;
	Array< Scalar, 3 > localVariance;
	
	// implicit curve representation
	Array< Scalar, 3 > implicitH;

	// mask background points
	Array< Scalar, 3 > mask;

	// global energy
	Array< Scalar, 3 > energy;
	Scalar totalEnergy;
	Array< Scalar, 1 > energyVector;
	int energyCounter;
	Scalar energyGlobalMinima;

	// smoothness term weight
	Scalar alpha;

};

template< class Pixel, class Scalar >
nbfFastMarchingCurvature3D< Pixel, Scalar > :: nbfFastMarchingCurvature3D( Array< Pixel, 3 > & weight )
: nbfFastMarching3D< Pixel, Scalar >( weight )
{
	this->nsize = 5;
	this->mean.resize( weight.shape() );
	this->localMean.resize( weight.shape() );
	this->localVariance.resize( weight.shape() );

	// stopping criteria
	this->energy.resize( weight.shape() );
	this->totalEnergy = 0;
	this->energyVector.resize( weight.numElements() );
	this->energyCounter = 0;
	this->energyGlobalMinima = numeric_limits<Scalar>::max();

	this->implicitH.resize( weight.shape() );
	this->implicitH = where( weight < numeric_limits< Scalar > :: max(), 1, 0 );
	this->alpha = .25;
	this->foregroundSize = 0;
	this->previousSizeUpdate = this->foregroundSize;

	// initialize background mask and size
	this->mask.resize( weight.shape() );
	//this->mask = 1;
	//Array< Scalar, 3 > :: iterator iterM = this->mask.begin();
	//Array< Pixel, 3 > :: iterator iterW = this->weight.begin();
	//while ( iterM != this->mask.end() ){
	//	if ( dot(*iterW,*iterW) > 0 ){
	//		(*iterM) = 1;
	//	}
	//	else{
	//		(*iterM) = 0;
	//	}
	//	++iterM; ++iterW;
	//}
	this->mask = where( weight < numeric_limits< Scalar > :: max(), 1, 0 );
}

template< class Pixel, class Scalar >
void nbfFastMarchingCurvature3D< Pixel, Scalar > :: initialize()
{
	// starting foreground size is number of points in initial alive set
	this->foregroundSize = this->aliveSetPositions.size();

	// enlarge arrays used for stencils computations
	BordStrategyMirrorSimple< Pixel, 3 > bsW( this->weight, this->nsize );
	bsW.refresh();
	BordStrategyMirrorSimple< Scalar, 3 > bsI( this->implicitH, this->nsize );
	bsI.refresh();
	BordStrategyMirrorSimple< Pixel, 3 > bsM( this->localMean, this->nsize );
	bsM.refresh();

	// LOCAL MEAN COMPUTATION:

    // assume implicitH = 1 for computing *full* image statistics
	this->localMean = nbfFastMarchingCurvature3Dmean(this->weight,this->implicitH,this->nsize);

	//nbfMatlabWriter bWriter;
	//bWriter.setFileName("w.array");
	//bWriter.write(this->localMean);

	// LOCAL VARIANCE COMPUTATION:

	// Temporarily store voxel-wise variance in: this->energy
	// this->energy = modSqr( this->weight - this->localMean );
	Array< Pixel, 3 > :: iterator iterW = this->weight.begin();
	Array< Pixel, 3 > :: iterator iterM = this->localMean.begin();
	Array< Scalar, 3 > :: iterator iterE = this->energy.begin();
	while ( iterW != this->weight.end() ){
		(*iterE) = nbfMetric::modSqr( (*iterW), (*iterM) );
		++iterW; ++iterM; ++iterE;
	}

	// enlarge for stencil computations
	BordStrategyMirrorSimple< Scalar, 3 > bsE( this->energy, this->nsize );
	bsE.refresh();

	//nbfImageWriter writer;
	//writer.setFileName("variance.vtk");
	//writer.write(this->energy);

	// compute local covariance throughout image (local average of variances at individual voxels)
	this->localVariance = nbfFastMarchingCurvature3Dmean(this->energy,this->implicitH,this->nsize);

	//bWriter.write(this->localVariance);

	// INSIDE MEAN COMPUTATION

	// set inside to 1, 0 outside
	this->implicitH = 0;
	this->implicitH[ this->aliveSetPositions ] = 1;

	// compute mean inside
	this->mean = 0;
	this->mean[ this->aliveSetPositions ] = nbfFastMarchingCurvature3Dmean(this->weight,this->implicitH,this->nsize);

	// update global statistics
	this->updateGlobalStatistics();

	// continue with regular FastMarching initialization
	nbfFastMarching< Pixel, 3, Scalar > :: initialize();
}

template< class Pixel, class Scalar >
TinyVector< int, 3 > nbfFastMarchingCurvature3D< Pixel, Scalar > :: execute( Array< Scalar, 3 > & D )
{
	// run regular FastMarching
	TinyVector< int, 3 > res = nbfFastMarching< Pixel, 3, Scalar > :: execute(D);
	
	// change output to 'state'
	this->distance = where( this->implicitH == 1, 1, 0 );

//	this->distance = this->energy;

	cout << "global minima E = " << this->energyGlobalMinima << endl;
	return res;
}

template< class Pixel, class Scalar >
void nbfFastMarchingCurvature3D< Pixel, Scalar > :: updatePoint( TinyVector< int, 3 > & current, TinyVector< int, 3 > & alive )
{
#if 1
	vector< TinyVector< int, 3 > > lCurrent;
	lCurrent.push_back(current);
	this->distance[lCurrent] = length3D(this->implicitH);
	this->distance( current ) = this->alpha * this->distance( current ) + this->weight(current);
	this->updateQueue( current );
#else
	// retrieve local (pre-computed) statistics
	Pixel     mu = this->localMean(current);
	Scalar sigma = this->localVariance(current);

	//mu = this->weight(current);

	// compute statistical energy term
	Scalar stat = log( this->sigmaInside / this->sigmaOutside ) + 
		nbfMetric::modSqr( mu, this->mean(alive) ) / this->sigmaInside -
		nbfMetric::modSqr( mu, this->meanOutside ) / this->sigmaOutside +
		sigma / this->sigmaInside -
		sigma / this->sigmaOutside;

	// stat = nbfMetric::modSqr( mu, this->mean(alive) ) - nbfMetric::modSqr( mu, this->meanOutside );

	stat = this->weight(current);
	//stat = nbfMetric::modSqr( this->weight(current), this->weight(alive) );

	// force advance when region still too small
	if ( this->foregroundSize < 1000 ){
		//stat = -10;
		stat = this->weight(current);
	}

	vector< TinyVector< int, 3 > > lCurrent;
	lCurrent.push_back(current);

	// temporarily store smooth term change
	Scalar d = this->distance(current);
	this->distance[lCurrent] = length3D(this->implicitH);

	// assign new energy change
	//this->distance( current ) = min( this->alpha * this->distance( current ) + stat, numeric_limits< Scalar > :: max() );
	this->distance( current ) = min( this->alpha * this->distance( current ) + stat, d );
	this->distance( current ) = this->alpha * this->distance( current ) + stat;

	//this->distance( current ) = this->weight(current);

	// update pqueue
    this->updateQueue( current );
#endif
}

template< class Pixel, class Scalar >
void nbfFastMarchingCurvature3D< Pixel, Scalar > :: updateGlobalStatistics()
{
	// compute background size
	// float backgroundSize = this->state.numElements() - this->foregroundSize;

	// Compute mean outside:
	Array< Scalar, 3 > :: iterator iH = this->implicitH.begin();
	Array< Scalar, 3 > :: iterator iM = this->mask.begin();
	Array<  Pixel, 3 > :: iterator iW = this->weight.begin();
	//this->meanOutside = 0;
	//while ( iH != this->implicitH.end() ){
	//	if ( (*iH) == 0 ){ // background
	//		this->meanOutside += (*iW);
	//	}
	//	++iH; ++iW;
	//}
	//this->meanOutside /= backgroundSize;

	//Scalar sum = 0;
	this->meanOutside = 0;
	this->backgroundSize = 0;
	while ( iH != this->implicitH.end() ){
		if ( ( (*iH) == 0 ) && ( (*iM) == 1 ) ){ // valid background
			//if ( sum == 0 ){
			//	this->meanOutside = (*iW);
			//}
			//else{
			//	Pixel tmp = (*iW) / sum;
			//	this->meanOutside = sum / ( sum + 1 ) * nbfMetric::sum( this->meanOutside, tmp );
			//}
			this->meanOutside = nbfMetric::sum( this->meanOutside, (*iW) );
			this->backgroundSize++;
		}
		++iH; ++iW; ++iM;
	}
	this->meanOutside /= this->backgroundSize;

	//cout << meanOutside << endl;
	//cout << this->backgroundSize << endl;
	//cout << sum( this->mask ) << endl;


	// compute variance (scalar) inside and outside simultaneously
	Scalar accInside = 0, accOutside = 0;
	Array< Scalar, 3 > :: iterator iterH = this->implicitH.begin();
	Array< Pixel, 3 > :: iterator iterW = this->weight.begin();
	Array< Pixel, 3 > :: iterator iterM = this->mean.begin();
	Array< Scalar, 3 > :: iterator iterMask = this->mask.begin();
	while ( iterH != this->implicitH.end() ){
		if ( (*iterH) == 1 ){ // inside
			accInside += nbfMetric::modSqr( (*iterW), (*iterM) );
		}
		else{ // outside
			if ( (*iterMask) == 1 ){
				accOutside += nbfMetric::modSqr( (*iterW), this->meanOutside );
			}
		}
		++iterH; ++iterW; ++iterM; ++iterMask;
	}
	this->sigmaInside = accInside / this->foregroundSize;
	this->sigmaOutside = accOutside / this->backgroundSize;
}


template< class Pixel, class Scalar >
void nbfFastMarchingCurvature3D< Pixel, Scalar > :: updateLocalStatistics( TinyVector< int, 3 > & current )
{
#if 1
	this->implicitH(current) = 1;
	this->foregroundSize++;
#else
	// differentially update sigma inside
	this->sigmaInside *= this->foregroundSize;

	// build list of affected points
	vector< TinyVector< int, 3 > > list;
	TinyVector< int, 3 > position;
	for( int i = - this->nsize; i <= this->nsize; i++ ){
		for( int j = - this->nsize; j <= this->nsize; j++ ){
			for( int k = - this->nsize; k <= this->nsize; k++ ){
				position = current + TinyVector<int,3>(i,j,k); 
				if ( this->weight.isInRange(position) ){
					if ( this->implicitH(position) == 1 ){
						// if inside
						list.push_back(position);
						// remove affected points from the inside variance estimate
						this->sigmaInside -= nbfMetric::modSqr( this->weight(position), this->mean(position) );
					}
				}
			}
		}
	}

	// change point to inside
	this->implicitH(current) = 1;
	this->foregroundSize++;

	// add current point to list
	list.push_back(current);

	// update mean inside
	this->mean[ list ] = nbfFastMarchingCurvature3Dmean( this->weight, this->implicitH, this->nsize );

	// add back affected points (and current one) to inside variance estimate
	vector< TinyVector< int, 3 > > :: iterator iter = list.begin();
	while( iter != list.end() ){
		this->sigmaInside += nbfMetric::modSqr( this->weight( *iter ), this->mean( *iter ) );
		++iter;
	}
	this->sigmaInside /= this->foregroundSize;
#endif
}

template< class Pixel, class Scalar >
void nbfFastMarchingCurvature3D< Pixel, Scalar > :: loopUpdate( TinyVector< int, 3 > & current )
{
	// Before refreshing, differentially update OUTSIDE statistics

	//float backgroundSize = this->state.numElements() - this->foregroundSize;
	this->sigmaOutside = ( this->sigmaOutside * this->backgroundSize - nbfMetric::modSqr( this->weight(current), this->meanOutside ) ) / ( this->backgroundSize - 1 );
	
	//this->meanOutside = ( this->meanOutside * backgroundSize - this->weight(current) ) / ( backgroundSize - 1 );
	// update with new metric
	this->meanOutside = this->meanOutside * this->backgroundSize;
	this->meanOutside = nbfMetric::sub( this->meanOutside, this->weight(current) ) / ( this->backgroundSize - 1 );
	this->backgroundSize--;

	// Update local statistics and promote point to inside
	this->updateLocalStatistics( current );

	// Now energy stuff
	this->totalEnergy += this->distance(current);
	this->energy( current ) = this->totalEnergy;
	this->energyVector( this->energyCounter ) = this->totalEnergy;
	this->energyCounter++;
	//cout << "E = " << this->totalEnergy << endl;
	if ( this->totalEnergy < this->energyGlobalMinima ){
		this->energyGlobalMinima = this->totalEnergy;
		//cout << "E = " << this->energyGlobalMinima << endl;
		//cout << current << endl;
		//cout << this->meanOutside << ", " << nbfMetric::mod(this->meanOutside) << ", " << sqrt(this->sigmaOutside) << ", " << sum(implicitH) <<endl;
		//cout << this->mean(current) << ", " << nbfMetric::mod(this->mean(current)) << ", " << sqrt(this->sigmaInside) << endl;
	}
	else
	{
	//	cout << "E = " << this->totalEnergy << endl;
	//	cout << current << endl;
	//	cout << this->meanOutside << ", " << nbfMetric::mod(this->meanOutside) << ", " << sqrt(this->sigmaOutside) << ", " << sum(implicitH) <<endl;
	//	cout << this->mean(current) << ", " << nbfMetric::mod(this->mean(current)) << ", " << sqrt(this->sigmaInside) << endl;
#if 0
		nbfImageWriter writer;
		writer.setFileName("implicit.vtk");
		Array< float, 3 > save( this->implicitH.shape() );
		save = this->implicitH - .5;
#endif
		//Array< float, 3 > P( save.shape() );
		//P = save;
		//BordStrategyMirrorDouble< float, 3 > bsForP( P, 1 );
		//bsForP.refresh();
		//// gaussian smoothing to compute gradient
		//for ( int i = 0; i < 1; i++ ){ // -> 5
		//	save = Laplacian3D(P);
		//	P = P + .1 * save;
		//	bsForP.refresh();
		//}
		//writer.writeFast(save);


	}

	cout << "f = " << this->foregroundSize << endl;

	if ( fmod( this->foregroundSize, 10000.0f ) == 0 ){
		//exit(0);
		nbfImageWriter writer;
		writer.setFileName("implicit.vtk");
		Array< float, 3 > save( this->implicitH.shape() );
		save = this->implicitH - .5;
		Array< float, 3 > P( save.shape() );
		P = save;
		BordStrategyMirrorDouble< float, 3 > bsForP( P, 2 );
		bsForP.refresh();
		// gaussian smoothing to compute gradient
		for ( int i = 0; i < 2; i++ ){ // -> 5
			save = Laplacian3D(P);
			P = P + .1 * save;
			bsForP.refresh();
		}
		save = P;
		writer.writeFast(save);
	}

	//nbfBlitzWriter bWriter;
	//bWriter.setFileName("w.array");
	//Array< Pixel, 3 > tmp( this->mean.shape() );
	//tmp = this->weight * this->implicitH;
	//bWriter.write(tmp);

	//int h = sum(implicitH); 
	//if ( ( h % 100 ) == 0 ){
	//	nbfImageWriter writer;
	//	writer.setFileName("implicit.vtk");
	//	writer.write(this->implicitH);

	//	//nbfBlitzWriter bWriter;
	//	//bWriter.setFileName("w.array");
	//	//bWriter.write(this->mean);
	//}
}


template< class Pixel, class Scalar >
void nbfFastMarchingCurvature3D< Pixel, Scalar > :: getNeighbors( TinyVector< int, 3 > & currentPoint,
												      vector< TinyVector< int, 3 > > & neighbors )
{
	nbfArray< Pixel, 3 > :: getFullNeighbors( currentPoint, neighbors );

	this->checkNeighbors( neighbors );

	//vector< TinyVector< int, 2 > > :: iterator iter = neighbors.begin();
	//while ( iter != neighbors.end() ){
	//	if ( this->domain(*iter) == 0 ){
	//		neighbors.erase(iter);
	//	}
	//	else{
	//		++iter;
	//	}
	//}
}

BZ_DECLARE_STENCIL_OPERATOR1(length3D,A)
	//sqrt( pow2( A(1,0,0) - A(0,0,0) ) + pow2( A(0,1,0) - A(0,0,0) ) + pow2( A(0,0,1) - A(0,0,0) ) ) +
	//sqrt( pow2( A(0,0,0) - A(-1,0,0) ) + pow2( A(-1,1,0) - A(-1,0,0) ) + pow2( A(-1,0,1) - A(-1,0,0) ) ) +
	//sqrt( pow2( A(1,-1,0) - A(0,-1,0) ) + pow2( A(0,0,0) - A(0,-1,0) ) + pow2( A(0,-1,1) - A(0,-1,0) ) ) +
	//sqrt( pow2( A(1,0,-1) - A(0,0,-1) ) + pow2( A(0,1,-1) - A(0,0,-1) ) + pow2( A(0,0,0) - A(0,0,-1) ) ); 
	// A(0,0,0) = 0
	return
	sqrt( pow2( A(1,0,0) - 1.0 ) + pow2( A(0,1,0) - 1.0 ) + pow2( A(0,0,1) - 1.0 ) )
	+ sqrt( pow2( 1.0 - A(-1,0,0) ) + pow2( A(-1,1,0) - A(-1,0,0) + 0.0 ) + pow2( A(-1,0,1) - A(-1,0,0) + 0.0 ) )
	+ sqrt( pow2( A(1,-1,0) - A(0,-1,0) + 0.0 ) + pow2( 1.0 - A(0,-1,0) ) + pow2( A(0,-1,1) - A(0,-1,0) + 0.0 ) )
	+ sqrt( pow2( A(1,0,-1) - A(0,0,-1) + 0.0 ) + pow2( A(0,1,-1) - A(0,0,-1) + 0.0 ) + pow2( 1.0 - A(0,0,-1) + 0.0 ) )
	- sqrt( pow2( A(1,0,0) - 0.0 ) + pow2( A(0,1,0) - 0.0 ) + pow2( A(0,0,1) - 0.0 ) )
    - sqrt( pow2( 0.0 - A(-1,0,0) ) + pow2( A(-1,1,0) - A(-1,0,0) + 0.0 ) + pow2( A(-1,0,1) - A(-1,0,0) + 0.0 ) )
	- sqrt( pow2( A(1,-1,0) - A(0,-1,0) + 0.0 ) + pow2( 0.0 - A(0,-1,0) ) + pow2( A(0,-1,1) - A(0,-1,0) + 0.0 ) )
	- sqrt( pow2( A(1,0,-1) - A(0,0,-1) + 0.0 ) + pow2( A(0,1,-1) - A(0,0,-1) + 0.0 ) + pow2( 0.0 - A(0,0,-1) ) ); 
BZ_END_STENCIL_OPERATOR

// allow direct call of stencil operator
BZ_ET_STENCIL(length3D,P_numtype)

// generic mean

template<class T1, class T2>
inline _bz_typename T1::T_numtype
nbfFastMarchingCurvature3Dmean( T1& A, T2& S, int nsize) {
	float sum = 0;
	T1::T_numtype res = 0;
	T1::T_numtype tmp;
	for ( int i = -nsize; i <= nsize; i++ ){								
		for ( int j = -nsize; j <= nsize; j++ ){							
			for ( int k = -nsize; k <= nsize; k++ ){
				if ( S(i,j,k) == 1 ){
					//res += A(i,j,k);
					//sum += 1.0;
					//tmp = A(i,j,k) / sum;
					//res = sum / ( sum + 1 ) * nbfMetric::sum( res, tmp );
					res = nbfMetric::sum( res, A(i,j,k) );
					sum++;
				}
			}														
		}															
	}
	res = res / sum;
	return res;
}

BZ_ET_STENCILV2(nbfFastMarchingCurvature3Dmean,3)


#endif /* FILE_nbfFastMarchingCurvature3D */