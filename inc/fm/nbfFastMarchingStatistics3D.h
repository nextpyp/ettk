#ifndef FILE_nbfFastMarchingStatistics3D
#define FILE_nbfFastMarchingStatistics3D

// Class nbfFastMarchingStatistics.
//
// Implements 3D fast marching method with arbitrary weights.
// 

#include <nbf-stencil-et.h>
#include <fm/nbfFastMarching3D.h>
#include <bs/nbfBordStrategyMirror.h>
#include <nbfDifferentials.h>

#include <io/nbfImageWriter.h>

template< class Pixel, class Scalar = Pixel >
class nbfFastMarchingStatistics3D : public nbfFastMarching3D< Pixel, Scalar >
{
public:

	// constructor takes weight array as input
	nbfFastMarchingStatistics3D( Array< Pixel, 3 > & );
	~nbfFastMarchingStatistics3D(){};

	// override from FastMarching 
	// (this is done only to return the 'state' as the output and not the 'distance' values)
	virtual TinyVector< int, 3 > execute( Array< Scalar, 3 > & );

	// set the smoothness term weight (defaults to .25)
	void setAlpha( Scalar a ){ this->alpha = a; }

	// set mean neighborhood size (defaults to 5)
	void setMeanNeighborhoodSize( int s ){ this->msize = s; }
	void setVarianceNeighborhoodSize( int s ){ this->ssize = s; }

	void getNeighbors( TinyVector< int, 3 > &, vector< TinyVector< int, 3 > > & );

protected:

	// override from FastMarching
	void updatePoint( TinyVector< int, 3 > &, TinyVector< int, 3 > & );

	// override from FastMarching
	virtual void initialize();

	// update mean and variance in the neighborhood around a point
	void updateLocalStatistics( TinyVector< int, 3 > & );

	// override from FastMarching
	void loopUpdate( TinyVector< int, 3 > & );

	//
	void updateSigmaOutside( TinyVector< int, 3 > & );

	// Statistical members

	// neighborhood size for mean computation
	int msize;

	// neighborhood size for variance computation
	int ssize;

	// Mean: inside (piecewise smooth), outside (constant)
	Array< Pixel, 3 > mean;

	// Mean outside ( note that mean(outside) = meanOutside )
	Pixel meanOutside;

	// variance inside/outside
	Scalar sigmaInside;
	Scalar sigmaInsideWeights;
	Array< Scalar, 3 > sigmaOutside;

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
nbfFastMarchingStatistics3D< Pixel, Scalar > :: nbfFastMarchingStatistics3D( Array< Pixel, 3 > & weight )
: nbfFastMarching3D< Pixel, Scalar >( weight )
{
	this->msize = 1;
	this->ssize = 5;

	this->mean.resize( weight.shape() );
	this->localMean.resize( weight.shape() );
	this->localVariance.resize( weight.shape() );

	this->sigmaOutside.resize( weight.shape() );

	// stopping criteria
	this->energy.resize( weight.shape() );
	this->totalEnergy = 0;
	this->energyVector.resize( weight.numElements() );
	this->energyCounter = 0;
	this->energyGlobalMinima = numeric_limits<Scalar>::max();

	this->implicitH.resize( weight.shape() );
	this->implicitH = 1;
	this->alpha = .25;
	this->foregroundSize = 0;
	this->previousSizeUpdate = this->foregroundSize;

	// THIS DOES NOT WORK FOR THE SCALAR CASE
	//// initialize background mask and size
	this->mask.resize( weight.shape() );
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
void nbfFastMarchingStatistics3D< Pixel, Scalar > :: initialize()
{
	// starting foreground size is number of points in initial alive set
	this->foregroundSize = this->aliveSetPositions.size();

	// enlarge arrays used for stencils computations
	BordStrategyMirrorSimple< Pixel, 3 > bsW( this->weight, this->msize );
	bsW.refresh();
	BordStrategyMirrorSimple< Scalar, 3 > bsI( this->implicitH, this->msize );
	bsI.refresh();
	BordStrategyMirrorSimple< Pixel, 3 > bsM( this->localMean, this->msize );
	bsM.refresh();

	// LOCAL MEAN COMPUTATION:

    // assume implicitH = 1 for computing *full* image statistics
	 this->localMean = nbfFastMarchingStatistics3Dmean(this->weight,this->implicitH,this->nsize);

	//nbfBlitzWriter bWriter;
	//bWriter.setFileName("w.array");
	//bWriter.write(this->weight);

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
	BordStrategyMirrorSimple< Scalar, 3 > bsE( this->energy, this->ssize );
	bsE.refresh();

	//nbfMatlabWriter writer;
	//writer.setFileName("ipath");
	//writer.write(this->energy);

	 // compute local covariance throughout image (local average of variances at individual voxels)
	 this->localVariance = nbfFastMarchingStatistics3Dmean(this->energy,this->implicitH,this->nsize);

	//writer.write(this->localVariance);

	// INSIDE MEAN COMPUTATION

	// set inside to 1, 0 outside
	this->implicitH = 0;
	this->implicitH[ this->aliveSetPositions ] = 1;
	writer.write(implicitH);

	// compute mean inside
	// this->mean = 0;
	this->mean[ this->aliveSetPositions ] = nbfFastMarchingStatistics3Dmean(this->weight,this->implicitH,this->msize);

	// compute variance (scalar) inside
	Scalar accInside = 0;

	// initialize sigma inside
	this->sigmaInsideWeights = 0;
	vector< TinyVector< int, 3 > > :: iterator iAlives = this->aliveSetPositions.begin();
	while ( iAlives != this->aliveSetPositions.end() ){
		// variance inside
		Scalar sigmaWeight = 1 - sqrt( nbfMetric::modSqr( this->weight(*iAlives), this->mean(*iAlives) ) );
		
		// SCALAR - OVERRIDE TO UNIFORM WEIGHTING
		sigmaWeight = 1;

		accInside += ( sigmaWeight * nbfMetric::modSqr( this->weight(*iAlives), this->mean(*iAlives) ) );
		this->sigmaInsideWeights += sigmaWeight;
		// variance outside
        //this->updateSigmaOutside(*iAlives);
		++iAlives;
	}

	// sigma inside
	this->sigmaInside = accInside / this->foregroundSize / this->sigmaInsideWeights;

	// continue with regular FastMarching initialization
	nbfFastMarching< Pixel, 3, Scalar > :: initialize();
}

template< class Pixel, class Scalar >
TinyVector< int, 3 > nbfFastMarchingStatistics3D< Pixel, Scalar > :: execute( Array< Scalar, 3 > & D )
{
	// run regular FastMarching
	TinyVector< int, 3 > res = nbfFastMarching< Pixel, 3, Scalar > :: execute(D);
	
	// change output to 'state'
	this->distance = where( this->implicitH == 1, 1, 0 );

//	this->distance = this->energy;

	//cout << "global minima E = " << this->energyGlobalMinima << endl;
	return res;
}

template< class Pixel, class Scalar >
void nbfFastMarchingStatistics3D< Pixel, Scalar > :: updatePoint( TinyVector< int, 3 > & current, TinyVector< int, 3 > & alive )
{
	// if first time, compute mean and variance outside
	if ( this->state(current) == FAR_SET ){
		vector< TinyVector< int, 3 > > list;
		list.push_back( current );
		this->mean[ list ] = nbfFastMarchingStatistics3Dmean( this->weight, this->implicitH, this->msize );
		this->updateSigmaOutside( current );
	}

	// retrieve local (pre-computed) statistics
	Pixel     mu = this->localMean(current);
	Scalar sigma = this->localVariance(current);

	// compute statistical energy term
	//cout << this->weight(current) << endl;
	//cout << this->mean(current) << endl;

	// TENSOR
	//Scalar projectionSqr = nbfMetric::modSqr( this->weight(current), this->mean(current) );
	Scalar currentSigmaOutside = this->sigmaOutside(current);
	//Scalar stat = log( this->sigmaInside / currentSigmaOutside ) + 
	//	projectionSqr / this->sigmaInside -
	//	pow2( 1 - sqrt(projectionSqr) ) / currentSigmaOutside;
	
	// SCALAR
	Scalar stat = log( this->sigmaInside / currentSigmaOutside ) + 
		nbfMetric::modSqr( mu, this->mean(alive) ) / this->sigmaInside -
		nbfMetric::modSqr( mu, this->meanOutside ) / currentSigmaOutside +
		sigma / this->sigmaInside -
		sigma / currentSigmaOutside;

	//// force advance when region still too small
	//if ( this->foregroundSize < 10 ){
	//	// stat = -1;
	//	stat = nbfMetric::mod( this->weight(current), this->weight(alive) );
	//}

	// cout << current << ", " << projectionSqr << ", " << currentSigmaOutside << ", " << stat << endl;

	vector< TinyVector< int, 3 > > lCurrent;
	lCurrent.push_back(current);

	// temporarily store smooth term change
	Scalar d = this->distance(current);
	this->distance[lCurrent] = length3D(this->implicitH);
	this->distance( current ) = this->alpha * this->distance( current ) + stat;

	Scalar tmp;

	// new
	if ( this->foregroundSize > 2 ){
		this->distance[lCurrent] = nbfFastMarchingStatistics3Darea(this->implicitH);
		tmp = distance(current);
		if ( this->distance(current) == -1 ){
			this->distance(current) = numeric_limits<Scalar>::max();
		}
		if ( this->distance(current) == 1 ){
			if ( stat < 0 ){
				this->distance( current ) = stat;
			}
			else{
				this->distance(current) = - 1 / stat;
			}
		}
		if ( this->distance(current) == 0 ){
			this->distance( current ) = stat;
		}
	}

	// assign new energy change
	//this->distance( current ) = min( this->alpha * this->distance( current ) + stat, numeric_limits< Scalar > :: max() );
	//this->distance( current ) = min( this->alpha * this->distance( current ) + stat, d );

	// update pqueue
    this->updateQueue( current );
}


template< class Pixel, class Scalar >
void nbfFastMarchingStatistics3D< Pixel, Scalar > :: updateLocalStatistics( TinyVector< int, 3 > & current )
{
	// differentially update sigma inside
	this->sigmaInside *= this->foregroundSize;
	this->sigmaInside *= this->sigmaInsideWeights;

	// build list of affected points
	vector< TinyVector< int, 3 > > inside;
	vector< TinyVector< int, 3 > > trial;
	TinyVector< int, 3 > position;
	for( int i = - this->msize; i <= this->msize; i++ ){
		for( int j = - this->msize; j <= this->msize; j++ ){
			for( int k = - this->msize; k <= this->msize; k++ ){
				position = current + TinyVector<int,3>(i,j,k); 

				if ( this->weight.isInRange(position) ){
					// if inside
					if ( this->implicitH(position) == 1 ){
						inside.push_back(position);
						// remove affected points from the inside variance estimate
						Scalar sigmaWeight = 1 - sqrt( nbfMetric::modSqr( this->weight(position), this->mean(position) ) );
						this->sigmaInside -= ( sigmaWeight * nbfMetric::modSqr( this->weight(position), this->mean(position) ) );
						this->sigmaInsideWeights -= sigmaWeight;
					}
					// store trial points
					if ( this->state(position) == TRIAL_SET ){
						trial.push_back(position);
					}
				}
			}
		}
	}

	// change point to inside
	this->implicitH(current) = 1;
	this->foregroundSize++;

	// add current point to list
	inside.push_back(current);

	// update mean inside
	this->mean[ inside ] = nbfFastMarchingStatistics3Dmean( this->weight, this->implicitH, this->msize );

	// add back affected points (and current one) to inside variance estimate
	vector< TinyVector< int, 3 > > :: iterator iter = inside.begin();
	while( iter != inside.end() ){
		Scalar sigmaWeight = 1 - sqrt( nbfMetric::modSqr( this->weight( *iter ), this->mean( *iter ) ) );
		this->sigmaInside += ( sigmaWeight * nbfMetric::modSqr( this->weight( *iter ), this->mean( *iter ) ) );
		this->sigmaInsideWeights += sigmaWeight;
		++iter;
	}
	this->sigmaInside /= this->foregroundSize;
	this->sigmaInside /= this->sigmaInsideWeights;

	// update mean and sigma at trial points
	this->mean[ trial ] = nbfFastMarchingStatistics3Dmean( this->weight, this->implicitH, this->msize );
	for( iter = trial.begin(); iter != trial.end(); ++iter ){
		this->updateSigmaOutside(*iter);
	}

	// update sigma outside
	//this->updateSigmaOutside(current);
}


template< class Pixel, class Scalar >
void nbfFastMarchingStatistics3D< Pixel, Scalar > :: updateSigmaOutside( TinyVector< int, 3 > & current )
{
	Range I( max( current[firstDim] - this->ssize, this->weight.lbound(firstDim) ), 
		     min( current[firstDim] + this->ssize, this->weight.ubound(firstDim) ) );
	Range J( max( current[secondDim] - this->ssize, this->weight.lbound(secondDim) ), 
		     min( current[secondDim] + this->ssize, this->weight.ubound(secondDim) ) );
	Range K( max( current[thirdDim] - this->ssize, this->weight.lbound(thirdDim) ), 
		     min( current[thirdDim] + this->ssize, this->weight.ubound(thirdDim) ) );

	Array<  Pixel, 3 > viewW( this->weight(I,J,K) );
	Array< Scalar, 3 > viewI( this->implicitH(I,J,K) );
	Array< Scalar, 3 > viewM( this->mask(I,J,K) );

	Array<  Pixel, 3 > :: iterator iterW = viewW.begin();
	Array< Scalar, 3 > :: iterator iterI = viewI.begin();
	Array< Scalar, 3 > :: iterator iterM = viewM.begin();

	Scalar count = 0;
	Scalar sigma = 0;
	Scalar wsum = 0;

	while ( iterW != viewW.end() ){
		if ( ( (*iterI) == 0 ) & ( (*iterM) == 1 ) ){ // if valid background
			Scalar w = sqrt( nbfMetric::modSqr( (*iterW), this->mean(current) ) );
			wsum += w;
			sigma = sigma + w * pow2( 1 - sqrt( nbfMetric::modSqr( (*iterW), this->mean(current) ) ) );
			count++;
		}
		++iterW; ++iterI; ++iterM;
	}

	this->sigmaOutside( current ) = sigma / count / wsum;
}

template< class Pixel, class Scalar >
void nbfFastMarchingStatistics3D< Pixel, Scalar > :: loopUpdate( TinyVector< int, 3 > & current )
{
	//cout << current << endl;
	Scalar projectionSqr = nbfMetric::modSqr( this->weight(current), this->mean(current) );

	// Update local statistics and promote point to inside
	this->updateLocalStatistics( current );

	// Now energy stuff
	this->totalEnergy += this->distance(current);
	this->energy( current ) = this->totalEnergy;
	this->energyVector( this->energyCounter ) = this->totalEnergy;
	this->energyCounter++;
	//cout << "E = " << this->totalEnergy << ", sigma = " << sqrt(this->sigmaInside) << ", outside = " << sqrt(this->sigmaOutside(current)) << ", d = " << sqrt(projectionSqr) << endl;
	if ( this->totalEnergy <= this->energyGlobalMinima ){
		this->energyGlobalMinima = this->totalEnergy;
		//cout << "E = " << this->energyGlobalMinima << endl;
		//cout << current << endl;
		//cout << this->meanOutside << ", " << nbfMetric::mod(this->meanOutside) << ", " << sqrt(this->sigmaOutside) << ", " << sum(implicitH) <<endl;
		//cout << this->mean(current) << ", " << nbfMetric::mod(this->mean(current)) << ", " << sqrt(this->sigmaInside) << endl;
	}
	else
	{
		this->flagForceStop = true;

		//cout << "E = " << this->totalEnergy << endl;

		Array< Scalar, 3 > phi( this->implicitH.shape() );
		phi = where( this->implicitH == 1, 1, -1);
		BordStrategyMirrorSimple< Scalar, 3 > bsForPhi( phi, 1 );
		bsForPhi.refresh();

		//Array< Scalar, 3 > F( phi.shape() );
		//F = 1.0;

		// parameters
		float timeStep = .25;	// evolution time step
		int iterations = 3; //3;		// evolution iterations
		for ( int i = 0; i < iterations; i++ ){
			this->implicitH = 1.0 * curvature3D(phi) + 1.0 * nablaPlus3D(phi);
			phi = phi + timeStep * this->implicitH;
			bsForPhi.refresh();
		}

		nbfImageWriter writer;
		writer.setFileName("implicit_current.vtk");
		writer.write(phi);


		//Array< Scalar, 3 > D( this->implicitH.shape() );
		//Array< Scalar, 3 > A( this->implicitH.shape() );
		//Array< Scalar, 3 > I( this->implicitH.shape() );

		//// compute iso-distances to path first
		//A = 1;
		//nbfFastMarching3D< Scalar > fm3d(A);

		//vector< TinyVector< int, 3 > > vpositions;
		//vector< Scalar > vdistances;
		//Array< Scalar, 3 > :: iterator iter = this->implicitH.begin();
		//while ( iter != this->implicitH.end() ){
		//	if ( (*iter) == 1 ){
		//		vpositions.push_back( iter.position() );
		//		vdistances.push_back( 0 );
		//	}
		//	++iter;
		//}
		//fm3d.setAliveSet(vpositions,vdistances);
		//fm3d.setStopDistance(5);
		//fm3d.execute(D);
		//I = D;
		//cout << min(I) << ", " << max(I) << endl;

		//vpositions.clear(); vdistances.clear();
		//Array< Scalar, 3 > :: iterator iter2 = this->implicitH.begin();
		//while ( iter2 != this->implicitH.end() ){
		//	if ( (*iter2) == 0 ){
		//		vpositions.push_back( iter2.position() );
		//		vdistances.push_back( 0 );
		//	}
		//	++iter2;
		//}

		//nbfFastMarching3D< Scalar > fm3d2(A);
		//fm3d2.setAliveSet(vpositions,vdistances);
		//fm3d2.setStopDistance(5);
		//fm3d2.execute(D);
		//cout << min(D) << ", " << max(D) << endl;

		//I = where( this->implicitH == 1, -D, I );
		//cout << min(I) << ", " << max(I) << endl;

		//writer.write(I);

		//	cout << current << endl;
	//	cout << this->meanOutside << ", " << nbfMetric::mod(this->meanOutside) << ", " << sqrt(this->sigmaOutside) << ", " << sum(implicitH) <<endl;
	//	cout << this->mean(current) << ", " << nbfMetric::mod(this->mean(current)) << ", " << sqrt(this->sigmaInside) << endl;
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

	// update all points in narrow band
	for ( int i = 0; i < this->pqueue->GetNumberOfItems(); i++ ){
		int id = this->pqueue->Peek(i);
		TinyVector< int, 3 > currentPosition = this->int2array(id);
		this->updatePoint( currentPosition, currentPosition );
	}
}


template< class Pixel, class Scalar >
void nbfFastMarchingStatistics3D< Pixel, Scalar > :: getNeighbors( TinyVector< int, 3 > & currentPoint,
												      vector< TinyVector< int, 3 > > & neighbors )
{
	nbfArray< Pixel, 3 > :: getFullNeighbors( currentPoint, neighbors );

	this->checkNeighbors( neighbors );
}


template< class Pixel, class Scalar >
void nbfFastMarchingStatistics3D< Pixel, Scalar > :: updateGlobalStatistics()
{
	// compute background size
	// float backgroundSize = this->state.numElements() - this->foregroundSize;

	// Compute mean outside:
	Array< Scalar, 3 > :: iterator iH = this->implicitH.begin();
	Array< Scalar, 3 > :: iterator iM = this->mask.begin();
	Array<  Pixel, 3 > :: iterator iW = this->weight.begin();

	//Scalar sum = 0;
	this->meanOutside = 0;
	this->backgroundSize = 0;
	while ( iH != this->implicitH.end() ){
		if ( ( (*iH) == 0 ) & ( (*iM) == 1 ) ){ // valid background
			this->meanOutside = nbfMetric::sum( this->meanOutside, (*iW) );
			this->backgroundSize++;
		}
		++iH; ++iW; ++iM;
	}
	this->meanOutside /= this->backgroundSize;

	//// compute variance (scalar) inside and outside simultaneously
	//Scalar accInside = 0, accOutside = 0;
	//Array< Scalar, 3 > :: iterator iterH = this->implicitH.begin();
	//Array< Pixel, 3 > :: iterator iterW = this->weight.begin();
	//Array< Pixel, 3 > :: iterator iterM = this->mean.begin();
	//Array< Scalar, 3 > :: iterator iterMask = this->mask.begin();
	//while ( iterH != this->implicitH.end() ){
	//	if ( (*iterH) == 1 ){ // inside
	//		accInside += nbfMetric::modSqr( (*iterW), (*iterM) );
	//	}
	//	else{ // outside
	//		if ( (*iterMask) == 1 ){
	//			accOutside += nbfMetric::modSqr( (*iterW), this->meanOutside );
	//		}
	//	}
	//	++iterH; ++iterW; ++iterM; ++iterMask;
	//}
	//this->sigmaInside = accInside / this->foregroundSize;
	//this->sigmaOutside = accOutside / this->backgroundSize;
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


template<class T1>
inline _bz_typename T1::T_numtype
nbfFastMarchingStatistics3Darea( T1& A ) 
{
	// if only one alive neighbor, halt growing
	if ( ( A(0,0,1) + A(0,0,-1) + A(0,1,0) + A(0,-1,0) + A(0,1,1) + A(0,1,-1) + A(0,-1,1) + A(0,-1,-1) +
		   A(1,0,1) + A(1,0,-1) + A(1,1,0) + A(1,-1,0) + A(1,1,1) + A(1,1,-1) + A(1,-1,1) + A(1,-1,-1) + A(1,0,0) +
		   A(-1,0,1) + A(-1,0,-1) + A(-1,1,0) + A(-1,-1,0) + A(-1,1,1) + A(-1,1,-1) + A(-1,-1,1) + A(-1,-1,-1) + A(-1,0,0) ) < 3 ){
			 return -1;
		   }
	
		   // check for each coordinate

		   // x

		   T1::T_numtype a = A(1,0,0) + A(-1,0,0);
		   T1::T_numtype b = A(0,1,0) + A(0,-1,0);
		   T1::T_numtype c = A(1,1,0) + A(-1,-1,0);
		   T1::T_numtype d = A(1,-1,0) + A(-1,1,0);

		   if ( ( a + b + c + d > 6 ) |
			    ( ( a == 2 ) & ( b == 0 ) ) |
			    ( ( a == 0 ) & ( b == 2 ) ) |
				( ( a == 0 ) & ( b == 0 ) & ( ( c == 2 ) | ( d == 2 ) ) ) ){
					return 1;
				}

		   if ( ( ( A(-1,0,0) == 1 ) & ( A(1,0,0) == 0 ) & ( A(0,1,0) == 0 ) & ( A(0,-1,0) == 0 ) &
			      ( A(1,1,0) + A(1,-1,0) > 0 ) ) |
				( ( A(-1,0,0) == 0 ) & ( A(1,0,0) == 1 ) & ( A(0,1,0) == 0 ) & ( A(0,-1,0) == 0 ) &
				  ( A(-1,1,0) + A(-1,-1,0) > 0 ) ) |
				( ( A(0,-1,0) == 1 ) & ( A(0,1,0) == 0 ) & ( A(1,0,0) == 0 ) & ( A(-1,0,0) == 0 ) &
				  ( A(1,1,0) + A(-1,1,0) > 0 ) ) |
				( ( A(0,-1,0) == 0 ) & ( A(0,1,0) == 1 ) & ( A(1,0,0) == 0 ) & ( A(-1,0,0) == 0 ) &
				  ( A(1,-1,0) + A(-1,-1,0) > 0 ) ) ){
					return 1;
				}


		   // y

		   a = A(0,1,0) + A(0,-1,0);
		   b = A(0,0,1) + A(0,0,-1);
		   c = A(0,1,1) + A(0,-1,-1);
		   d = A(0,1,-1) + A(0,-1,1);

		   if ( ( a + b + c + d > 6 ) |
			    ( ( a == 2 ) & ( b == 0 ) ) |
			    ( ( a == 0 ) & ( b == 2 ) ) |
				( ( a == 0 ) & ( b == 0 ) & ( ( c == 2 ) | ( d == 2 ) ) ) ){
					return 1;
				}

		   if ( ( ( A(0,-1,0) == 1 ) & ( A(0,1,0) == 0 ) & ( A(0,0,1) == 0 ) & ( A(0,0,-1) == 0 ) &
			      ( A(0,1,1) + A(0,1,-1) > 0 ) ) |
				( ( A(0,-1,0) == 0 ) & ( A(0,1,0) == 1 ) & ( A(0,0,1) == 0 ) & ( A(0,0,-1) == 0 ) &
				  ( A(0,-1,1) + A(0,-1,-1) > 0 ) ) |
				( ( A(0,0,-1) == 1 ) & ( A(0,0,1) == 0 ) & ( A(0,1,0) == 0 ) & ( A(0,-1,0) == 0 ) &
				  ( A(0,1,1) + A(0,-1,1) > 0 ) ) |
				( ( A(0,0,-1) == 0 ) & ( A(0,0,1) == 1 ) & ( A(0,1,0) == 0 ) & ( A(0,-1,0) == 0 ) &
				( A(0,1,-1) + A(0,-1,-1) > 0 ) ) ){
					return 1;
				}

		   // z
		   a = A(1,0,0) + A(-1,0,0);
		   b = A(0,0,1) + A(0,0,-1);
		   c = A(1,0,1) + A(-1,0,-1);
		   d = A(1,0,-1) + A(-1,0,1);

		   if ( ( a + b + c + d > 6 ) |
			    ( ( a == 2 ) & ( b == 0 ) ) |
			    ( ( a == 0 ) & ( b == 2 ) ) |
				( ( a == 0 ) & ( b == 0 ) & ( ( c == 2 ) | ( d == 2 ) ) ) ){
					return 1;
				}

		   if ( ( ( A(-1,0,0) == 1 ) & ( A(1,0,0) == 0 ) & ( A(0,0,1) == 0 ) & ( A(0,0,-1) == 0 ) &
			      ( A(1,0,1) + A(1,0,-1) > 0 ) ) |
				( ( A(-1,0,0) == 0 ) & ( A(1,0,0) == 1 ) & ( A(0,0,1) == 0 ) & ( A(0,0,-1) == 0 ) &
				  ( A(-1,0,1) + A(-1,0,-1) > 0 ) ) |
				( ( A(0,0,-1) == 1 ) & ( A(0,0,1) == 0 ) & ( A(1,0,0) == 0 ) & ( A(-1,0,0) == 0 ) &
				  ( A(1,0,1) + A(-1,0,1) > 0 ) ) |
				( ( A(0,0,-1) == 0 ) & ( A(0,0,1) == 1 ) & ( A(1,0,0) == 0 ) & ( A(-1,0,0) == 0 ) &
				( A(1,0,-1) + A(-1,0,-1) > 0 ) ) ){
					return 1;
				}
		return 0;
}


BZ_ET_STENCIL(nbfFastMarchingStatistics3Darea,P_numtype)

// generic mean

template<class T1, class T2>
inline _bz_typename T1::T_numtype
nbfFastMarchingStatistics3Dmean( T1& A, T2& S, int nsize) {
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

BZ_ET_STENCILV2(nbfFastMarchingStatistics3Dmean,3)

#endif /* FILE_nbfFastMarchingStatistics3D */