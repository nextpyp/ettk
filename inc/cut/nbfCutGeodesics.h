#ifndef FILE_nbfCutGeodesics
#define FILE_nbfCutGeodesics

// Class nbfCutGeodesics.
//
// Find geodesic paths in 2D images.
// 

#include <fm/nbfFastMarchingFool.h>
#include <fm/nbfFastMarching2D8.h>
#include <bp/nbfGeodesicPath.h>

#include <vector>
#include <algorithm>

template< class Pixel >
class nbfCutGeodesics : public nbfGeodesicPath< Pixel >
{
public:

	// takes $g$ array and center point as input
	nbfCutGeodesics( Array< Pixel, 2 > &, TinyVector< int, 2 > & );
	nbfCutGeodesics( Array< Pixel, 2 > & );

	~nbfCutGeodesics(){};

	// get minimal circular path in implicit form and distance array
	void getCircularPath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, Array< Pixel, 2 > &, bool = true );
	void getCircularPath( Array< Pixel, 2 > &, Array< Pixel, 2 > & );
	void getNewCircularPath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, Array< Pixel, 2 > &, bool = true );
	void getNewCircularPath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, vector< TinyVector< int, 2 > > &, bool = true );

	// Once we computed the 3D distances, we look at the accumulated distances in the
	// slice on the other side of the cut. In the polar case, we pick the point with the minimum
	// distance and then start going 'horizontally' in both directions.
	// But on the cartesian domain, we need to go 'circularly' and this is not so simple
	// to program, so (restricted to the mentioned slice) we compute at once a geodesic
	// that should correspond to the boundary of the shape on this slice.
	//
	// This function then computes the minimal semi-circular path surrounding the center point.
	// It is semicircular because the center point is in the upper boundary of the image.
	// 
	// To make it robust, a forward and a backward geodesic are considered.
	// (this is because they can give quite different geodesics), we compute both and keep the one
	// with the minimum lenght.
	void getSemiCircularPath( Array< Pixel, 2 > &, vector< TinyVector< int, 2 > > & );
	void getConstrainedSemiCircularPath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, TinyVector< int, 2 > &, TinyVector< int, 2 > &, vector< TinyVector< int, 2 > > & );
	void getFullCircularPath( Array< Pixel, 2 > &, vector< TinyVector< Pixel, 2 > > &, bool = true, Pixel = 1e-1 );
	void getConstrainedCircularPath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, TinyVector< int, 2 > &, vector< TinyVector< int, 2 > > & );

	void curvature( Array< Pixel, 2 > &, int = 10 );

protected:

	TinyVector< int, 2 > getMinimaOnTheBoundary( Array< Pixel, 2 > & );

};


template< class Pixel >
nbfCutGeodesics< Pixel > :: nbfCutGeodesics( Array< Pixel, 2 > & weight )
: nbfGeodesicPath< Pixel >( weight )
{
}

template< class Pixel >
nbfCutGeodesics< Pixel > :: nbfCutGeodesics( Array< Pixel, 2 > & weight,
										     TinyVector< int, 2 > & center )
											 : nbfGeodesicPath< Pixel >( weight, center )
{
	this->updateTheta();
}


template< class Pixel >
void nbfCutGeodesics< Pixel > :: getCircularPath( Array< Pixel, 2 > & weights,
													 TinyVector< int, 2 > & start,
													 Array< Pixel, 2 > & ipath,
													 bool closed = true )
{
	this->weight.reference( weights );

	// find a circular path between cut-sides of a 2d image.
	// The circular constraint is enforced by restricting the start and
	// end points of the path to be circular neighbors, i.e. to have the
	// same coordinate in the appropriate dimension.

	// back propagation from start until we reach a point with 0 distance
	vector< TinyVector< int, 2 > > path;
	this->getForwardPath( weights, start, path );

	this->getImplicitPath( path, ipath );

	nbfMatlabWriter writer;
	writer.setFileName("ipath");
	writer.write(ipath);

	// if last point in the path "coincides" with start, we are done
	// else, we need to find a path from a different starting point
	if ( closed & ( abs( path[ path.size() - 1 ](firstDim) - start(firstDim) ) > 1 ) ){
		
		// re-set start point
		int newStart = path[ path.size() - 1 ](firstDim);

		// the new starting point should be +/- 1 pixel closer to the computed end position
		// either we pick the closest one to the original starting point
		// or we take the exact same position as the end of the computed path.
#if 0
		if ( start(firstDim) > newStart ){
			start(firstDim) = newStart + 1;
		} 
		else {
			if ( start(firstDim) < newStart ){
				start(firstDim) = newStart - 1;
			}
			else {
				start(firstDim) = newStart;
			}
		}
#else
		start(firstDim) = newStart;
#endif

		cout << path[ path.size() - 1 ] << endl;
		cout << start << endl;

#ifdef NBF_DEBUG
		cout << "re-start = " << start << endl;
#endif
		// re-compute path
		this->getForwardPath( weights, start, path );
	}

	// get implicit representation
	this->getImplicitPath( path, ipath );
}


template< class Pixel >
void nbfCutGeodesics< Pixel > :: getCircularPath( Array< Pixel, 2 > & weights, 
													 Array< Pixel, 2 > & ipath )
{
	TinyVector< int, 2 > start;

	// look for starting point with the smallest distance value
	Array< Pixel, 1 > line( weights( Range( this->center[firstDim] + 1, toEnd ), 
		                                    this->center[secondDim] + 1 ) );

	// if stating point is close to the center coordinates
	if ( ( min(line) == numeric_limits<Pixel>::max() ) ) // ||
//		( start(firstDim) - this->center[firstDim] < 0 ) )
	{
		ipath = numeric_limits<Pixel>::max();
#ifdef NBF_DEBUG
		cout << "bad start = " << start << endl;
#endif

	}
	else
	{
		start(firstDim) = ( minIndex( line ) )(firstDim) + this->center[firstDim] + 1;
		start(secondDim) = this->center[secondDim] + 1;

		cout << start << endl;

		this->getCircularPath( weights, start, ipath );
	}
}


template< class Pixel >
void nbfCutGeodesics< Pixel > :: getNewCircularPath( Array< Pixel, 2 > & weights, 
													 TinyVector< int, 2 > & start,
													 vector< TinyVector< int, 2 >  > & path,
													 bool polar )
{
	Array< Pixel, 2 > w( weights.shape() );
	w = weights - min(weights);
	if ( polar == true ){
		firstIndex i;
		secondIndex j;
		w = w / sqrt( pow2( center[0] - i + 0.0 ) + pow2( center[1] - j + 0.0 ) );
		w = w + 1e-10;
	}
	else{
        w = w + 1;
	}

	nbfFastMarching2D8< Pixel > fastMarching( w, center );

	TinyVector< int, 2 > aliveP( start[firstDim], start[secondDim] - 1 );
	Pixel aliveD = 0;
	fastMarching.setAliveSet( aliveP, aliveD );
	Array< Pixel, 2 > d( w.shape() );

	fastMarching.setStopPoint( start );
	fastMarching.execute( d );

	//nbfMatlabWriter writer;
	//writer.setFileName("ipath");
	//writer.write(d);

	// end point
	this->getPath( d, start, path );
	//this->weight = d + 1e-10;
	//vector< TinyVector< int, 2 > > path;
	//this->getForwardPath( d, start, path );
	//this->getImplicitPath( path, ipath );
}

template< class Pixel >
void nbfCutGeodesics< Pixel > :: getNewCircularPath( Array< Pixel, 2 > & weights, 
													 TinyVector< int, 2 > & start,
													 Array< Pixel, 2 > & ipath,
													 bool polar )
{
	vector< TinyVector< int, 2 > > path;
	this->getNewCircularPath( weights, start, path, polar );
	this->getImplicitPath( path, ipath );
}


template< class Pixel >
TinyVector< int, 2 > nbfCutGeodesics< Pixel> :: getMinimaOnTheBoundary( Array< Pixel, 2 > & input )
{
	// look for minima in four image sides
	
	Array< Pixel, 1 > lineUp = input( 0, Range::all() );
	Array< Pixel, 1 > lineDown = input( input.ubound(firstDim), Range::all() );
	Array< Pixel, 1 > lineLeft = input( Range::all(), 0 );
	Array< Pixel, 1 > lineRight = input( Range::all(), input.ubound(secondDim) );

	TinyVector< int, 1 > mUp = minIndex( lineUp );
	TinyVector< int, 1 > mDown = minIndex( lineDown );
	TinyVector< int, 1 > mLeft = minIndex( lineLeft );
	TinyVector< int, 1 > mRight = minIndex( lineRight );

	TinyVector< int, 2 > pUp( 0, mUp[0] );
	TinyVector< int, 2 > pDown( input.ubound(firstDim), mDown[0] );
	TinyVector< int, 2 > pLeft( mLeft[0], 0 );
	TinyVector< int, 2 > pRight( mRight[0], input.ubound(secondDim) );

	Pixel dUp = lineUp( mUp[0] );
	Pixel dDown = lineDown( mDown[0] );
	Pixel dLeft = lineLeft( mLeft[0] );
	Pixel dRight = lineRight( mRight[0] );

	// initialize
	Pixel minDistance = dUp;
	TinyVector< int, 2 > pEdge = pUp;

	if ( ( dLeft < minDistance ) & ( dLeft > 0 ) )
	{
		minDistance = dLeft;
		pEdge = pLeft;
	}
	if ( ( dRight < minDistance ) & ( dRight > 0 ) )
	{
		minDistance = dRight;
		pEdge = pRight;
	}
	// this should not be the case
	if ( ( dDown < minDistance ) & ( dDown > 0 ) )
	{
		cout << "**WARNING** - this should not happen" << endl;
		minDistance = dDown;
		pEdge = pDown;
	}

	return pEdge;
}

template< class Pixel >
void nbfCutGeodesics< Pixel> :: getConstrainedCircularPath( Array< Pixel, 2 > & input,
													        TinyVector< int, 2 > & top,
															TinyVector< int, 2 > & bottom,
															vector< TinyVector< int, 2 > > & path )
{
	//nbfMatlabWriter writer;
	//writer.setFileName("ipath");

	// assume center is baricenter between given points
	//TinyVector< int, 2 > center = ( top + bottom ) / 2.0;
	//this->center = ( top + bottom ) / 2.0;

	// transform weights to radial
	Array< Pixel, 2 > w ( input.shape() );
	w = input - min(input);
	firstIndex i;
	secondIndex j;
	w = w / sqrt( pow2( this->center[0] - i + 0.0 ) + pow2( this->center[1] - j + 0.0 ) );
	w = w + 1e-10;

	//writer.write(w);
	//cout << min(w) << ", " << max(w) << endl;

	nbfFastMarching2D8< Pixel > fastMarching( w, this->center );
	fastMarching.setCutDimension( this->getCutDimension() );
	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;

	cout << top << endl;
	cout << bottom << endl;

	aliveP.push_back( bottom );
	aliveD.push_back(0);
	fastMarching.setAliveSet( aliveP, aliveD );
	Array< Pixel, 2 > distances( w.shape() );
	fastMarching.execute( distances );

	//writer.write(distances);

	// left semi-path
	this->getPath( distances, top, path );	

	//cout << path.size() << endl;

	// right semi-path
	vector< TinyVector< int, 2 > > rightPath;
	this->getPath( distances, TinyVector< int, 2 >( top[firstDim] + 1, top[secondDim] ), rightPath );
	for ( int i = rightPath.size() - 2; i >= 0; i-- )
	{
		path.push_back( rightPath[i] );
	}

	//cout << path.size() << endl;
}


template< class Pixel >
void nbfCutGeodesics< Pixel> :: getFullCircularPath( Array< Pixel, 2 > & input, 
													 vector< TinyVector< Pixel, 2 > > & path,
													 bool polar,
													 Pixel smooth )
{
	nbfMatlabWriter writer;
	writer.setFileName("ipath");

	// transform weights to radial
	Array< Pixel, 2 > w ( input.shape() );
	w = input;

	//writer.write(w);
	//Pixel minW = min(input);
	//w = where( w < numeric_limits<Pixel>::max(), input - minW, w );

	Array< Pixel, 2 > rho ( input.shape() );
	firstIndex indexI;
	secondIndex indexJ;
	rho = sqrt( pow2( this->center[0] - indexI + 0.0 ) + pow2( this->center[1] - indexJ + 0.0 ) );
	rho = where( rho < 3, 0, rho );
	rho = rho / max(rho);

	if ( polar == true ){
		w = where( input < numeric_limits<Pixel>::max(), w / rho + smooth, input );
		w = where( rho > 0, w, numeric_limits<Pixel>::max() );
	}
	else{
		w = where( input > .99, 2, input );
		w = where( input < numeric_limits<Pixel>::max(), w + 1e-1, input );
		w = where( rho > 0, w, numeric_limits<Pixel>::max() );
	}

	//w = where( rho > 3, numeric_limits<Pixel>::max(), w );
	writer.write(w);

	nbfFastMarching2D8< Pixel > fastMarching( w, this->center );
	fastMarching.setCutDimension( this->getCutDimension() );
	//fastMarching.setStopBorder(secondDim);
	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;

	Array< Pixel, 2 > distances( input.shape() );

	// d = 0 on up-left side
	for ( int i = this->center[firstDim]; i < w.rows(); i++ )
	{
		TinyVector< int, 2 > current( i, this->center[secondDim] );
		if ( input( current ) < numeric_limits<Pixel>::max() ){
			aliveP.push_back(current);
			aliveD.push_back(0);
		}
	}			

	fastMarching.setAliveSet( aliveP, aliveD );
	//fastMarching.setStopPoint( TinyVector<int,2>(w.ubound(firstDim)-5,this->center[secondDim]+1) );
	fastMarching.execute( distances );

	writer.write( distances );

	Array< Pixel, 1 > lineA( distances.extent(firstDim) - this->center[firstDim] );
	lineA = distances( Range( this->center[firstDim], toEnd ), center[secondDim] + 1 );
	
	aliveP.clear();
	aliveD.clear();

	// d = 0 on up-right side
	for ( int i = this->center[firstDim]; i < w.rows(); i++ )
	{
		TinyVector< int, 2 > current( i, this->center[secondDim] + 1 );
		if ( input( current ) < numeric_limits<Pixel>::max() ){
			aliveP.push_back(current);
			aliveD.push_back(0);
		}
	}			

	fastMarching.setAliveSet( aliveP, aliveD );
	fastMarching.execute( distances );

	//writer.write(lineA);
	lineA = lineA + distances( Range( this->center[firstDim], toEnd ), center[secondDim] );
	writer.write(lineA);

	TinyVector< int, 1 > mA = minIndex( lineA );

	//cout << mA << endl;

	//cout << this->center << endl;
	//writer.write(lineA);

	// it can happen that a closed path does not exist (e.g. holes in the domain)
	// in this case the returned path is an empty one
	//if ( mA(firstDim) < 0 ){
	//	return;
	//}

	TinyVector< int, 2 > pA( center[firstDim] + mA[0], center[secondDim] );

	//Array< Pixel, 2 > ipath( distances.shape() );

	//TinyVector< int, 2 > st( pA[firstDim], pA[secondDim] + 1 );
	//cout << st << endl;

	//this->getSimplePath( distances, st, path );
	//this->getImplicitPath( path, ipath);
	//writer.write(ipath);

	//// if out of band, look one pixel up or down
	//int i = 1;
	//while ( w( pA ) == numeric_limits<Pixel>::max() ){
	//	if( w( pA[0] - i, pA[1] ) < numeric_limits<Pixel>::max() ){
	//		pA[0] = pA[0] - i;
	//		//cout << pA << ", " << w( pA ) << endl;
	//	}
	//	else{
	//		if( w( pA[0] + i, pA[1] ) < numeric_limits<Pixel>::max() ){
	//			pA[0] = pA[0] + i;
	//		}
	//	}
	//	i++;
	//}

	//if ( w( pA ) == numeric_limits<Pixel>::max() ){
	//	cout << "ERROR" << endl;
	//}

	// set this point as alive
	aliveP.clear();
	aliveD.clear();
	aliveP.push_back( pA );
	aliveD.push_back(0);
	fastMarching.setAliveSet( aliveP, aliveD );

	TinyVector< int, 2 > start( pA[firstDim], pA[secondDim] + 1 );

	//cout << pA << endl;
	//cout << start << endl;

	fastMarching.setStopPoint( start );
	fastMarching.unSetStopBorder();
	fastMarching.execute( distances );

	writer.write(distances);

	this->getSimplePath( distances, start, path );
	//this->getImplicitPath( path, ipath);
	//writer.write(ipath);
}

template< class Pixel >
void nbfCutGeodesics< Pixel> :: getSemiCircularPath( Array< Pixel, 2 > & input, 
													 vector< TinyVector< int, 2 > > & path )
{
	// transform weights to radial
	Array< Pixel, 2 > w ( input.shape() );
	w = input - min(input);
	firstIndex i;
	secondIndex j;
	w = w / sqrt( pow2( this->center[0] - i + 0.0 ) + pow2( this->center[1] - j + 0.0 ) );
	w = w + 1e-10;

	nbfFastMarching2D8< Pixel > fastMarching( w );
	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;

	Array< Pixel, 2 > distancesFor( input.shape() );
	Array< Pixel, 2 > distancesBack( input.shape() );

	// forward

	// d = 0 on up-left side
	for ( int j = 0; j < center[secondDim]; j++ )
	{
		aliveP.push_back( TinyVector< int, 2 >(0,j) );
		aliveD.push_back(0);
	}			

	fastMarching.setAliveSet( aliveP, aliveD );
	fastMarching.execute( distancesFor );

	// look for first arriving point on the up-right side
	Array< Pixel, 1 > lineAfor = distancesFor( 0, Range( center[secondDim] + 1, toEnd ) );
	TinyVector< int, 1 > mAfor = minIndex( lineAfor );
	TinyVector< int, 2 > pAfor( 0, mAfor[0] + center[secondDim] + 1 );

	// d = 0 on the first arriving point on the up-right side
	aliveP.clear();
	aliveD.clear();
	aliveP.push_back( pAfor );
	aliveD.push_back(0);

	fastMarching.setAliveSet( aliveP, aliveD );
	fastMarching.execute( distancesFor );

	// look for first arriving point on the up-left side
	Array< Pixel, 1 > lineAback = distancesFor( 0, Range( 0, center[secondDim] - 1 ) );
	TinyVector< int, 1 > mAback = minIndex( lineAback );
	TinyVector< int, 2 > pAback( 0, mAback[0] );

	// try the geodesic in the other direction

	aliveP.clear();
	aliveD.clear();
	aliveP.push_back( pAback );
	aliveD.push_back(0);
	fastMarching.setAliveSet( aliveP, aliveD );
	fastMarching.setStopPoint( pAfor );
	fastMarching.execute( distancesBack );

	// use the one that gives minimum distance
	nbfGeodesicPath< Pixel > geodesic2D( distancesFor );
    if ( distancesFor( pAback ) < distancesBack( pAfor ) ){
		geodesic2D.getPath( distancesFor, pAback, path );
	}
	else{
		geodesic2D.getPath( distancesBack, pAback, path );
	}
	
#if 0 // DO NOT CONSIDER HALFWAY OBJECTS - ONLY OBJECTS COMPLETELY INSIDE VOLUME

	// B case (halfway out - left)

	// B - forward

	// d = 0 on left side
	aliveP.clear();
	aliveD.clear();

	for ( int i = 0; i < input.rows(); i++ )
	{
		aliveP.push_back( TinyVector< int, 2 >(i,0) );
		aliveD.push_back(0);
	}			

	fastMarching.setAliveSet( aliveP, aliveD );
	fastMarching.execute( distances );

	Array< Pixel, 1 > lineBfor = distances( 0, Range( center[secondDim] + 1, toEnd ) );
	TinyVector< int, 1 > mBfor = minIndex( lineBfor );
	TinyVector< int, 2 > pBfor( 0, mBfor[0] + center[secondDim] + 1 );
	Pixel dBfor = lineBfor( mBfor[0] );

    bool stateBfor = this->getBackwardPath( pBfor, pathFor );
	Pixel pathBfor = pathFor.size();

	//ipath = 1;
	//for ( int i = 0; i < pathFor.size(); i++ )
	//{
	//	ipath(pathFor[i]) = 0;
	//}
	//writer.write(ipath);

	// B - backward

	// d = 0 on up-right side
	aliveP.clear();
	aliveD.clear();
	for ( int j = center[secondDim] + 1; j < input.cols(); j++ )
	{
	
		aliveP.push_back( TinyVector< int, 2 >( 0, j ) );
		aliveD.push_back(0);
	}			

	fastMarching.setAliveSet( aliveP, aliveD );
	fastMarching.execute( distances );

	// arrival
	Array< Pixel, 1 > lineBback = distances( Range::all(), 0 );
	TinyVector< int, 1 > mBback = minIndex( lineBback );
	TinyVector< int, 2 > pBback( mBback[0], 0 );
	Pixel dBback = lineBback( mBback[0] );

    bool stateBback = this->getForwardPath( pBback, pathBack );
	Pixel pathBback = pathBack.size();

	//ipath = 1;
	//for ( int i = 0; i < pathBack.size(); i++ )
	//{
	//	ipath(pathBack[i]) = 0;
	//}
	//writer.write(ipath);

	cout << "f = " << pathBfor << ", b = " << pathBback << endl;

	// keep the shortest path

	if ( stateBfor & ( pathBfor < minimalPathMeasure ) )
	{
		minimalPathMeasure = pathBfor;
		path.clear();
		path.resize( pathFor.size() );
		copy(pathFor.begin(),pathFor.end(),path.begin());
	}

	if ( stateBback & ( pathBback < minimalPathMeasure ) )
	{
		minimalPathMeasure = pathBback;
		path.clear();
		path.resize( pathBack.size() );
		copy(pathBack.begin(),pathBack.end(),path.begin());
	}

    // C case (halfway out - right)

	// C - forward

	// d = 0 on up-left side
	aliveP.clear();
	aliveD.clear();

	for ( int j = 0; j < center[secondDim]; j++ )
	{
		aliveP.push_back( TinyVector< int, 2 >(0,j) );
		aliveD.push_back(0);
	}			

	fastMarching.setAliveSet( aliveP, aliveD );
	fastMarching.execute( distances );

	Array< Pixel, 1 > lineCfor = distances( Range::all(), distances.ubound(secondDim) );
	TinyVector< int, 1 > mCfor = minIndex( lineCfor );
	TinyVector< int, 2 > pCfor( mCfor[0], distances.ubound(secondDim) );
	Pixel dCfor = lineCfor( mCfor[0] );

    bool stateCfor = this->getBackwardPath( pCfor, pathFor );
	Pixel pathCfor = pathFor.size();

	//ipath = 1;
	//for ( int i = 0; i < pathFor.size(); i++ )
	//{
	//	ipath(pathFor[i]) = 0;
	//}
	//writer.write(ipath);

	// C - backward

	// d = 0 on right side
	aliveP.clear();
	aliveD.clear();
	for ( int i = 0; i < input.rows(); i++ )
	{
	
		aliveP.push_back( TinyVector< int, 2 >( i, input.ubound(secondDim) ) );
		aliveD.push_back(0);
	}			

	fastMarching.setAliveSet( aliveP, aliveD );
	fastMarching.execute( distances );

	// arrival
	Array< Pixel, 1 > lineCback = distances( 0, Range( 0, center[secondDim] - 1 ) );
	TinyVector< int, 1 > mCback = minIndex( lineCback );
	TinyVector< int, 2 > pCback( 0, mCback[0] );
	Pixel dCback = lineCback( mCback[0] );

    bool stateCback = this->getForwardPath( pCback, pathBack );
	Pixel pathCback = pathBack.size();

	//ipath = 1;
	//for ( int i = 0; i < pathBack.size(); i++ )
	//{
	//	ipath(pathBack[i]) = 0;
	//}
	//writer.write(ipath);

	cout << "f = " << pathCfor << ", b = " << pathCback << endl;

	// keep the shortest path

	if ( stateCfor & ( pathCfor < minimalPathMeasure ) )
	{
		minimalPathMeasure = pathCfor;
		path.clear();
		path.resize( pathFor.size() );
		copy(pathFor.begin(),pathFor.end(),path.begin());
	}

	if ( stateCback & ( pathCback < minimalPathMeasure ) )
	{
		path.clear();
		path.resize( pathBack.size() );
		copy(pathBack.begin(),pathBack.end(),path.begin());
	}
#endif
}

template< class Pixel >
void nbfCutGeodesics< Pixel> :: getConstrainedSemiCircularPath( Array< Pixel, 2 > & input,
													            TinyVector< int, 2 > & center,
															    TinyVector< int, 2 > & start,
															    TinyVector< int, 2 > & end,
															    vector< TinyVector< int, 2 > > & path )
{
	//nbfMatlabWriter writer;
	//writer.setFileName("ipath");

	//cout << center << endl;
	//cout << start << endl;
	//cout << end << endl;

	// transform weights to radial
	Array< Pixel, 2 > w ( input.shape() );
	w = input - min(input);
	firstIndex i;
	secondIndex j;
	w = w / sqrt( pow2( center[0] - i + 0.0 ) + pow2( center[1] - j + 0.0 ) );
	w = w + 1e-10;

	nbfFastMarching2D8< Pixel > fastMarching(w);
	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;

	aliveP.push_back(start);
	aliveD.push_back(0);
	fastMarching.setAliveSet( aliveP, aliveD );
	fastMarching.setStopPoint(end);
	Array< Pixel, 2 > distancesFor( w.shape() );
	fastMarching.execute( distancesFor );

	//writer.write( distancesFor );

	aliveP.clear();
	aliveD.clear();
	aliveP.push_back(end);
	aliveD.push_back(0);
	fastMarching.setAliveSet( aliveP, aliveD );
	fastMarching.setStopPoint(start);
	Array< Pixel, 2 > distancesBack( w.shape() );
	fastMarching.execute( distancesBack );

	//writer.write( distancesBack );

	nbfCutGeodesics< Pixel > cutGeodesic( distancesBack );
	if ( distancesFor(end) < distancesBack(start) ){
        cutGeodesic.getPath( distancesFor, end, path );	
	}
	else{
        cutGeodesic.getPath( distancesBack, start, path );	
	}
}

template< class Pixel >
void nbfCutGeodesics< Pixel> :: curvature( Array< Pixel, 2 > & input, int iterations )
{
	Pixel timeStep = .25;
	BordStrategyMirrorDouble< Pixel, 2 > BSforInput( input, 1 );
	Array< Pixel, 2 > R( input.shape() );
	for ( int i = 0; i < iterations; i++ ){
		R = curvature2D(input);
		input = input + timeStep * R;
		BSforInput.refresh();
	}
}

#endif /* FILE_nbfCutGeodesics */