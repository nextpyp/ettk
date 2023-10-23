#ifndef FILE_nbfGeodesicPath
#define FILE_nbfGeodesicPath

// Class nbfGeodesicPath.
//
// Find geodesic paths in 2D images.
// 

#include <fm/nbfFastMarchingFool.h>
#include <fm/nbfFastMarching2D8.h>
#include <fm/nbfFastFastMarching2D.h>
#include <nbfLinearInterpolator.h>

#include <vector>
#include <algorithm>

template< class Pixel >
class nbfGeodesicPath : public nbfArray< Pixel, 2 >
{
public:

	// takes $g$ array and center point as input
	nbfGeodesicPath( Array< Pixel, 2 > &, TinyVector< int, 2 > & );
	nbfGeodesicPath( Array< Pixel, 2 > & );

	~nbfGeodesicPath(){};

	// get geodesic curve given start point
	// the end point is assumed to have zero distance
	// return true if succesful, false otherwise

	// (weights, start, path)
	bool getPath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, vector< TinyVector< int, 2 > > & );
	
	// (weights, start, implicit path)
	bool getPath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, Array< Pixel, 2 > & );

	// (weights, start, path)
	bool getGoodPath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, vector< TinyVector< Pixel, 2 > > & );
	bool getSimplePath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, vector< TinyVector< Pixel, 2 > > & );

	void getPointGradient( Array< Pixel, 2 > &, TinyVector< int, 2 > &, TinyVector< Pixel, 2 > & );
	void getGradient( Array< Pixel, 2 > &, Array< Pixel, 2 > &, Array< Pixel, 2 > & );

	// restrict path to rounded (clock or counter-clock wise)
	bool getForwardPath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, vector< TinyVector< int, 2 > > & );	
	bool getBackwardPath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, vector< TinyVector< int, 2 > > & );

	// get circular (same starting and end height) geodesic
	Pixel getCircularPath( Array< Pixel, 2 > &, vector< TinyVector< Pixel, 2 > > & );

	// gives the first arriving point on the right side of the image.
	// we set to alive all point in the left-hand side, and let fast marching find the
	// first arriving point on the right side
	void getFirstPointOnRightSide( Array< Pixel, 2 > &, TinyVector< Pixel,2 > &, TinyVector< Pixel, 2 > & );

	// get path from one side to the other (not neccesarily circular)
	bool getPathBetweenSides( Array< Pixel, 2 > &, vector< TinyVector< Pixel, 2 > > & );

	void getForwardLinearPath( Array< Pixel, 2 > &, TinyVector< int, 2 > &, vector< TinyVector< int, 2 > > & );

	// get path as implicit function (only closed paths)
	void getImplicitPath( vector< TinyVector< int, 2 > > &, Array< Pixel, 2 > & );
	void getImplicitPath( vector< TinyVector< Pixel, 2 > > &, Array< Pixel, 2 > &, Pixel = numeric_limits<Pixel>::max() );

	// side-to-side open paths
	void getImplicitOpenPath( vector< TinyVector< int, 2 > > &, Array< Pixel, 2 > & );

protected:

	// check if path valid (halt if too long for example)
	bool checkPath( vector< TinyVector< int, 2 > > & );

	// update theta zones
	void updateTheta();

	// weight array
	Array< Pixel, 2 > weight;

	// theta values
	Array< Pixel, 2 > theta;
};


template< class Pixel >
nbfGeodesicPath< Pixel > :: nbfGeodesicPath( Array< Pixel, 2 > & weight )
: nbfArray< Pixel, 2 >( weight )
{
}

template< class Pixel >
nbfGeodesicPath< Pixel > :: nbfGeodesicPath( Array< Pixel, 2 > & weight,
												   TinyVector< int, 2 > & center )
												  : nbfArray< Pixel, 2 >( weight, center )
{
	this->weight.reference( weight );
	this->center = center;
	this->updateTheta();
}


template< class Pixel >
void nbfGeodesicPath< Pixel > :: updateTheta()
{
	this->getTheta( this->theta );

	// convert to 0 - 2PI => convert to 0 - 8 => convert to 0.5 - 8.5, => quantize
	this->theta = floor( where( this->theta < 0, this->theta + 2.0 * vtkMath::Pi(), this->theta ) / vtkMath::Pi() / 2.0 * 8.0 + .5 );
	
	// make periodic
	this->theta = where( this->theta == 8.0, 0.0, this->theta );
}


template< class Pixel >
bool nbfGeodesicPath< Pixel> :: getPath( Array< Pixel, 2 > & weights, 
									     TinyVector< int, 2 > & start,
										 vector< TinyVector< int, 2 > > & path )
{	
	// reset if not empty
	path.clear();

	TinyVector< int, 2 > next;
	
	vector< TinyVector< int, 2 > > neighbors;

	// back propagation - minimum distance neighbor
	next = start;
	path.push_back( next );
	
	Pixel minDistance;

	// stop when reach null distance
	while ( weights( next ) != 0 )
	{
		bool hasMultipleMinima = false;

		this->getFullNeighbors( next, neighbors );

		TinyVector< int, 2 > current = next;
		
		if ( neighbors.size() < 1 ){
			break;
		}
		else
		{
			minDistance = weights(neighbors[0]);
			next = neighbors[0];
			for ( int i = 1; i < neighbors.size(); i++ ){
				if ( weights( neighbors[i] ) < minDistance ){
					minDistance = weights( neighbors[i] );
					next = neighbors[i];
					hasMultipleMinima = false;
				}
				else if ( weights( neighbors[i] ) == minDistance ){ 
					hasMultipleMinima = true;
					//cout << "Multiple minima = " << minDistance << endl;
				}
			}

			if ( hasMultipleMinima == true ){
				Pixel minD = numeric_limits<Pixel>::max();
				Pixel minIndex = 0;
				for ( int i = 0; i < neighbors.size(); i++ ){
					Pixel d = pow2( neighbors[i](firstDim) - current(firstDim) ) + pow2( neighbors[i](secondDim) - current(secondDim) );
					d = ( weights( neighbors[i] ) - weights( current ) ) / sqrt( d );
					if ( d < minD ){
						minD = d;
						minIndex = i;
					}
				}
				next = neighbors[minIndex];
			}

			if ( !this->checkPath( path ) ){
				return false;
			}			
			path.push_back( next );
		}
	}
	return true;
}


template< class Pixel >
bool nbfGeodesicPath< Pixel> :: getPath( Array< Pixel, 2 > & weights, 
									     TinyVector< int, 2 > & start,
										 Array< Pixel, 2 > & ipath )
{	
	vector< TinyVector< int, 2 > > path;
	bool state = this->getPath( weights, start, path );
	this->getImplicitPath( path, ipath );
	return state;
}


// Here we use a true back propagation, computing the gradient with bilinear interpolation
// As distances are computed with 8-neighbors we can get stucked (numerical errors) so
// a procedure to handle those cases is provided, if we ever go back to a position already
// in the path, we replace it with the 8-neighbor with the smallest distance,
// and this is equivalent to switching to the first order back propagation.
template< class Pixel >
bool nbfGeodesicPath< Pixel> :: getGoodPath( Array< Pixel, 2 > & weights, 
									         TinyVector< int, 2 > & start,
										     vector< TinyVector< Pixel, 2 > > & path )
{	
	// reset path if not empty
	path.clear();

	// store next point in geodesic
	TinyVector< Pixel, 2 > next;
	
	// store current point neighbors
	vector< TinyVector< Pixel, 2 > > neighbors;

	// push first point into path
	next = start;
	path.push_back( next );
	
	Pixel minDistance;

	// interpolator to get bilinear gradient estimation
	nbfLinearInterpolator< Pixel, 2 > interp( weights );

	// this is the time step, we use a fixed value. 
	// Gradients (this is the incremental update) are normalized to [0,1]
	// So 'delta' is telling us that we move at most delta in the direction
	// of the gradient, the grid spacing being 1.
	Pixel delta = .5;

	// flag for signaling when the backprop gets stuck.
	bool stucked = false;

	// loop until we reach zero distance
	while ( weights( (int)next(firstDim), (int)next(secondDim) ) > 0 )
	{
		cout << next << endl;
		// made last position coordinates handy
		Pixel x = next[firstDim];
		Pixel y = next[secondDim];

		// get square box where the point is contained
		int lx = floor( x ); int ux = ceil( x );
		int ly = floor( y ); int uy = ceil( y );
		
		// compute gradient value from bilinear interpolation
		TinyVector< Pixel, 2 > gradient = interp.interpolateSingleGradient( next );

		// Switch to upwing gradient on grid lines/points.
		// when we are close enough to a grid point (1e-3), we switch to this approximation.

		// store position where we got stuck
		TinyVector< int, 2 > stuckPosition;

		// x dimension
		stuckPosition = next;

		if ( ( next[firstDim] - lx ) < 1e3 ){
			if ( ( lx > 0 ) & ( lx < weights.ubound(firstDim) ) ){ // if within array limits
				if ( ( weights(lx+1,ly) < weights(lx,ly) ) & ( weights(lx-1,ly) < weights(lx,ly) ) ){
					if ( weights(lx-1,ly) > weights(lx+1,ly) ){
						gradient(firstDim) = weights(lx+1,ly) - weights(lx,ly);
					}
					else{
						gradient(firstDim) = weights(lx,ly) - weights(lx-1,ly);
					}
				}
				else{
					if ( weights(lx+1,ly) < weights(lx,ly) ){
						gradient(firstDim) = weights(lx+1,ly) - weights(lx,ly);
					}
					else{
						if ( weights(lx-1,ly) < weights(lx,ly) ){
							gradient(firstDim) = weights(lx,ly) - weights(lx-1,ly);
						}
						else{
							stucked = true;
						}
					}
				}
			}
			else{
                if ( lx == 0 ){
					if ( weights(lx+1,ly) < weights(lx,ly) ){
						gradient(firstDim) = weights(lx+1,ly) - weights(lx,ly);
					}
					else{
						stucked = true;
					}
				}
				if ( lx == weights.ubound(firstDim) ){
					if ( weights(lx-1,ly) < weights(lx,ly) ){
						gradient(firstDim) = weights(lx,ly) - weights(lx-1,ly);
					}
					else{
						stucked = true;
					}
				}
			}
		}
		// y dimension
		if ( ( next[secondDim] - ly ) < 1e3 ){
			if ( ( ly > 0 ) & ( ly < weights.ubound(secondDim) ) ){
				if ( ( weights(lx,ly-1) < weights(lx,ly) ) & ( weights(lx,ly+1) < weights(lx,ly) ) ){
					if ( weights(lx,ly-1) > weights(lx,ly+1) ){
						gradient(secondDim) = weights(lx,ly+1) - weights(lx,ly);
					}
					else{
						gradient(secondDim) = weights(lx,ly) - weights(lx,ly-1);
					}
				}
				else{
					if ( weights(lx,ly-1) < weights(lx,ly) ){
						gradient(secondDim) = weights(lx,ly) - weights(lx,ly-1);
					}
					else{
						if (  weights(lx,ly+1) < weights(lx,ly) ){
							gradient(secondDim) = weights(lx,ly+1) - weights(lx,ly);
						}
						else{
							stucked = true;
						}
					}
				}
			}
			else{
				if ( ly == 0 ){
					if ( weights(lx,ly+1) < weights(lx,ly) ){
						gradient(secondDim) = weights(lx,ly+1) - weights(lx,ly);
					}
					else{
						stucked = true;
					}
				}
                if ( ly == weights.ubound(secondDim) ){
					if ( weights(lx,ly-1) < weights(lx,ly) ){
						gradient(secondDim) = weights(lx,ly) - weights(lx,ly-1);
					}
					else{
						stucked = true;
					}
				}
			}
		}

		// Check if the bilinear gradient approximation is congruent with the stepest descent direction
		// given by the neighboring pixel with the minimum distance value.
		// We check that the gradient direction and the stepest direction are in the same direction,
		// if that is not the case, we substitute the bilinear approximation by the
		// stepest descent direction. This direction is given by the vector from the current point
		// to the neighboring pixel with the minimum distance and guarantees that we march towards the 
		// global minima.
		if ( ( weights(lx,ly) < weights(lx,uy) ) & 
			 ( weights(lx,ly) < weights(ux,ly) ) & 
			 ( weights(lx,ly) < weights(ux,uy) ) ){
				 TinyVector< Pixel, 2 > minDir( lx - x, ly - y );
				 if ( gradient[firstDim] * minDir[firstDim] + gradient[secondDim] * minDir[secondDim] >= 0 ){
					gradient = - minDir;
				 }
			 }
		if ( ( weights(lx,uy) < weights(lx,ly) ) & 
			 ( weights(lx,uy) < weights(ux,ly) ) & 
			 ( weights(lx,uy) < weights(ux,uy) ) ){
				 TinyVector< Pixel, 2 > minDir( lx - x, uy - y );
				 if ( gradient[firstDim] * minDir[firstDim] + gradient[secondDim] * minDir[secondDim] >= 0 ){
					gradient = - minDir;
				 }
			 }
		if ( ( weights(ux,ly) < weights(ux,uy) ) & 
			 ( weights(ux,ly) < weights(lx,ly) ) & 
			 ( weights(ux,ly) < weights(lx,uy) ) ){
				 TinyVector< Pixel, 2 > minDir( ux - x, ly - y );
				 if ( gradient[firstDim] * minDir[firstDim] + gradient[secondDim] * minDir[secondDim] >= 0 ){
					gradient = - minDir;
				 }
			 }
		if ( ( weights(ux,uy) < weights(ux,ly) ) & 
			 ( weights(ux,uy) < weights(lx,ly) ) & 
			 ( weights(ux,uy) < weights(lx,uy) ) ){
				 TinyVector< Pixel, 2 > minDir( ux - x, uy - y );
				 if ( gradient[firstDim] * minDir[firstDim] + gradient[secondDim] * minDir[secondDim] >= 0 ){
					gradient = - minDir;
				 }
			 }

		// normalize gradient approximation to [0,1]
		Pixel norm = sqrt( pow2( gradient(firstDim) ) + pow2( gradient(secondDim) ) );
		gradient(firstDim) /= norm;
		gradient(secondDim) /= norm;

		//cout << gradient << endl;

		// The path can oscillate because of the fixed delta. When we are close to the geodesic the value
		// of delta can be too big, so we go to the oposite side of the geodesic.
		// We check that both the x and y axis are not crossed simultaneously. When this happens
		// it means that we jumped one voxel, potentially missing the stepest direction.
		// In this case we set the next point to be the 'jumped' voxel.
		if ( ( lx != floor( lx - delta * gradient(firstDim) ) ) &   // detect x crossing
			 ( ly != floor( ly - delta * gradient(secondDim) ) ) ){ // detect y crossing
				 stucked = true;
				 // find stuck position
				 if ( lx - delta * gradient(firstDim) > lx ){
						stuckPosition[firstDim] = lx + 1;
					}
				 else{
						stuckPosition[firstDim] = lx;
				 }
				 if ( ly - delta * gradient(secondDim) > ly ){
					stuckPosition[secondDim] = ly + 1;
					}
				 else{
					stuckPosition[secondDim] = ly;
				 }
			 }
		else{
			// advance 'delta' in the direction contrary to the gradient
			next(firstDim) = next(firstDim) - delta * gradient(firstDim);
			next(secondDim) = next(secondDim) - delta * gradient(secondDim);
		}

		// Check if we got stuck
		if ( stucked == true ){
			// reset flag
			stucked = false;

			// look for stepest descent direction (neighbor with smallest distance value)
			vector< TinyVector< int, 2 > > neighbors;
			path.push_back( stuckPosition );
			this->getFullNeighbors( stuckPosition, neighbors );
			Pixel minDistance = weights(neighbors[0]);
			next = neighbors[0];
			for ( int i = 1; i < neighbors.size(); i++ ){
				if ( weights( neighbors[i] ) < minDistance ){
					minDistance = weights( neighbors[i] );
					next = neighbors[i];
				}
			}
		}

		// it should never get stucked, but provide a stoping criteria if failure
		if ( path.size() > 5000 ){
			// DEBUG
			nbfMatlabWriter writer;
			writer.setFileName("ipath");
			writer.write(weights);
			for ( int i = 0; i < path.size(); i++ ){
				cout << path[i] << endl;
			}
			cout << weights( Range( next[0] - 3, next[0] + 3 ), Range( next[1] - 3, next[1] + 3 ) ) << endl;
			return false;
		}			

		path.push_back( next );
	}
	next(firstDim) = floor(next(firstDim));
	next(secondDim) = floor(next(secondDim));
	path.push_back( next );
	return true;
}

template< class Pixel >
void nbfGeodesicPath< Pixel> :: getPointGradient( Array< Pixel, 2 > & weights, 
												  TinyVector< int, 2 > & position,
												  TinyVector< Pixel, 2 > & gradient )
{	
	//if ( weights(position) == numeric_limits<Pixel>::max() ){
	//	int h = 1;
	//}

	vector< TinyVector< int, 2 > > neighbors;
	this->getFullNeighbors( position, neighbors );

	// look for neighbor position with greatest gradient
	TinyVector< int, 2 > next = neighbors[0];
	Pixel minGradient = ( weights(next) - weights(position) ) / sqrt( pow2(next[0]-position[0])+pow2(next[1]-position[1]) + 0.0 );
	for ( int i = 1; i < neighbors.size(); i++ ){
		Pixel gradient = ( weights(neighbors[i]) - weights(position) ) / sqrt(pow2(neighbors[i][0]-position[0])+pow2(neighbors[i][1]-position[1]) + 0.0);
		if ( gradient < minGradient ){
			minGradient = gradient;
			next = neighbors[i];
		}
	}

	TinyVector< int, 2 > pA; // grid neighbor point
	TinyVector< int, 2 > pB; // diagonal point

	TinyVector< int, 2 > p, q;

	// if diagonal point
	if ( ( position[firstDim] != next[firstDim] ) && ( position[secondDim] != next[secondDim] ) ){
		pB = next;
		p = position + TinyVector< int, 2 >( pB[firstDim] - position[firstDim], 0 );
		q = position + TinyVector< int, 2 >( 0, pB[secondDim] - position[secondDim] );
		if ( weights(p) < weights(q) ){
			pA = p;
		}
		else{
			pA = q;
		}
	}
	// look for smallest corner point
	else{
		pA = next;
		if ( position[firstDim] == pA[firstDim] ){
			p = pA + TinyVector< int, 2 >(1,0);
			q = pA + TinyVector< int, 2 >(-1,0);
		}
		else{
			p = pA + TinyVector< int, 2 >(0,1);
			q = pA + TinyVector< int, 2 >(0,-1);
		}
		if ( weights.isInRange(p) && ( this->isAtCut(position,p) == false ) ){
			if ( weights.isInRange(q) && ( this->isAtCut(position,q) == false ) ){
				if ( weights(q) < weights(p) ){
					pB = q;
				}
				else{
					pB = p;
				}
			}
			else{
				pB = p;
			}
		}
		else{
			pB = q;
		}
	}

	//if ( ( position[firstDim] == 25 ) & ( position[secondDim] == 250 ) ){
	//	cout << position << " = " << weights(position) << endl;
	//	cout << pA << " = " << weights(pA) << endl;
	//	cout << pB << " = " << weights(pB) << endl;
	//}

	// compute gradient from computed triangle

	if ( weights(pA) > weights(position) ){
		gradient = position - pB;
	}
	else{
		// compute gradient estimation
		if ( weights(pA) <= weights(pB) ){
			gradient = position - pA;
		}
		else{
			if ( weights(pB) <= ( 2 * weights(pA) - weights(position) ) ){
				gradient = position - pB;
			}
			else{
				if ( position[secondDim] == pA[secondDim] ){
					if ( position[firstDim] > pA[firstDim] ){
						gradient[firstDim] = weights(position) - weights(pA);
					}
					else{
						gradient[firstDim] = weights(pA) - weights(position);
					}
					if ( pA[secondDim] > pB[secondDim] ){
						gradient[secondDim] = weights(pA) - weights(pB);
					}
					else{
						gradient[secondDim] = weights(pB) - weights(pA);
					}
				}
				else{
					if ( position[secondDim] > pA[secondDim] ){
						gradient[secondDim] = weights(position) - weights(pA);
					}
					else{
						gradient[secondDim] = weights(pA) - weights(position);
					}
					if ( pA[firstDim] > pB[firstDim] ){
						gradient[firstDim] = weights(pA) - weights(pB);
					}
					else{
						gradient[firstDim] = weights(pB) - weights(pA);
					}
				}
			}
		}
	}

	// normalize gradient vector
	Pixel norm = sqrt( pow2( gradient[firstDim] ) + pow2( gradient[secondDim] ) );
	gradient = gradient / norm;

	//if ( position[secondDim] == 80 ){
	//	cout << neighbors.size() << endl;
	//	cout << gradient << endl;
	//	cout << weights(pA) << ", " << weights(position) << ", " << weights(pB) << endl;
	//}

	// handle the origin
	if ( ( minGradient == 0 ) || ( ( weights(pA) > weights(position) ) & ( weights(pB) > weights(position) ) ) ){
//		cout << weights( Range(position[firstDim]-1,position[firstDim]+1), Range(position[secondDim]-1,position[secondDim]+1) ) << endl;
		gradient = 0;
	}
}

template< class Pixel >
void nbfGeodesicPath< Pixel> :: getGradient( Array< Pixel, 2 > & weights, 
									         Array< Pixel, 2 > & dx,
											 Array< Pixel, 2 > & dy )
{	
	vector< TinyVector< int, 2 > > neighbors;
	
	TinyVector< Pixel, 2 > gradient;
		
	Array< Pixel, 2 > :: iterator iter = weights.begin();
	
	while( iter != weights.end() ){
		
		TinyVector< int, 2 > position = iter.position();

		// compute gradient from computed triangle
		this->getPointGradient( weights, position, gradient );

		// assign values
		dx( position ) = gradient[firstDim];
		dy( position ) = gradient[secondDim];

		++iter;
	}
}

// This back propagation (state-of-the-art) first computes a gradient image
// and then interpolates to get the gradient at intermediate steps.
// The gradient is computed assuming that distances were computed with 8-neighbor fast
// marching, so the gradient direction at each pixel is that from where the distance
// was updated.
// We now bilinearly interpolate the x and y components of the gradient image
// to get the advancing direction at arbitrary point in the image domain.

template< class Pixel >
bool nbfGeodesicPath< Pixel> :: getSimplePath( Array< Pixel, 2 > & weights, 
									           TinyVector< int, 2 > & start,
										       vector< TinyVector< Pixel, 2 > > & path )
{	
	// reset path if not empty
	path.clear();

	// store next point in geodesic
	TinyVector< Pixel, 2 > next;
	
	// push first point into path
	next = start;
	path.push_back( next );
	
	// posta
	Array< Pixel, 2 > dx( weights.shape() );
	Array< Pixel, 2 > dy( weights.shape() );

	dx = numeric_limits<Pixel>::max();
	dy = numeric_limits<Pixel>::max();
	
	//this->getGradient( weights, dx, dy );

	//nbfMatlabWriter writer;
	//writer.setFileName("ipath");
	//writer.write(weights);
	//writer.write(dy);

	nbfLinearInterpolator< Pixel, 2 > interpX( dx );
	nbfLinearInterpolator< Pixel, 2 > interpY( dy );

	TinyVector< Pixel, 2 > gridGradient;

	vector< TinyVector< int, 2 > > neighbors;

	while ( weights( next ) > 0 ){

		// made last position coordinates handy
		Pixel x = next[firstDim];
		Pixel y = next[secondDim];

		// get square box where the point is contained
		int lx = floor( x ); int ux = ceil( x );
		int ly = floor( y ); int uy = ceil( y );

		// Before interpolating the gradient, make sure all adjacent points 
		// have a gradient value (costly operation)
		neighbors.clear();
		if ( ( lx == ux ) & ( ly == uy ) ){
			// do all 9 points
			neighbors.push_back( TinyVector< int, 2 >( lx - 1, ly - 1 ) );
			neighbors.push_back( TinyVector< int, 2 >( lx - 1, ly) );
			neighbors.push_back( TinyVector< int, 2 >( lx - 1, ly + 1 ) );
			neighbors.push_back( TinyVector< int, 2 >( lx, ly - 1 ) );
			neighbors.push_back( TinyVector< int, 2 >( lx, ly) );
			neighbors.push_back( TinyVector< int, 2 >( lx, ly + 1 ) );
			neighbors.push_back( TinyVector< int, 2 >( lx + 1, ly - 1 ) );
			neighbors.push_back( TinyVector< int, 2 >( lx + 1, ly) );
			neighbors.push_back( TinyVector< int, 2 >( lx + 1, ly + 1 ) );
			this->restrictNeighbors( TinyVector<int,2>(lx,ly), neighbors );
			//cout << neighbors.size() << endl;
		}
		else{
			if ( lx == ux ){
				// 6 point neighborhood
				neighbors.push_back( TinyVector< int, 2 >( lx - 1, ly ) );
				neighbors.push_back( TinyVector< int, 2 >( lx - 1, uy ) );
				neighbors.push_back( TinyVector< int, 2 >( lx, ly ) );
				neighbors.push_back( TinyVector< int, 2 >( lx, uy ) );
				neighbors.push_back( TinyVector< int, 2 >( lx + 1, ly ) );
				neighbors.push_back( TinyVector< int, 2 >( lx + 1, uy ) );
			}
			else{
				if ( ly == uy ){
					// 6 point neighborhood
					neighbors.push_back( TinyVector< int, 2 >( lx, ly - 1 ) );
					neighbors.push_back( TinyVector< int, 2 >( ux, ly - 1 ) );
					neighbors.push_back( TinyVector< int, 2 >( lx, ly ) );
					neighbors.push_back( TinyVector< int, 2 >( ux, ly ) );
					neighbors.push_back( TinyVector< int, 2 >( lx, ly + 1 ) );
					neighbors.push_back( TinyVector< int, 2 >( ux, ly + 1 ) );
				}
				else{
					neighbors.push_back( TinyVector< int, 2 >( lx, ly ) );
					neighbors.push_back( TinyVector< int, 2 >( lx, uy ) );
					neighbors.push_back( TinyVector< int, 2 >( ux, ly ) );
					neighbors.push_back( TinyVector< int, 2 >( ux, uy ) );
				}
			}
		}

		vector< TinyVector< int, 2 > > :: iterator iter = neighbors.begin();
		while ( iter != neighbors.end() ){
			if ( dx.isInRange(*iter) == true ){
				if ( dx( *iter ) == numeric_limits<Pixel>::max() ){
					if ( weights(*iter) < numeric_limits<Pixel>::max() ){
						this->getPointGradient( weights, *iter, gridGradient );
					}
					else{
						gridGradient = 0;
					}
					//cout << *iter << ", g = " << gridGradient << endl;
					dx( *iter ) = gridGradient[firstDim];
					dy( *iter ) = gridGradient[secondDim];
				}
			}
			++iter;
		}

		TinyVector< Pixel, 2 > gradient( interpX.interpolateSingle( next ), 
			                             interpY.interpolateSingle( next ) );

		// cout << next << ", g = " << gradient << endl;
		
		// re-normalize gradient vector
		Pixel norm = sqrt( pow2( gradient[firstDim] ) + pow2( gradient[secondDim] ) );

		if ( norm > 0 ){
			gradient = gradient / norm;
		}
		else{
			cout << next << endl;
			cout << gradient << endl;
			for ( int i = 0; i < neighbors.size(); i++ ){
				cout << neighbors[i] << ", [" << dx(neighbors[i]) << "," << dy(neighbors[i]) << "]\n";
			}
			return false;
		}

		// es muy poco sensible al tamaño del paso, achicarlo solo genera discrepancias minimas
		Pixel lambda = .25;
		
		next[firstDim] = next[firstDim] - lambda * gradient[firstDim];
		next[secondDim] = next[secondDim] - lambda * gradient[secondDim];

		// keep inside domain
		if ( next[firstDim] < weights.lbound(firstDim) ){
			next[firstDim] = weights.lbound(firstDim);
		}
		else{
			if ( next[firstDim] > weights.ubound(firstDim) ){
				next[firstDim] = weights.ubound(firstDim);
			}
		}

		if ( next[secondDim] < weights.lbound(secondDim) ){
			next[secondDim] = weights.lbound(secondDim);
		}
		else{
			if ( next[secondDim] > weights.ubound(secondDim) ){
				next[secondDim] = weights.ubound(secondDim);
			}
		}

		//if ( this->hasCut & this->isAtCut( next, path[ path.size() - 1 ] ) ){
		//	cout  << path[ path.size() - 1 ] << end;
		//	cout  << next << end;
		//	path.push_back( next );
		//	return true;
		//}

		int nfx = floor(next[firstDim]);
		int ncx = ceil(next[firstDim]);
		int nfy = floor(next[secondDim]);
		int ncy = ceil(next[secondDim]);

		if ( weights( nfx, nfy ) == 0 ){
			path.push_back( next );
			next = TinyVector< int, 2 >(nfx,nfy);
		}
		else{
			if ( weights( nfx, ncy ) == 0 ){
				path.push_back( next );
				next = TinyVector< int, 2 >(nfx,ncy);
			}
			else{
				if ( weights( ncx, nfy ) == 0 ){
					path.push_back( next );
					next = TinyVector< int, 2 >(ncx,nfy);
				}
				else{
					if ( weights( ncx, ncy ) == 0 ){
						path.push_back( next );
						next = TinyVector< int, 2 >(ncx,ncy);
					}
				}
			}
		}

		path.push_back( next );

		// check if got stuck
		if ( path.size() > 150000 ){
			cout << next << endl;
			return false;
		}
	}
}

//template< class Pixel >
//bool nbfGeodesicPath< Pixel> :: getBilinearPath( Array< Pixel, 2 > & weights, 
//									             TinyVector< int, 2 > & start,
//										         vector< TinyVector< Pixel, 2 > > & path )
//{	
//	// reset path if not empty
//	path.clear();
//
//	// store next point in geodesic
//	TinyVector< Pixel, 2 > next;
//	
//	// store current point neighbors
//	vector< TinyVector< Pixel, 2 > > neighbors;
//
//	// push first point into path
//	next = start;
//	path.push_back( next );
//	
//	Pixel minDistance;
//
//	// interpolator to get bilinear gradient estimation
//	nbfLinearInterpolator< Pixel, 2 > interp( weights );
//
//	// this is the time step, we use a fixed value. 
//	// Gradients (this is the incremental update) are normalized to [0,1]
//	// So 'delta' is telling us that we move at most delta in the direction
//	// of the gradient, the grid spacing being 1.
//	Pixel delta = .5;
//
//	// flag for signaling when the backprop gets stuck.
//	bool stucked = false;
//
//	// loop until we reach zero distance
//	while ( weights( (int)next(firstDim), (int)next(secondDim) ) > 0 )
//	{
//		//cout << next << endl;
//		// made last position coordinates handy
//		Pixel x = next[firstDim];
//		Pixel y = next[secondDim];
//
//		// get square box where the point is contained
//		int lx = floor( x ); int ux = ceil( x );
//		int ly = floor( y ); int uy = ceil( y );
//
//		TinyVector< int, 2 > minimo;
//		Pixel minDistance = numeric_limits<Pixel>::max();
//
//		// grid point (look in all neighbors)
//		if ( ( lx == ux ) & ( ly == uy ) ){
//			vector< TinyVector< int, 2 > > neighbors;
//			this->getFullNeighbors( TinyVector< int, 2 >(lx,ly), neighbors );
//			minDistance = weights(neighbors[0]);
//			minimo = neighbors[0];
//			for ( int i = 1; i < neighbors.size(); i++ ){
//				if ( weights( neighbors[i] ) < minDistance ){
//					minDistance = weights( neighbors[i] );
//					minimo = neighbors[i];
//				}
//			}
//		}
//		else{
//			// x grid line (look right and left)
//			if ( lx == ux ){
//				if ( lx < weights.ubound(firstDim) ){
//					p = TinyVector< int, 2 >(lx+1,ly);
//					d = weights(p);
//					if ( d < minDistance ){
//						minimo = p;
//						minDistance = d;
//					}
//					p = TinyVector< int, 2 >(lx+1,uy);
//					d = weights(p);
//					if ( d < minDistance ){
//						minimo = p;
//						minDistance = d;
//					}
//				}
//				if ( lx > weights.lbound(firstDim) ){
//					p = TinyVector< int, 2 >(lx-1,ly);
//					d = weights(p);
//					if ( d < minDistance ){
//						minimo = p;
//						minDistance = d;
//					}
//					p = TinyVector< int, 2 >(lx-1,uy);
//					d = weights(p);
//					if ( d < minDistance ){
//						minimo = p;
//						minDistance = d;
//					}
//				}
//			}
//			else{
//				// y grid line (look up and down)
//				if ( ly == uy ){
//					if ( ly < weights.ubound(secondDim) ){
//						p = TinyVector< int, 2 >(lx,ly+1);
//						d = weights(p);
//						if ( d < minDistance ){
//							minimo = p;
//							minDistance = d;
//						}
//						p = TinyVector< int, 2 >(ux,ly+1);
//						d = weights(p);
//						if ( d < minDistance ){
//							minimo = p;
//							minDistance = d;
//						}
//					}
//					if ( ly > weights.lbound(secondDim) ){
//						p = TinyVector< int, 2 >(lx,ly-1);
//						d = weights(p);
//						if ( d < minDistance ){
//							minimo = p;
//							minDistance = d;
//						}
//						p = TinyVector< int, 2 >(ux,ly-1);
//						d = weights(p);
//						if ( d < minDistance ){
//							minimo = p;
//							minDistance = d;
//						}
//					}
//				}
//				// generic point
//				else{
//					minimo = next;
//				}
//			}
//		}
//
//		// take average between current point and minimum neighbor to get point inside square
//		TinyVector< Pixel, 2 > center = ( minimo + next ) / 2.0;
//
//		Pixel cx = center[firstDim];
//		Pixel cy = center[secondDim];
//		int lcx = floor( cx ); int ucx = ceil( cx );
//		int lcy = floor( cy ); int ucy = ceil( cy );
//
//		TinyVector< Pixel, 2 > gradient;
//
//		if ( lcx == ucx ){
//			gradient[firstDim] = 0;
//			gradient[secondDim] = 1;
//		}
//		else{
//			if ( lcy == ucy ){
//				gradient[firstDim] = 0;
//				gradient[secondDim] = 1;
//			}
//			else{
//				gradient = interp.interpolateSingleGradient( next );
//			}
//		}
//
//		// Switch to upwing gradient on grid lines/points.
//		// when we are close enough to a grid point (1e-3), we switch to this approximation.
//
//		// store position where we got stuck
//		TinyVector< int, 2 > stuckPosition;
//
//}

template< class Pixel >
bool nbfGeodesicPath< Pixel> :: getForwardPath( Array< Pixel, 2 > & input,
											    TinyVector< int, 2 > & start,
												vector< TinyVector< int, 2 > > & path )
{	
	path.clear();

	TinyVector< int, 2 > next;

	next = start;
	path.push_back( next );

	Pixel t1, t2, t3;
	int x1, x2, x3, y1, y2, y3;

	// while distance not zero
	while ( input( next ) != 0 ){

		int x = next(firstDim); 
		int y = next(secondDim);
		int currentTheta = this->theta( next );
		switch ( currentTheta ){

			case 0:
				x1 = x + 1; y1 = y + 1;
				x2 = x;		y2 = y + 1;
				x3 = x - 1; y3 = y + 1;
				break;
			case 1:
				x1 = x;		y1 = y + 1;
				x2 = x - 1;	y2 = y + 1;
				x3 = x - 1;	y3 = y;
				break;
			case 2:
				x1 = x - 1;	y1 = y + 1;
				x2 = x - 1; y2 = y;
				x3 = x - 1;	y3 = y - 1;
				break;
			case 3:
				x1 = x - 1;	y1 = y;
				x2 = x - 1;	y2 = y - 1;
				x3 = x;		y3 = y - 1;
				break;
			case 4:
				x1 = x - 1;	y1 = y - 1;
				x2 = x;		y2 = y - 1;
				x3 = x + 1;	y3 = y - 1;
				break;
			case 5:
				x1 = x;		y1 = y - 1;
				x2 = x + 1;	y2 = y - 1;
				x3 = x + 1;	y3 = y;
				break;
			case 6:
				x1 = x + 1;	y1 = y - 1;
				x2 = x + 1;	y2 = y;
				x3 = x + 1;	y3 = y + 1;
				break;
			case 7:
				x1 = x + 1;	y1 = y;
				x2 = x + 1;	y2 = y + 1;
				x3 = x;		y3 = y + 1;
				break;
		} // switch
		
		if ( input.isInRange(x1,y1) ){
			t1 = input( x1, y1 );
		}
		else {
			t1 = numeric_limits<Pixel>::max();

		}
		if ( input.isInRange(x2,y2) ){
			t2 = input( x2, y2 );
		}
		else{
			t2 = numeric_limits<Pixel>::max();
		}
		if ( input.isInRange(x3,y3) ){
			t3 = input( x3, y3);
		}
		else{
			t3 = numeric_limits<Pixel>::max();
		}

		if ( ( t1 < t2 ) & ( t1 < t3 ) ){
			next(firstDim) = x1;
			next(secondDim) = y1;
		}
		else if ( ( t3 < t1 ) & ( t3 < t2 ) ){
			next(firstDim) = x3;
			next(secondDim) = y3;
		}
		else{
			next(firstDim) = x2;
			next(secondDim) = y2;
		}

		if ( !( input.isInRange(x1,y1) | input.isInRange(x2,y2) | input.isInRange(x3,y3) ) ){
			// ** WARNING ** - path didn't end at distance = 0
			return false;
		}

		path.push_back( next );

		// check if cut reached
		if ( ( next(secondDim) == (start(secondDim) - 1) ) & ( next(firstDim) > center(firstDim) ) ) 
		{
			cout << "start = " << start << endl;
			cout << "next = " << next << endl;
			return true;
		}

		if ( !this->checkPath( path ) ){
			return false;
		}

	}
	return true;
}


template< class Pixel >
bool nbfGeodesicPath< Pixel> :: getBackwardPath( Array< Pixel, 2 > & input,
												 TinyVector< int, 2 > & start,
												 vector< TinyVector< int, 2 > > & path )
{
	// we need these modifications to do the forward path in the adequate direction (couter-clock-wise)

	input.reverseSelf(secondDim);

	// mirror transform center
	this->center(secondDim) = input.ubound(secondDim) - center(secondDim);

	// we changed the center, so we recompute thetas
	this->updateTheta();

	TinyVector< int, 2 > newStart;
	newStart(firstDim) = start(firstDim);
	newStart(secondDim) = input.ubound(secondDim) - start(secondDim);
	
	bool state = this->getForwardPath( input, newStart, path );

	// undo transformation
	this->weight.reverseSelf(secondDim);

	this->center(secondDim) = input.ubound(secondDim) - center(secondDim);
	this->updateTheta();

	for ( int i = 0; i < path.size(); i++ )
	{
		path[i](secondDim) = input.ubound(secondDim) - path[i](secondDim);
	}
	reverse( path.begin(), path.end() );

	return state;
}

template< class Pixel >
Pixel nbfGeodesicPath< Pixel > :: getCircularPath( Array< Pixel, 2 > & weights, vector< TinyVector< Pixel, 2 > > & path )
{
	//nbfMatlabWriter writer;
	//writer.setFileName("ipath.joder");

	nbfFastMarching2D8< Pixel > fm( weights );
	//nbfFastFastMarching2D< Pixel > fm( weights );
	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;

	// set alive (d = 0) all points on the left side of the image
	for ( int i = 0; i < weights.rows(); i++ )
	{
		TinyVector< int, 2 > current( i, 0 );
		if ( weights( current ) < numeric_limits<Pixel>::max() ){
			aliveP.push_back(current);
			aliveD.push_back(0);
		}
	}			

	fm.setAliveSet( aliveP, aliveD );
	Array< Pixel, 2 > distances( weights.shape() );

	// stop when opposite side reached
	// fm.setStopBorder(secondDim);

	// compute distances
	TinyVector< int, 2 > first = fm.execute( distances );

	//writer.write(weights);
	//writer.write(distances);

	//cout << first << endl;

	// new way of computing circular paths

	Array< Pixel, 1 > rightLine( distances.rows() );
	rightLine = distances( Range::all(), distances.ubound(secondDim) );

	// check if path exists, if not, exit
	if ( min( rightLine ) == numeric_limits< Pixel > :: max() ){
		return -1;
	}

	aliveP.clear();
	aliveD.clear();
	for ( int i = 0; i < weights.rows(); i++ )
	{
		TinyVector< int, 2 > current( i, weights.ubound(secondDim) );
		if ( weights( current ) < numeric_limits<Pixel>::max() ){
			aliveP.push_back(current);
			aliveD.push_back(0);
		}
	}			

	fm.setAliveSet( aliveP, aliveD );
	fm.execute(distances);

	rightLine += distances( Range::all(), 0 );
	TinyVector< int, 1 > point = minIndex( rightLine );

	//cout << point << endl;

	aliveP.clear();
	aliveD.clear();
	aliveP.push_back( TinyVector< int, 2>( point[0], 0 ) );
	aliveD.push_back(0);
	fm.setAliveSet( aliveP, aliveD );
	TinyVector< int, 2 > endPoint( point[0], distances.ubound(secondDim) );
	fm.setStopPoint( endPoint );
	TinyVector< int, 2 > d = fm.execute(distances);

	//writer.write(distances);

	this->getSimplePath(distances,endPoint,path);

	//// back-propagation
	//vector< TinyVector< Pixel, 2 > > pathi;
	////this->getForwardLinearPath(distances,first,pathi);
	//this->getSimplePath(distances,first,pathi);

	//// enforce circularity with computed end points
	//if ( pathi.size() > 0 ){
	//	aliveP.clear();
	//	aliveP.push_back( pathi[ pathi.size() - 1 ] );
	//	aliveD.clear();
	//	aliveD.push_back(0);
	//	fm.setAliveSet( aliveP, aliveD );
	//	TinyVector< int, 2 > endPoint( pathi[ pathi.size() - 1 ](firstDim), distances.ubound(secondDim) );
	//	//cout << endPoint << endl;
	//	fm.setStopPoint(endPoint);
	//	fm.execute(distances);

	//	writer.write(distances);
	//	//this->getForwardLinearPath(distances,endPoint,path);
	//	this->getSimplePath(distances,endPoint,path);
	//}

	return distances(d);
}


template< class Pixel >
void nbfGeodesicPath< Pixel > :: getFirstPointOnRightSide( Array< Pixel, 2 > & weights, TinyVector< Pixel,2 > & start, TinyVector< Pixel, 2 > & end )
{
	nbfFastMarching2D8< Pixel > fm( weights );
	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;

	// d = 0 on the left side
	for ( int i = 0; i < weights.rows(); i++ )
	{
		TinyVector< int, 2 > current( i, 0 );
		if ( weights( current ) < numeric_limits<Pixel>::max() ){
			aliveP.push_back(current);
			aliveD.push_back(0);
		}
	}			

	fm.setAliveSet( aliveP, aliveD );
	Array< Pixel, 2 > distances( weights.shape() );
	fm.setStopBorder(secondDim);
	end = fm.execute( distances );
	vector< TinyVector< Pixel, 2 > > pathi;
	TinyVector< int, 2 > stopFM = floor( end );
	this->getSimplePath(distances,stopFM,pathi);
	start = pathi[ pathi.size() - 1 ];
}


template< class Pixel >
bool nbfGeodesicPath< Pixel > :: getPathBetweenSides( Array< Pixel, 2 > & weights, vector< TinyVector< Pixel, 2 > > & path )
{
	//nbfMatlabWriter writer;
	//writer.setFileName("ipath");

	nbfFastMarching2D8< Pixel > fm( weights );
	vector< TinyVector< int, 2 > > aliveP;
	vector< Pixel > aliveD;

	// d = 0 on the left side
	for ( int i = 0; i < weights.rows(); i++ )
	{
		TinyVector< int, 2 > current( i, weights.lbound(secondDim) );
		if ( weights( current ) < numeric_limits<Pixel>::max() ){
			aliveP.push_back(current);
			aliveD.push_back(0);
		}
	}			

	fm.setAliveSet( aliveP, aliveD );
	Array< Pixel, 2 > distances( weights.shape() );
	fm.setStopBorder(secondDim);
	TinyVector< int, 2 > first = fm.execute( distances );

	this->getSimplePath(distances,first,path);

	//writer.write(weights);
	//writer.write(distances);

	//aliveP.clear();
	//aliveP.push_back( first );
	//aliveD.clear();
	//aliveD.push_back(0);
	//fm.setAliveSet( aliveP, aliveD );
	//fm.setStopBorder(secondDim);
	//TinyVector< int, 2 > last = fm.execute(distances);

	//writer.write(distances);

	////this->getForwardLinearPath(distances,last,path);
	//this->getSimplePath(distances,last,path);
	//Array< Pixel, 2 > ipath( distances.shape() );
	//this->getImplicitPath(path,ipath);
	//writer.write(ipath);

	return true;
}

template< class Pixel >
bool nbfGeodesicPath< Pixel > :: checkPath( vector< TinyVector< int, 2 > > & path )
{
	// check if too many points
	if ( path.size() > 500 ){
		return false;
	}
	return true;
}

template< class Pixel >
void nbfGeodesicPath< Pixel > :: getImplicitPath( vector< TinyVector< int, 2 > > & path,
												  Array< Pixel, 2 > & A )
{
	if ( path.size() == 0 ){
		A = numeric_limits<Pixel>::max();
	}
	else{
		// if open path then switch to simple point representation
		if ( ( path[0](firstDim) - path[ path.size() - 1 ](firstDim) > 1 ) |
             ( path[0](secondDim) - path[ path.size() - 1 ](secondDim) > 1 ) )
		{
			A = -1;
			for ( int i = 0; i < path.size(); i++ )
			{
				A(path[i]) = 1;
			}
		}
		else
		{
			Array< Pixel, 2 > D( A.shape() );

			// compute iso-distances to path first
			A = 1;
			//nbfFastFastMarching2D< Pixel > fm2d(A);
			nbfFastMarching2D< Pixel > fm2d(A);

			vector< TinyVector< int, 2 > > positions;
			vector< Pixel > distances;
			for ( int i = 0; i < path.size(); i++ ){
				positions.push_back( path[i] );
				distances.push_back( 0 );
			}
			fm2d.setAliveSet(positions,distances);
			fm2d.setStopDistance(10);
			fm2d.execute(D);

			// now compute inside/outside
			A[path] = numeric_limits<Pixel>::max();
			TinyVector< int, 2 > alive( this->center(firstDim), this->center(secondDim) );
			nbfFastMarchingFool2D< Pixel > fm(A);
			fm.setAliveSet( alive, 0 );
			Array< Pixel, 2 > S( A.shape() );
			fm.execute(S);
			// if not sufficient points inside
			if ( sum( S == numeric_limits<Pixel>::max() ) < 5 )
			{
				A = 1;
			}
			else
			{
				A = where( S < numeric_limits<Pixel>::max(), -D, D );
			}
		}
	}
}

template< class Pixel >
void nbfGeodesicPath< Pixel > :: getImplicitOpenPath( vector< TinyVector< int, 2 > > & path,
													  Array< Pixel, 2 > & A )
{
	A = 10 * blitz::minmax::max( A.rows(), A.cols() );

	for ( int i = 0; i < path.size(); i++ ){
		Pixel x = path[i](firstDim);
		Pixel y = path[i](secondDim);

		firstIndex fi;
		A( Range( fromStart, x), y ) = where( A( Range( fromStart, x), y ) > fi - x, fi - x, A( Range( fromStart, x), y ) );
		A( Range( x + 1, toEnd ), y ) = where( A( Range( x + 1, toEnd ), y ) > fi, fi + 1, A( Range( x + 1, toEnd ), y ) );
	}
}

//template< class Pixel >
//void nbfGeodesicPath< Pixel > :: getImplicitPath( vector< TinyVector< Pixel, 2 > > & path,
//												  Array< Pixel, 2 > & A )
//{
//	A = numeric_limits<Pixel>::max();
//
//	Pixel radius = 5;
//
//	int lboundX, uboundX, lboundY, uboundY;
//
//	for ( int i = 0; i < path.size(); i++ )
//	{
//		Pixel cX = path[i](firstDim);
//		Pixel cY = path[i](secondDim);
//		lboundX = blitz::minmax::max( A.lbound(firstDim), floor(cX) - radius );
//		uboundX = blitz::minmax::min( A.ubound(firstDim), ceil(cX) + radius );
//		lboundY = blitz::minmax::max( A.lbound(secondDim), floor(cY) - radius );
//		uboundY = blitz::minmax::min( A.ubound(secondDim), ceil(cY) + radius );
//		Array< Pixel, 2 > Ai( A( Range(lboundX,uboundX), Range(lboundY,uboundY) ) );
//		Array< Pixel, 2 > Di( Ai.shape() );
//		firstIndex indexX;
//		secondIndex indexY;
//		Di = sqrt( pow2(indexX-cX+lboundX) + pow2(indexY-cY+lboundY) );
//		Ai = where( Di < Ai, Di, Ai );
//	}
//}


template< class Pixel >
void nbfGeodesicPath< Pixel > :: getImplicitPath( vector< TinyVector< Pixel, 2 > > & path,
												  Array< Pixel, 2 > & A,
												  Pixel maxDistance )
{
	if ( path.size() == 0 ){
		A = numeric_limits<Pixel>::max();
	}
	else{
		Array< Pixel, 2 > D( A.shape() );
		D = numeric_limits<Pixel>::max();

		Array< Pixel, 2 > S( A.shape() );
		S = 1;

		TinyVector< Pixel, 2 > position;
		
		TinyVector< Pixel, 2 > gradient;

		// build inside/outside regions
		for ( int i = 0; i < path.size(); i++ ){

			Pixel cX = path[i](firstDim);
			Pixel cY = path[i](secondDim);

			int lx = floor( cX ); int ux = ceil( cX );
			int ly = floor( cY ); int uy = ceil( cY );
	
			if ( i < path.size() - 1 ){
				gradient = path[i+1] - path[i];
			}

			if ( D.isInRange(lx,ly) ){
				D( lx, ly ) = min( D( lx, ly ), sqrt( pow2(cX-lx) + pow2(cY-ly) ) );
				// handle this appart
				if ( ( lx == ux ) | ( ly == uy ) ){
					S(lx,ly) = numeric_limits<Pixel>::max();
				}
				else{
					position = TinyVector< Pixel, 2 >(lx,ly) - path[i];
					if ( gradient[1] * position[0] - gradient[0] * position[1] < 0 ){
						S( lx, ly ) = numeric_limits<Pixel>::max();
					}
				}
			}
			if ( D.isInRange(lx,uy) ){
				D( lx, uy ) = min( D( lx, uy ), sqrt( pow2(cX-lx) + pow2(cY-uy) ) );
				position = TinyVector< Pixel, 2 >(lx,uy) - path[i];
				if ( gradient[1] * position[0] - gradient[0] * position[1] < 0){
					S( lx, uy ) = numeric_limits<Pixel>::max();
				}
			}
			if ( D.isInRange(ux,ly) ){
				D( ux, ly ) = min( D( ux, ly ), sqrt( pow2(cX-ux) + pow2(cY-ly) ) );
				position = TinyVector< Pixel, 2 >(ux,ly) - path[i];
				if ( gradient[1] * position[0] - gradient[0] * position[1] < 0 ){
					S( ux, ly ) = numeric_limits<Pixel>::max();
				}
			}
			if ( D.isInRange(ux,uy) ){
				D( ux, uy ) = min( D( ux, uy ), sqrt( pow2(cX-ux) + pow2(cY-uy) ) );
				position = TinyVector< Pixel, 2 >(ux,uy) - path[i];
				if ( gradient[1] * position[0] - gradient[0] * position[1] < 0 ){
					S( ux, uy ) = numeric_limits<Pixel>::max();
				}
			}
			
			//if ( ( lx == ux ) | ( ly == uy ) ){
			//	S(lx,ly) = numeric_limits<Pixel>::max();
			//}
			//else{
			//	if ( i < path.size() - 1 ){
			//		TinyVector< Pixel, 2 > gradient( path[i+1] - path[i] );
			//		//cout << gradient << endl;
			//		position = TinyVector< Pixel, 2 >(lx,ly) - path[i];
			//		//cout << position << endl;
			//		if ( gradient[1] * position[0] - gradient[0] * position[1] < 0 ){
			//			S( lx, ly ) = numeric_limits<Pixel>::max();
			//		}
			//		position = TinyVector< Pixel, 2 >(lx,uy) - path[i];
			//		//cout << position << endl;
			//		if ( gradient[1] * position[0] - gradient[0] * position[1] < 0){
			//			S( lx, uy ) = numeric_limits<Pixel>::max();
			//		}
			//		position = TinyVector< Pixel, 2 >(ux,ly) - path[i];
			//		//cout << position << endl;
			//		if ( gradient[1] * position[0] - gradient[0] * position[1] < 0 ){
			//			S( ux, ly ) = numeric_limits<Pixel>::max();
			//		}
			//		position = TinyVector< Pixel, 2 >(ux,uy) - path[i];
			//		//cout << position << endl;
			//		if ( gradient[1] * position[0] - gradient[0] * position[1] < 0 ){
			//			S( ux, uy ) = numeric_limits<Pixel>::max();
			//		}
			//	}
			//}
		}

		if ( maxDistance > 1 ){
			// compute iso-distances to path first
			A = 1;
			nbfFastMarching2D< Pixel > fm2d(A);
			//nbfFastFastMarching2D< Pixel > fm2d(A);

			vector< TinyVector< int, 2 > > positions;
			vector< Pixel > distances;

			Array< Pixel, 2 > :: iterator iter = D.begin();
			while( iter != D.end() ){
				if ( (*iter) < numeric_limits<Pixel>::max() ){
					positions.push_back( iter.position() );
					distances.push_back( *iter );
				}
				++iter;
			}
			fm2d.setAliveSet(positions,distances);
			fm2d.setStopDistance(maxDistance);
			fm2d.execute(D);
		}

		// now compute inside/outside

		TinyVector< int, 2 > alive( 0, 0 );
		nbfFastMarchingFool2D< Pixel > fm(S);
		fm.setAliveSet( alive, 0 );
		Array< Pixel, 2 > sign( A.shape() );
		fm.execute(sign);

		A = where( sign < numeric_limits<Pixel>::max(), D, -D );
	}	
}


template< class Pixel >
void nbfGeodesicPath< Pixel> :: getForwardLinearPath( Array< Pixel, 2 > & distance,
												      TinyVector< int, 2 > & start,
												      vector< TinyVector< int, 2 > > & path )
{	
	path.clear();

	TinyVector< int, 2 > next;

	// back propagation - minimum distance neighbor
	next = start;
	path.push_back( next );

	int x, x1, x2, x3;
	int y, y1, y2, y3;
	//int x4, x5;
	//int y4, y5;
	while ( next(secondDim) != 0 ){

			Pixel t1, t2, t3;
			//Pixel t4, t5;

			x = next(firstDim); 
			y = next(secondDim);

			x1 = x - 1; y1 = y - 1;
			x2 = x;		y2 = y - 1;
			x3 = x + 1; y3 = y - 1;

			//x4 = x - 1; y4 = y;
			//x5 = x + 1; y5 = y;

			if ( distance.isInRange(x1,y1) ){
				t1 = distance( x1, y1 );
			}
			else {
				t1 = numeric_limits<Pixel>::max();
			}
			if ( distance.isInRange(x2,y2) ){
				t2 = distance( x2, y2 );
			}
			else{
				t2 = numeric_limits<Pixel>::max();
			}
			if ( distance.isInRange(x3,y3) ){
				t3 = distance( x3, y3);
			}
			else{
				t3 = numeric_limits<Pixel>::max();
			}
			//if ( distance.isInRange(x4,y4) ){
			//	t4 = distance( x4, y4);
			//}
			//else{
			//	t4 = numeric_limits<Pixel>::max();
			//}
			//if ( distance.isInRange(x5,y5) ){
			//	t5 = distance( x5, y5);
			//}
			//else{
			//	t5 = numeric_limits<Pixel>::max();
			//}

			//if ( ( t1 < t2 ) & ( t1 < t3 ) & ( t1 < t4 ) & ( t1 < t5 ) ){
			if ( ( t1 < t2 ) & ( t1 < t3 ) ){
				next(firstDim) = x1;
				next(secondDim) = y1;
			}
			//else if ( ( t2 < t1 ) & ( t2 < t3 ) & ( t2 < t5 ) & ( t2 < t5 ) ){
			else if ( ( t3 < t1 ) & ( t3 < t2 ) ){
				next(firstDim) = x3;
				next(secondDim) = y3;
			}
			//else if ( ( t3 < t1 ) & ( t3 < t2 ) & ( t3 < t4 ) & ( t3 < t5 ) ){
			else{
				next(firstDim) = x2;
				next(secondDim) = y2;
			}
			//else if ( ( t4 < t1 ) & ( t4 < t2 ) & ( t4 < t3 ) & ( t4 < t5 ) ){
			//	next(firstDim) = x4;
			//	next(secondDim) = y4;
			//}
			//else if ( ( t5 < t1 ) & ( t5 < t2 ) & ( t5 < t3 ) & ( t5 < t4 ) ){
			//	next(firstDim) = x5;
			//	next(secondDim) = y5;
			//}
			//else{
			//	next(firstDim) = next(firstDim);
			//	next(secondDim) = next(secondDim) - 1;
			//}

			path.push_back( next );
	}

	// add last point
	path.push_back( next );
}

#endif /* FILE_nbfGeodesicPath */