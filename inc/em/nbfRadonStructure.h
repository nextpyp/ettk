#pragma once

#include <em/nbfRadonBlock.h>
#include <random/discrete-uniform.h>

/** Top level data structure to handle Radon projections and backprojections.
	This class provides a framework so all iterative algorithms ART-type, 
	SIRT-type, etc. have a trivial implementation.
	A data structure is assembled that contains coeficients (matrix A) and references
	to pixel values (vector x) in which all computations are based. This structure
	has to be computed only once, since A remains fixed throughout the process.
*/
template< class Pixel >
class nbfRadonStructure
{
public:

	// Default constructor
	nbfRadonStructure(){ this->geometryReady = false; };

	~nbfRadonStructure(){};

	// Set output image (where resulting image is stored)
	void setImage( Array< Pixel, 2 > & );

	// Set target projections.
	void setProjections( Array< Pixel, 2 > & );

	// Set projection angles.
	// WARNING:
	// For +/-45 angles the projection operation becomes ill posed resulting in projection artifacts.
	// The main reason why this happens is because 2 of the 4 points representing a single pixel
	// have coincident contributions.
	// To remedy this, we slightly change the projection angle closest to 45 degrees.
	void setAngles( Array< Pixel, 1 > & );

	// ART iterative technique. 
	// Input image is used as initial state and number of iterations is specified in argument.
	void art(int);

	// SIRT iterative technique.
	// Input image is used as initial state and number of iterations is specified in argument.
	void sirt(int);

protected:

	// Build projection framework, equivalent to build the matrix A in the system Ax=b.
	void buildRadonFramework();

	// Project and accumulate all blocks: this->outImage += \sum_i P_i( this->inImage )
	void project();

	// Project single block: this->outImage += P_i( this->inImage )
	void project(int);

	// Get projection at given angle
	void getProjections( Array< Pixel, 2 > & );

	// Backproject and accumulate all blocks: this->outImage += \sum_i BP_i( offset )
	void backProject();

	// Backproject and accumulate single block: this->outImage += BP_i( offset )
	void backProject(int);

	// Flag to track the status of the data structure.
	bool geometryReady;

	// List of equation blocks.
	vector< nbfRadonBlock< Pixel > > blocks;

	// Store projection angles.
	Array< Pixel, 1 > angles;

	// Store image (i.e. vector x)
	Array< Pixel, 2 > inImage;

	// Store output image
	Array< Pixel, 2 > outImage;

	// Store image projections
	Array< Pixel, 2 > projections;

	// Store normalization matrix for computed projections
	Array< Pixel, 2 > circle;
};

template< class Pixel >
void nbfRadonStructure< Pixel > :: setImage( Array< Pixel, 2 > & image )
{
	this->inImage.reference( image );
	this->outImage.resize( this->inImage.shape() );
}


template< class Pixel >
void nbfRadonStructure< Pixel > :: setAngles( Array< Pixel, 1 > & angles )
{
	this->angles.resize( angles.shape() );
	this->angles = angles;
	//this->angles.reference(angles);
}

template< class Pixel >
void nbfRadonStructure< Pixel > :: getProjections( Array< Pixel, 2 > & proj )
{
	typedef typename vector< nbfRadonBlock< Pixel > > :: iterator blockIterator;
	blockIterator iBlock = this->blocks.begin();
	int i = 0;
	while ( iBlock != this->blocks.end() ){
		typedef typename vector< nbfRadonEquation< Pixel > > :: iterator equationIterator;
		equationIterator iEquation = iBlock->equations.begin();
		int j = 0;
		while ( iEquation != iBlock->equations.end() ){
			proj(j,i) = iEquation->proj;
			iEquation++;
			j++;
		}
		iBlock++;
		i++;
	}
}

template< class Pixel >
void nbfRadonStructure< Pixel > :: setProjections( Array< Pixel, 2 > & in )
{
	this->projections.reference( in );

	//// normalize projection values
	//Pixel R = ( this->inImage.rows() - 1 ) / 2.0;
	//this->circle.resize( this->inImage.rows() );
	//firstIndex i;
	//this->circle = ( 2 * sqrt( R*R - pow2( R - i + 0.0 ) ) + 1 ) / (float)this->inImage.rows(); 
	//for ( int i = 0; i < this->projections.cols(); i++ ){
	//	this->projections( Range::all(), i ) *= this->circle;
	//	// correct for infinite plane specimen
	//	//this->projections( Range::all(), i ) *= ( circle * cos( this->angles(i) ) );
	//}

	// assign b values to equations if data structure already in place
	for ( unsigned int i = 0; i < this->blocks.size(); i++ ){
		Array< Pixel, 1 > dummy( this->projections( Range::all(), i ) );
		this->blocks[i].setProjections( dummy );
	}
}

template< class Pixel >
void nbfRadonStructure< Pixel > :: buildRadonFramework()
{
	// calculate image center
	TinyVector< Pixel, 2 > center = ( this->inImage.shape() - 1 ) / 2.0;

	// Calculate number of useful equations per block (i.e. for each angle).
	// This is done by computing the distance of a corner point to the center.
	Array< int, 1 > projectionBounds( this->angles.size() );
	
	//// Compute usable bounds (the next line is kept as reference, a simplyfied version is actually used):
	//// TinyVector< int, 2 > c1( this->inImage.lbound(firstDim), this->inImage.ubound(secondDim) );
	//// projectionBounds = ceil( 2 * ( center(firstDim) - c1(firstDim) ) * sin( abs(this->angles) ) + ( center(secondDim) - c1(secondDim) ) * cos( abs(this->angles) );
	//projectionBounds = ceil( 2 * ( center(firstDim) * sin( abs(this->angles) ) - center(secondDim) * cos( abs(this->angles) ) ) - 1 );

	//// eliminate all angles that include insuficient information
	//projectionBounds = where( fabs( this->angles ) > atan2( center[1], center[0] ), 0, abs( projectionBounds ) );

	projectionBounds = this->projections.rows();
	//cout << projectionBounds << endl;

	// create blocks of equations (one for each projection angle) with only the neccesary number of equations.
	for ( int i = 0; i < this->angles.size(); i++ ){
		Array< Pixel, 1 > dummy( this->projections( Range::all(), i ) );
		this->blocks.push_back( nbfRadonBlock< Pixel >( projectionBounds(i), dummy ) );
	}

	Array< Pixel, 1 >   sines( sin( this->angles ) );
	Array< Pixel, 1 > cosines( cos( this->angles ) );

	// iterate through all pixels in input image
	typedef typename Array< Pixel, 2 > :: iterator arrayIterator2D;
	arrayIterator2D iImage = this->inImage.begin();
	arrayIterator2D oImage = this->outImage.begin();
	while ( iImage != this->inImage.end() ){

		// for all projection angles

		Pixel offX = iImage.position()[0] - center[0];
		Pixel offY = iImage.position()[1] - center[1];

		// only consider points in central circle
		if ( sqrt( offX * offX + offY * offY ) > center[1] ){
			++iImage;
			++oImage;
			continue;
		}

		//if ( ( abs( offX ) > 256 ) || ( abs( offY ) > 256 ) ){
		//	++iImage;
		//	++oImage;
		//	continue;
		//}

		nbfRadonBlock< Pixel > * pBlock;

		for ( int i = 0; i < this->angles.size(); i++ ){

			//bool specialPositive45 = false;
			//Pixel epsilon = .25 * vtkMath::DegreesToRadians();
			//// set flag if current angle is close to 45 degrees
			//if ( fabs( angles(i) - 45 * vtkMath::DegreesToRadians() ) < epsilon ){
			//	specialPositive45 = true;
			//}

			//bool specialNegative45 = false;
			//// set flag if current angle is close to -45 degrees
			//if ( fabs( angles(i) + 45 * vtkMath::DegreesToRadians()  ) < epsilon ){
			//	specialNegative45 = true;
			//}

			//// ignore these cases
			//if ( ( specialPositive45 == true ) || ( specialNegative45 == true ) ){
			//	continue;
			//}

			//// define epsilon to separate points in ill-posed case
			//Pixel epsilon = gEpsilon / 10.0;

			Pixel sine   =   sines(i);
			Pixel cosine = cosines(i);
			
			Pixel sinpos =   sine * offX + ( projectionBounds(i) - 1.0 ) / 2.0;
			Pixel cospos = cosine * offY;
			
			// create current pixel contribution element in three consecutive equations
			Pixel distanceCenter = sinpos + cospos;
			int fCenter = ceil(distanceCenter);
			if ( fCenter - distanceCenter > .5 ){
				fCenter--;
			}

			if ( ( fCenter > -2 ) && ( fCenter <= projectionBounds(i) ) ){
				pBlock = &(this->blocks[i]);
				nbfRadonElement< Pixel > p( &(*iImage), &(*oImage) );
				if ( ( fCenter > -1 ) && ( fCenter < projectionBounds(i) ) ){
					pBlock->equations[ fCenter     ].elements.push_back( p );
				}
				if ( fCenter > 0 ){
					pBlock->equations[ fCenter - 1 ].elements.push_back( p );
				}
				if ( fCenter < projectionBounds(i) - 1 ){
					pBlock->equations[ fCenter + 1 ].elements.push_back( p );
				}
			}
			else{
				continue;
			}

			sinpos += .25 * sine;
			cospos += .25 * cosine;
			Pixel sinneg = sinpos - .5 * sine;
			Pixel cosneg = cospos - .5 * cosine;

			Pixel distance = sinpos + cospos;

			int fDistance = distance;

			Pixel delta = .25 * ( distance - fDistance );

			// set cumulative weights for both positions
			if ( ( fDistance > -1 ) && ( fDistance < projectionBounds(i) ) ){
				pBlock->equations[ fDistance ].elements.back().incrementWeight( .25 - delta );
			}
			if ( ( fDistance > -2 ) && ( fDistance < projectionBounds(i) - 1 ) ){
				pBlock->equations[ fDistance + 1 ].elements.back().incrementWeight( delta );
			}

			distance = sinpos + cosneg;
			fDistance = distance;
			delta = .25 * ( distance - fDistance );

			// set cumulative weights for both positions
			if ( ( fDistance > -1 ) && ( fDistance < projectionBounds(i) ) ){
				pBlock->equations[ fDistance ].elements.back().incrementWeight( .25 - delta );
			}
			if ( ( fDistance > -2 ) && ( fDistance < projectionBounds(i) - 1 ) ){
				pBlock->equations[ fDistance + 1 ].elements.back().incrementWeight( delta );
			}

			distance = sinneg + cospos;
			fDistance = distance;
			delta = .25 * ( distance - fDistance );

			//if ( specialPositive45 == true ){
			//	// displace point in the positive direction
			//	delta += epsilon;
			//}

			// set cumulative weights for both positions
			if ( ( fDistance > -1 ) && ( fDistance < projectionBounds(i) ) ){
				pBlock->equations[ fDistance ].elements.back().incrementWeight( .25 - delta );
			}
			if ( ( fDistance > -2 ) && ( fDistance < projectionBounds(i) - 1 ) ){
				pBlock->equations[ fDistance + 1 ].elements.back().incrementWeight( delta );
			}

			distance = sinneg + cosneg;
			fDistance = distance;
			delta = .25 * ( distance - fDistance );

			//if ( specialNegative45 == true ){
			//	// displace point in the positive direction
			//	delta += epsilon;
			//}

			// set cumulative weights for both positions
			if ( ( fDistance > -1 ) && ( fDistance < projectionBounds(i) ) ){
				pBlock->equations[ fDistance ].elements.back().incrementWeight( .25 - delta );
			}
			if ( ( fDistance > -2 ) && ( fDistance < projectionBounds(i) - 1 ) ){
				pBlock->equations[ fDistance + 1 ].elements.back().incrementWeight( delta );
			}

			// check for dummy points (w==0)
		    nbfRadonElement< Pixel > * element;
			//if ( ( fCenter > -1 ) && ( fCenter < projectionBounds(i) ) ){
			//	element = &pBlock->equations[ fCenter     ].elements.back();
			//	if ( element->w == 0 ){
			//		pBlock->equations[ fCenter     ].elements.pop_back();
			//	}
			//}
			if ( fCenter > 0 ){
				element = &pBlock->equations[ fCenter - 1 ].elements.back();
				if ( element->w == 0 ){
					pBlock->equations[ fCenter - 1 ].elements.pop_back();
				}
			}
			if ( fCenter < projectionBounds(i) - 1 ){
				element = &pBlock->equations[ fCenter + 1 ].elements.back();
				if ( element->w == 0 ){
					pBlock->equations[ fCenter + 1 ].elements.pop_back();
				}
			}

		}
		++iImage;
		++oImage;
	}

	// Compute normalization along rows of matrix A.
	// Various options may be available:
	// - As if we were to solve the system Ax=b.
	// - Use a value proportional to thickness of tilted specimen
	// - Use a constant value.
	//this->outImage = 0;

	// initialize normalization matrix for real projections
	this->circle.resize( this->blocks[0].equations.size(), this->blocks.size() );

	for ( int i = 0; i < this->blocks.size(); i++ ){
		nbfRadonBlock< Pixel > * pBlock = &( this->blocks[i] );
		int bSize = pBlock->equations.size();
		for ( int j = 0; j < bSize; j++ ){
			Pixel norm = 0;
			nbfRadonEquation< Pixel > * pEquation = &( pBlock->equations[j] );
			Pixel eSize = pEquation->elements.size();
			for ( int k = 0; k < eSize; k++ ){
				//norm += pow2( pEquation->elements[k].w );
				norm += pEquation->elements[k].w;
				//norm += 1;
			}
			//pEquation->norm = sqrt(norm);
			pEquation->norm = norm;

			// store in normalization matrix
			this->circle( j, i ) = norm;

			//pEquation->norm = 1.0 * (float)this->inImage.rows() / cos( this->angles(i) );
			//pEquation->norm = 1.0 * (float)this->inImage.rows();
		}
	}
}

template< class Pixel >
void nbfRadonStructure< Pixel > :: project()
{
	typedef typename vector< nbfRadonBlock< Pixel > > :: iterator blockIterator;
	blockIterator iBlock = this->blocks.begin();
	while ( iBlock != this->blocks.end() ){
		iBlock->project();
		iBlock++;
	}
}

template< class Pixel >
void nbfRadonStructure< Pixel > :: project( int angle )
{
	this->blocks[ angle ].project();
}

template< class Pixel >
void nbfRadonStructure< Pixel > :: backProject()
{
	this->outImage = 0;

	typedef typename vector< nbfRadonBlock< Pixel > > :: iterator blockIterator;
	blockIterator iBlock = this->blocks.begin();
	while ( iBlock != this->blocks.end() ){
		iBlock->backProject();
		iBlock++;
	}
}

template< class Pixel >
void nbfRadonStructure< Pixel > :: art( int maxIters )
{
	// merge input and output image in single reference (for computational efficiency)
	this->outImage.reference( this->inImage );

	// build geometry if neccesary
	if ( this->geometryReady != true ){
		this->buildRadonFramework();
		this->geometryReady = true;
	}

	// Sequential order ART
	for ( int iter = 0; iter < maxIters; iter++ ){
		for ( int block = 0; block < this->blocks.size(); block++ ){		
		//for ( int block = this->blocks.size() - 1; block >= 0; block-- ){		
			// compute forward projection at current angle
			this->blocks[block].project();
			this->blocks[block].backProject( 1.0 );
		}
	//	//this->inImage = where( this->inImage > 2, 2, this->inImage );
	//	//this->inImage = where( this->inImage < 0, 0, this->inImage );
	//	cout << "iter = " << iter << endl;
	//	//Pixel m = mean( this->inImage );
	//	//Pixel std = 2*m;
	//	//cout << "[" << min( this->inImage ) << "," << max(this->inImage) << "], m = " << mean(this->inImage) << endl;
	//	//this->inImage = where( this->inImage > m + std, m + std, this->inImage );
	//	////this->inImage = where( this->inImage < 0, 0, this->inImage );
	}

	//// Random order ART
	//for ( int iter = 0; iter < maxIters; iter++ ){
	//	// build random ordered indexes
	//	vector< int > indexes( this->blocks.size() );
	//	for ( int i = 0; i < indexes.size(); i++ ){
	//		indexes[i] = i;
	//	}
	//	// traverse index vector
	//	while ( indexes.size() > 0 ){
	//		ranlib::DiscreteUniform<int> randomGen( indexes.size() );
	//		randomGen.seed((unsigned int)time(0));
	//		int next = randomGen.random();
	//		vector< int > :: iterator iter = indexes.begin();
	//		iter = iter + next;
	//		this->blocks[*iter].project();
	//		//this->blocks[*iter].backProject(.25);	
	//		this->blocks[*iter].backProject(1.0);	
	//		indexes.erase(iter);
	//	}
	//}
}

template< class Pixel >
void nbfRadonStructure< Pixel > :: sirt( int maxIters )
{
	// build geometry if neccesary
	if ( this->geometryReady != true ){
		this->buildRadonFramework();
		this->geometryReady = true;
	}

	//Array< Pixel, 2 > V( this->outImage.shape() );
	//BordStrategyMirrorDouble< Pixel, 2 > bsForV( V, 1 );
	//float lambda = .0;

	for ( int iter = 0; iter < maxIters; iter++ ){

		this->project();
		this->backProject();

		//V = this->inImage;
		//bsForV.refresh();
		//this->inImage += 3.0 / this->blocks.size() * this->outImage - lambda * curvature2D(V);

		this->inImage += 2.0 / this->blocks.size() * this->outImage;

		// Positivity constraint
		//this->inImage = where( this->inImage < 0, 0, this->inImage );
		
		//// Density constraint
		// Pixel density = sum( this->projections( floor( ( this->projections.rows() - 1 ) / 2.0 ), Range::all() ) );
		// this->inImage = this->inImage * density / sum( this->inImage );
	}
}