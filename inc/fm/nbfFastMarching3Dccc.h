#pragma once

// Class nbfFastMarching3Dccc.
//
// Implements 3D fast marching method with arbitrary weights.
// 

#include <fm/nbfFastMarching3D.h>

template< class Pixel >
class nbfFastMarching3Dccc : public nbfFastMarching3D< Pixel >
{
public:

	// constructor takes weight array as input
	nbfFastMarching3Dccc( Array< Pixel, 3 > & );

	/** Set template image to compute CCC function at each pixel.
	*/
	void setTemplate( Array< Pixel, 3 > & );

	// apply template and return CCC volume
	void applyTemplate( Array< Pixel, 3 > & );

protected:

	/** Handle weight computation on demand.
		When there is no need to compute distances in the entire domain
		and the weight function is computationally expensive to get, this
		method provides significantly computational savings, by only computing
		weight values as they are needed.
	*/
	Pixel getWeight( TinyVector< int, 3 > & );

	Array< Pixel, 3 > templateImage;
	int maxX;
	int maxY;
	int maxZ;
	int templateX;
	int templateY;
	int templateZ;

	Array< Pixel, 3 > weightsOnFly;
};


template< class Pixel >
nbfFastMarching3Dccc< Pixel > :: nbfFastMarching3Dccc( Array< Pixel, 3 > & weight )
: nbfFastMarching3D< Pixel >( weight )
{
	this->weightsOnFly.resize( weight.shape() );
	this->weightsOnFly = where( weight == numeric_limits< Pixel > :: max(), numeric_limits< Pixel > :: max(), 0 );
}

template< class Pixel >
void nbfFastMarching3Dccc< Pixel > :: applyTemplate( Array< Pixel, 3 > & p )
{
	p.resize( this->weight.shape() );
	Array< Pixel, 3 > :: iterator iter = this->weight.begin();
	Array< Pixel, 3 > :: iterator iterP = p.begin();
	while( iter != this->weight.end() ){
		TinyVector< int, 3 > pos = iter.position();
		(*iterP) = this->getWeight( pos );
		++iter;
		++iterP;
	}
}

template< class Pixel >
Pixel nbfFastMarching3Dccc< Pixel > :: getWeight( TinyVector< int, 3 > & p )
{
	// compute weight only if neccesary
	if ( this->weightsOnFly( p ) == 0 ){
		
		int minXw = max( p(firstDim) - templateX, 0 );
		int maxXw = min( p(firstDim) + templateX, maxX );
		Range I( minXw, maxXw );
		
		int minYw = max( p(secondDim) - templateY, 0 );
		int maxYw = min( p(secondDim) + templateY, maxY );
		Range J( minYw, maxYw );
		
		int minZw = max( p(thirdDim) - templateZ, 0 );
		int maxZw = min( p(thirdDim) + templateZ, maxZ );
		Range K( minZw, maxZw );

		Array< Pixel, 3 > wv( this->weight( I, J, K ) );
		
		// make copy
		Array< Pixel, 3 > I1( wv.shape() );
		I1 = wv;

		Array< Pixel, 3 > mask( I1.shape() );
		mask = where( I1 < numeric_limits< Pixel > :: max(), 1, 0 );
		Pixel elements = sum( mask );
		I1 = where( I1 < numeric_limits< Pixel > :: max(), I1, 0.0 );
		I1 = I1 - sum(I1) / elements;
		Pixel norm1 = sqrt( 1.0 * sum( pow2( I1 ) ) );

		int r = I1.ubound(firstDim);
		int c = I1.ubound(secondDim);
		int d = I1.ubound(thirdDim);

		Range It;
		if ( minXw == 0 ){
			It = Range(this->templateImage.ubound(firstDim)-r,this->templateImage.ubound(firstDim));
		}
		else if ( maxXw == maxX ){
			It = Range(0,r);
		}
		else{
			It = Range(0,this->templateImage.ubound(firstDim));
		}

		Range Jt;
		if ( minYw == 0 ){
			Jt = Range(this->templateImage.ubound(secondDim)-c,this->templateImage.ubound(secondDim));
		}
		else if ( maxYw == maxY ){
			Jt = Range(0,c);
		}
		else{
			Jt = Range(0,this->templateImage.ubound(secondDim));
		}

		Range Kt;
		if ( minZw == 0 ){
			Kt = Range(this->templateImage.ubound(thirdDim)-d,this->templateImage.ubound(thirdDim));
		}
		else if ( maxZw == maxZ ){
			Kt = Range(0,d);
		}
		else{
			Kt = Range(0,this->templateImage.ubound(thirdDim));
		}

		Array< Pixel, 3 > I2( I1.shape() );
		I2 = this->templateImage( It, Jt, Kt ) * mask;

		I2 = I2 - sum(I2) / elements;
		Pixel norm2 = sqrt( 1.0 * sum( pow2( I2 ) ) );

		// compute and store CCC
		this->weightsOnFly( p ) = 10.0 * fabs( 1.0 - sum( I1 * I2 ) / norm1 / norm2 ) + 1e-5;
		if ( this->weightsOnFly( p ) <= 0 ){
			cout << this->weightsOnFly( p ) << " negative metric?" << endl;
			assert(0);
		}
	}
	return this->weightsOnFly( p );
}

template< class Pixel >
void nbfFastMarching3Dccc< Pixel > :: setTemplate( Array< Pixel, 3 > & T )
{
	this->templateImage.reference( T );

	this->maxX = this->weight.ubound(firstDim);
	this->maxY = this->weight.ubound(secondDim);
	this->maxZ = this->weight.ubound(thirdDim);

	this->templateX = floor( ( T.rows()  - 1 ) / 2.0 );
	this->templateY = floor( ( T.cols()  - 1 ) / 2.0 );
	this->templateZ = floor( ( T.depth() - 1 ) / 2.0 );
}