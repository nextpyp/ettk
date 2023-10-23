#ifndef FILE_nbfLinearInterpolator
#define FILE_nbfLinearInterpolator

// Class nbfLinearInterpolator.
//
// Implements 2D and 3D multilinear interpolation.
// 
// INFTY values are considered as being out of the domain. Therefore,
// when doing the interpolation if any of the neighbors is INFTY
// the interpolated value will be INFTY.

template< class Pixel, int const Dim >
class nbfLinearInterpolator
{
public:

	// constructor takes source values array as input
	nbfLinearInterpolator( Array< Pixel, Dim > & );

	~nbfLinearInterpolator(){};

	// bi-linear interpolation
	// arguments: x and y interpolation positions and result
	void interpolate( Array< Pixel, Dim > &, Array< Pixel, Dim > &, Array< Pixel, Dim > &, Array< bool, Dim > & );

	// tri-linear interpolation
	// arguments: x, y and z interpolation positions and result
	void interpolate( Array< Pixel, Dim > &, Array< Pixel, Dim > &, Array< Pixel, Dim > &, Array< Pixel, Dim > &, Array< bool, Dim > & );

	// bi-linear single point interpolation
	// return INFTY if outside domain
	Pixel interpolateSingle( Pixel &, Pixel &, bool & );
	Pixel interpolateSingle( TinyVector< Pixel, 2 > & );

	TinyVector< Pixel, 2 > interpolateSingleGradient( Pixel &, Pixel & );
	TinyVector< Pixel, 2 > interpolateSingleGradient( TinyVector< Pixel, 2 > & );

	// tri-linear single point interpolation
	// return INFTY if outside domain
	Pixel interpolateSingle( Pixel &, Pixel &, Pixel &, bool & );
	Pixel interpolateSingle( TinyVector< Pixel, 3 > & );

	Pixel interpolateSingleClosest( Pixel &, Pixel &, Pixel &, bool & ); 
	Pixel interpolateSingleClosest( TinyVector< Pixel, 3 > &, bool & );

private:

	// source image
	Array< Pixel, Dim > * input;
	//Array< Pixel, Dim > inputShape;
};

template< class Pixel, int const Dim >
nbfLinearInterpolator< Pixel, Dim > :: nbfLinearInterpolator( Array< Pixel, Dim > & A )
{
	this->input = &A;
	//this->inputShape.resize( A.shape() );
//	this->enlarged = (*this->input);
}

template< class Pixel, int const Dim >
Pixel nbfLinearInterpolator< Pixel, Dim > :: interpolateSingle( Pixel & x, Pixel & y, bool & inside )
{
	//inside = true;
	//int l = floor( x ); int c = ceil( x );
	//int k = floor( y ); int d = ceil( y );
	//Pixel a = x - l;
	//Pixel b = y - k;

	//// check if position valid in the original array
	//if ( this->input.isInRange(l,k) && this->input.isInRange(c,d) ){
	//	if ( ( l == c ) & ( k == d ) ){
	//		return ( (*this->input)(l,k) );
	//	}
	//	Array< Pixel, 2 > A( (*this->input)( Range(l,l+1), Range(k,k+1) ) );
	//	// if at least one neighbor is INFTY, then we are out of the domain
	//	if ( sum( where( A == numeric_limits<Pixel>::max(), 1, 0 ) ) ){
	//		return numeric_limits<Pixel>::max();
	//	}
	//	else{
	//		return ( (1-a) * (1-b) * A(0,0) + a * ( 1 - b ) * A(1,0) +
	//			(1-a) * b * A(0,1) + a * b * A(1,1) );
	//	}
	//}
	//else{
	//	inside = false;
	//}

	// new

	inside = true;

	// local interpolation coordinates
	Pixel a, b;

	// interpolate corners
	int xF, xC, yF, yC;

	// x - check bounds
	if ( x <= this->input->lbound(firstDim) ){
		xF = this->input->lbound(firstDim);
		a = 0;
		inside = false;
	}
	else if ( x >= this->input->ubound(firstDim) ){
		xF = this->input->ubound(firstDim) - 1;
		a = 1;
		inside = false;
	}
	else{
		xF = floor(x);
		a = x - xF;
	}
	xC = xF + 1;
	
	// y - check bounds
	if ( y <= this->input->lbound(secondDim) ){
		yF = this->input->lbound(secondDim);
		b = 0;
		inside = false;
	}
	else if ( y >= this->input->ubound(secondDim) ){
		yF = this->input->ubound(secondDim) - 1;
		b = 1;
		inside = false;
	}
	else{
		yF = floor(y);
		b = y - yF;
	}
	yC = yF + 1;
	
	// get image values at corners
	Array< Pixel, 2 > A( (*this->input)( Range(xF,xC), Range(yF,yC) ) );

	// interpolate
	return ( (1-a) * (1-b) * A(0,0) + a * ( 1 - b ) * A(1,0) +
		(1-a) * b * A(0,1) + a * b * A(1,1) );

}

template< class Pixel, int const Dim >
Pixel nbfLinearInterpolator< Pixel, Dim > :: interpolateSingle( TinyVector< Pixel, 2 > & x )
{
	bool dummy;
	return this->interpolateSingle( x(firstDim), x(secondDim), dummy );
}
	
template< class Pixel, int const Dim >
TinyVector< Pixel, 2 > nbfLinearInterpolator< Pixel, Dim > :: interpolateSingleGradient( Pixel & x, Pixel & y )
{
	TinyVector< Pixel, 2 > gradient;
	int l = floor( x ); int c = ceil( x );
	int k = floor( y ); int d = ceil( y );
	Pixel a = x - l;
	Pixel b = y - k;

	if ( l == this->input.ubound(firstDim) ){
		gradient[firstDim] = (*this->input)(l,k) - (*this->input)(l-1,k);
	}
	else{
		Array< Pixel, 2 > A( (*this->input)( Range(l,l+1), Range(k,k+1) ) );
		gradient[firstDim] = - (1-b) * A(0,0) + (1-b) * A(1,0) - b * A(0,1) + b * A(1,1);
	}
	if ( k == this->input.ubound(secondDim) ){
		gradient[secondDim] = (*this->input)(l,k) - (*this->input)(l,k-1);
		//cout << l << " < " << x << " < " << c << endl;
		//cout << k << " < " << y << " < " << d << endl;
		//cout << (*this->input)(l,k) << " - " << (*this->input)(l,k-1) << " = " << gradient[secondDim] << endl;
	}
	else{
		Array< Pixel, 2 > A( (*this->input)( Range(l,l+1), Range(k,k+1) ) );
		gradient[secondDim] = - (1-a) * A(0,0) - a * A(1,0) + (1-a) * A(0,1) + a * A(1,1);
	}

	return gradient;
}

template< class Pixel, int const Dim >
TinyVector< Pixel, 2 > nbfLinearInterpolator< Pixel, Dim > :: interpolateSingleGradient( TinyVector< Pixel, 2 > & x )
{
	return this->interpolateSingleGradient( x(firstDim), x(secondDim) );
}

template< class Pixel, int const Dim >
Pixel nbfLinearInterpolator< Pixel, Dim > :: interpolateSingle( Pixel & x, Pixel & y, Pixel & z, bool & inside  )
{
	inside = true;

	// local interpolation coordinates
	Pixel a, b, c;

	// interpolate corners
	int xF, xC, yF, yC, zF, zC;

	// x - check bounds
	if ( x <= this->input->lbound(firstDim) ){
		xF = this->input->lbound(firstDim);
		a = 0;
		inside = false;
	}
	else if ( x >= this->input->ubound(firstDim) ){
		xF = this->input->ubound(firstDim) - 1;
		a = 1;
		inside = false;
	}
	else{
		xF = floor(x);
		a = x - xF;
	}
	xC = xF + 1;
	
	// y - check bounds
	if ( y <= this->input->lbound(secondDim) ){
		yF = this->input->lbound(secondDim);
		b = 0;
		inside = false;
	}
	else if ( y >= this->input->ubound(secondDim) ){
		yF = this->input->ubound(secondDim) - 1;
		b = 1;
		inside = false;
	}
	else{
		yF = floor(y);
		b = y - yF;
	}
	yC = yF + 1;
	
	// z - check bounds
	if ( z <= this->input->lbound(thirdDim) ){
		zF = this->input->lbound(thirdDim);
		c = 0;
		inside = false;
	}
	else if ( z >= this->input->ubound(thirdDim) ){
		zF = this->input->ubound(thirdDim) - 1;
		c = 1;
		inside = false;
	}
	else{
		zF = floor(z);
		c = z - zF;
	}
	zC = zF + 1;
	
	// get image values at corners
	Array< Pixel, 3 > A( (*this->input)( Range(xF,xC), Range(yF,yC), Range(zF,zC) ) );

	// interpolate
	return ( ( 1 - a ) * ( 1 - b ) * ( 1 - c ) * A(0,0,0) + 
		a * ( 1 - b ) * ( 1 - c ) * A(1,0,0) +
		( 1 - a ) * b * ( 1 - c ) * A(0,1,0) + 
		( 1 - a ) * ( 1 - b ) * c * A(0,0,1) +
		a * ( 1 - b ) * c * A(1,0,1) +
		( 1 - a ) * b * c * A(0,1,1) +
		a * b * ( 1 - c ) * A(1,1,0) +
		a * b * c * A(1,1,1) );

	//inside = this->input->isInRange(x,y,z);

	//int xF = floor( x ); int xC = ceil( x );
	//int yF = floor( y ); int yC = ceil( y );
	//int zF = floor( z ); int zC = ceil( z );

	//Pixel a = x - xF;
	//Pixel b = y - yF;
	//Pixel c = z - zF;

	//Range I,J,K;

	//// check if position valid in the original array
	//if ( this->input->isInRange(xF,yF,zF) && this->input->isInRange(xC,yC,zC) ){
	//	if ( xF < this->input->ubound(firstDim) ){
	//		I = Range(xF,xF+1);
	//	}
	//	else{
	//		I = Range(xF-1,xF);
	//	}
	//	if ( yF < this->input->ubound(secondDim) ){
	//		J = Range(yF,yF+1);
	//	}
	//	else{
	//		J = Range(yF-1,yF);
	//	}
	//	if ( zF < this->input->ubound(thirdDim) ){
	//		K = Range(zF,zF+1);
	//	}
	//	else{
	//		K = Range(zF-1,zF);
	//	}
	//	//if ( ( xF == xC ) & ( yF == yC ) & ( zF == zC ) ){
	//	//	return ( (*this->input)(xF,yF,zF) );
	//	//}
	//	//Array< Pixel, 3 > A( (*this->input)( Range(xF,xF+1), Range(yF,yF+1), Range(zF,zF+1) ) );
	//	Array< Pixel, 3 > A( (*this->input)( I, J, K ) );
	//	// if at least one neighbor is INFTY, then we are out of the domain
	//	if ( sum( where( A == numeric_limits<Pixel>::max(), 1, 0 ) ) ){
	//		return numeric_limits<Pixel>::max();
	//	}
	//	else{
	//		return ( ( 1 - a ) * ( 1 - b ) * ( 1 - c ) * A(0,0,0) + 
	//			a * ( 1 - b ) * ( 1 - c ) * A(1,0,0) +
	//			( 1 - a ) * b * ( 1 - c ) * A(0,1,0) + 
	//			( 1 - a ) * ( 1 - b ) * c * A(0,0,1) +
	//			a * ( 1 - b ) * c * A(1,0,1) +
	//			( 1 - a ) * b * c * A(0,1,1) +
	//			a * b * ( 1 - c ) * A(1,1,0) +
	//			a * b * c * A(1,1,1) );
	//	}
	//}
	//else{
	//	return numeric_limits<Pixel>::max();
	//}
}

template< class Pixel, int const Dim >
Pixel nbfLinearInterpolator< Pixel, Dim > :: interpolateSingle( TinyVector< Pixel, 3 > & x )
{
	bool dummy;
	return this->interpolateSingle( x(firstDim), x(secondDim), x(thirdDim), dummy );
}

template< class Pixel, int const Dim >
Pixel nbfLinearInterpolator< Pixel, Dim > :: interpolateSingleClosest( Pixel & x, Pixel & y, Pixel & z, bool & inside )
{
	int xF = floor( x ); int xC = ceil( x );
	int yF = floor( y ); int yC = ceil( y );
	int zF = floor( z ); int zC = ceil( z );

	if ( x - xF < .5 ){
		xC = xF;
	}

	if ( y - yF < .5 ){
		yC = yF;
	}

	if ( z - zF < .5 ){
		zC = zF;
	}

	inside = this->input->isInRange(xC,yC,zC);

	xC = extrema::max( xC, this->input->lbound(0) );
	xC = extrema::min( xC, this->input->ubound(0) );
	yC = extrema::max( yC, this->input->lbound(1) );
	yC = extrema::min( yC, this->input->ubound(1) );
	zC = extrema::max( zC, this->input->lbound(2) );
	zC = extrema::min( zC, this->input->ubound(2) );

	return (*this->input)( xC, yC, zC );
}

template< class Pixel, int const Dim >
Pixel nbfLinearInterpolator< Pixel, Dim > :: interpolateSingleClosest( TinyVector< Pixel, 3 > & x, bool & inside )
{
	return this->interpolateSingleClosest( x(firstDim), x(secondDim), x(thirdDim), inside );
}

template< class Pixel, int const Dim >
void nbfLinearInterpolator< Pixel, Dim > :: interpolate( Array< Pixel, Dim > & sourceX,
											   Array< Pixel, Dim > & sourceY,
											   Array< Pixel, Dim > & interp,
											   Array< bool, Dim > & inside )
{
	// Add one row of zeros in each dimension for adequate edge handling

	// this keeps the original array in place, only incrementing the size in 1
	this->input->resizeAndPreserve( this->input->shape() + 1 );

	// force to '0' newly created extension
	(*this->input)( this->input->ubound(firstDim), Range::all() ) = 0;
	(*this->input)( Range::all(), this->input->ubound(secondDim) ) = 0;

	// set output size
	interp.resize( sourceX.shape() );
	inside.resize( sourceX.shape() );

	// interpolate for each array position
	typename Array< Pixel, Dim > :: iterator iSourceX = sourceX.begin(),
		iSourceY = sourceY.begin();
	typename Array< bool, Dim > :: iterator iInside = inside.begin();
	typename Array< Pixel, Dim > :: iterator iInterp = interp.begin();

	while ( iSourceX != sourceX.end() ){
		(*iInterp) = this->interpolateSingle( *iSourceX, *iSourceY, *iInside );
		++iSourceX;
		++iSourceY;
		++iInterp;
		++iInside;
	}

	// remove the extra row of zeros in each dimension
	this->input->resizeAndPreserve( this->input->shape() );
}

template< class Pixel, int const Dim >
void nbfLinearInterpolator< Pixel, Dim > :: interpolate( Array< Pixel, Dim > & sourceX,
											   Array< Pixel, Dim > & sourceY,
											   Array< Pixel, Dim > & sourceZ,
											   Array< Pixel, Dim > & interp,
											   Array< bool, Dim > & inside )
{
	//// Add one row of zeros in each dimension for adequate edge handling

	//// this keeps the original array in place, only incrementing the size in 1
	//this->input->resizeAndPreserve( this->input.shape() + 1 );

	//// force to '0' newly created extension
	//(*this->input)( this->input->ubound(firstDim), Range::all(), Range::all() ) = 0;
	//(*this->input)( Range::all(), this->input->ubound(secondDim), Range::all() ) = 0;
	//(*this->input)( Range::all(), Range::all() , this->input->ubound(thirdDim) ) = 0;

	// set output size same as input
	interp.resize( sourceX.shape() );
	inside.resize( sourceX.shape() );

	// interpolate for each array position
	typename Array< Pixel, Dim > :: iterator iSourceX = sourceX.begin(),
		iSourceY = sourceY.begin(),
		iSourceZ = sourceZ.begin();
	typename Array< bool, Dim > :: iterator iInside = inside.begin();
	typename Array< Pixel, Dim > :: iterator iInterp = interp.begin();

	while ( iSourceX != sourceX.end() ){
		(*iInterp) = this->interpolateSingle( *iSourceX, *iSourceY, *iSourceZ, *iInside );
		++iSourceX;
		++iSourceY;
		++iSourceZ;
		++iInterp;
		++iInside;
	}

	//// remove the extra row of zeros in each dimension
	//this->input->resizeAndPreserve( this->input.shape() );
}

#endif /* FILE_nbfLinearInterpolator */