#ifndef FILE_nbfVeselnessFilter
#define FILE_nbfVeselnessFilter

#include <vtkMath.h>
#include <nbfArrayFilter.h>
#include <nbfDifferentials.h>
#include <nbfGaussianFilter.h>
#include <nbfLinearInterpolator.h>

template< class Pixel, int const Dim >
class nbfVeselnessFilter : public nbfArrayFilter< Pixel, Dim >
{
public:

	// constructor takes weight array as input
	nbfVeselnessFilter( Array< Pixel, Dim > & );

	~nbfVeselnessFilter(){};

	// 3D //
	void execute( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< Pixel, 3 > & );
	void execute( Array< Pixel, 3 > &, int );
	void execute( Array< Pixel, 3 > & );

	void nonMaxSupress( Array< Pixel, 3 > &, Array< Pixel, 3 > & );
	void nonNewMaxSupress( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< Pixel, 3 > & );
};

template< class Pixel, int const Dim >
nbfVeselnessFilter< Pixel, Dim > :: nbfVeselnessFilter( Array< Pixel, Dim > & input )
: nbfArrayFilter< Pixel, Dim >( input )
{
}

template< class Pixel, int const Dim >
void nbfVeselnessFilter< Pixel, Dim > :: execute( Array< Pixel, 3 > & T )
{
	T.resize( this->input->shape() );

	Array< Pixel, 3 > Tx( T.shape() );
	Array< Pixel, 3 > Ty( T.shape() );
	Array< Pixel, 3 > Tz( T.shape() );

	// compute measure and orientation of features
	applyStencil(veselness(), (*this->input), T, Tx, Ty, Tz );

	// non-maxima suppresion
	// applyStencil( nonMaxSupress3D(), T, Tx );
}

template< class Pixel, int const Dim >
void nbfVeselnessFilter< Pixel, Dim > :: execute( Array< Pixel, 3 > & T, int dim )
{
	// initialize intermediate array	
	Array< Pixel, 3 > output( this->input->shape() );

	Array< Pixel, 3 > Tx( output.shape() );
	Array< Pixel, 3 > Ty( output.shape() );
	Array< Pixel, 3 > Tz( output.shape() );

	T.resize( output.shape() );

	// compute measure and orientation of features
	applyStencil(veselness(), (*this->input), output, Tx, Ty, Tz );

	//output( output.lbound(firstDim), Range::all(), Range::all() ) = output( output.lbound(firstDim) + 1, Range::all(), Range::all() );
	//output( output.ubound(firstDim), Range::all(), Range::all() ) = output( output.ubound(firstDim) - 1, Range::all(), Range::all() );
	//output( Range::all(), output.lbound(secondDim), Range::all() ) = output( Range::all(), output.lbound(secondDim) + 1, Range::all() );
	//output( Range::all(), output.ubound(secondDim), Range::all() ) = output( Range::all(), output.ubound(secondDim) - 1, Range::all() );
	//output( Range::all(), Range::all(), output.lbound(thirdDim) ) = output( Range::all(), Range::all(), output.lbound(thirdDim) + 1 );
	//output( Range::all(), Range::all(), output.ubound(thirdDim) ) = output( Range::all(), Range::all(), output.ubound(thirdDim) - 1 );

	//cout << min(output) << ", " << max(output) << endl;

	T = output;
	T = T / max(T);

	//return;

	switch (dim){

	case firstDim:

		// non-max suppression in YZ plane
		//Tx = where( output != 0, atan(Tz/Ty), 0 );
		Tx = where( output != 0, Tx, 0 );

		for ( int i = 0; i < output.rows(); i++ ){
			Array< Pixel, 2 > In( output( i, Range::all(), Range::all() ) );
			Array< Pixel, 2 > Out( Tx( i, Range::all(), Range::all() ) );
			applyStencil( nonMaxSupress2D(), In, Out );
		}
		Tx = Tx / max(Tx);
		T = Tx;
		break;

	case secondDim:

		// non-max suppression in XZ plane
		//Tz = where( output != 0, atan(Tz/Tx), 0 );
		Ty = where( output != 0, Ty, 0 );

		for ( int j = 0; j < output.cols(); j++ ){
			Array< Pixel, 2 > In( output( Range::all(), j, Range::all() ) );
			Array< Pixel, 2 > Out( Ty( Range::all(), j, Range::all() ) );
			applyStencil( nonMaxSupress2D(), In, Out );
		}
		Ty = Ty / max(Ty);
		T = Ty;
		break;

	case thirdDim:

		// non-max suppression in XY plane
		Tz = where( output != 0, Tz, 0 );

		for ( int k = 0; k < output.depth(); k++ ){
			Array< Pixel, 2 > In( output( Range::all(), Range::all(), k ) );
			Array< Pixel, 2 > Out( Tz( Range::all(), Range::all(), k ) );
			applyStencil( nonMaxSupress2D(), In, Out );
		}
		Tz = Tz / max(Tz);
		T = Tz;
		break;
	}
}

template< class Pixel, int const Dim >
void nbfVeselnessFilter< Pixel, Dim > :: execute( Array< Pixel, 3 > & Tx, 
												  Array< Pixel, 3 > & Ty,
												  Array< Pixel, 3 > & Tz )
{
	// initialize intermediate array	
	Array< Pixel, 3 > output( this->input->shape() );

	Tx.resize( output.shape() );
	Ty.resize( output.shape() );
	Tz.resize( output.shape() );

	// compute measure and orientation of features
	applyStencil(veselness(), (*this->input), output, Tx, Ty, Tz );

	// perform non-maxima supression using the orientation above
	// get result in theta to optimize memory
	//this->nonMaxSupress(output,theta);

	// non-max suppression in XY plane
	//Ty = where( output != 0, atan(Ty/Tx), 0 );
	Tz = where( output != 0, Tz, 0 );
	
	for ( int k = 0; k < output.depth(); k++ ){
		Array< Pixel, 2 > In( output( Range::all(), Range::all(), k ) );
		Array< Pixel, 2 > Out( Tz( Range::all(), Range::all(), k ) );
		applyStencil( nonMaxSupress2D(), In, Out );
	}

	// non-max suppression in XZ plane
	//Tz = where( output != 0, atan(Tz/Tx), 0 );
	Ty = where( output != 0, Ty, 0 );
	
	for ( int j = 0; j < output.cols(); j++ ){
		Array< Pixel, 2 > In( output( Range::all(), j, Range::all() ) );
		Array< Pixel, 2 > Out( Ty( Range::all(), j, Range::all() ) );
		applyStencil( nonMaxSupress2D(), In, Out );
	}

	// non-max suppression in YZ plane
	//Tx = where( output != 0, atan(Tz/Ty), 0 );
	Tx = where( output != 0, Tx, 0 );
	
	for ( int i = 0; i < output.rows(); i++ ){
		Array< Pixel, 2 > In( output( i, Range::all(), Range::all() ) );
		Array< Pixel, 2 > Out( Tx( i, Range::all(), Range::all() ) );
		applyStencil( nonMaxSupress2D(), In, Out );
	}

	//// combine both into single image (if local maxima in both planes)
	//output = where( ( Ty > 0 ) & ( Tz > 0 ), 1, 0 );
	//output = Ty;
	Tx = Tx / max(Tx);
	Ty = Ty / max(Ty);
	Tz = Tz / max(Tz);
}


template< class Pixel, int const Dim >
void nbfVeselnessFilter< Pixel, Dim > :: nonMaxSupress( Array< Pixel, 3 > & M, Array< Pixel, 3 > & T )
{
	// non-maxima supression slice-wise

	for ( int k = 0; k < M.depth(); k++ ){
		Array< Pixel, 2 > In( M( Range::all(), Range::all(), k ) );
		Array< Pixel, 2 > Out( T( Range::all(), Range::all(), k ) );
		applyStencil( nonMaxSupress2D(), In, Out );
	}
}

template< class Pixel, int const Dim >
void nbfVeselnessFilter< Pixel, Dim > :: nonNewMaxSupress( Array< Pixel, 3 > & M, 
														  Array< Pixel, 3 > & output, 
														  Array< Pixel, 3 > & Tx, 
														  Array< Pixel, 3 > & Ty, 
														  Array< Pixel, 3 > & Tz )
{
	// non-maxima supression slice-wise
	output = 0;
	Pixel lambda, dForward, dBackward;
	nbfLinearInterpolator< Pixel, 3 > interpolator(M);
	TinyVector< Pixel, 3 > position;
	Array< Pixel, 3 > :: iterator iM = M.begin();
	Array< Pixel, 3 > :: iterator iT = output.begin();
	Array< Pixel, 3 > :: iterator iTx = Tx.begin();
	Array< Pixel, 3 > :: iterator iTy = Ty.begin();
	Array< Pixel, 3 > :: iterator iTz = Tz.begin();
	int i = 0;
	int x,y,z;
	while ( iM != M.end() ){
		if ( fabs(*iTx) + fabs(*iTy) + fabs(*iTz) > 0 ){ 
			// cout << (*iTx) << ", " << (*iTy) << ", " << (*iTz) << endl;
			lambda = 1 / minmax::max( minmax::max( fabs(*iTx), fabs(*iTy) ), fabs(*iTz) );
			//cout << iM.position() << endl;
			position = iM.position() + lambda * TinyVector< Pixel, 3 >(*iTx,*iTy,*iTz);
			//cout << position << endl;
			dForward = interpolator.interpolateSingleClosest(position);
			position = iM.position() - lambda * TinyVector< Pixel, 3 >(*iTx,*iTy,*iTz);
			//cout << position << endl;
			dBackward = interpolator.interpolateSingleClosest(position);
			// if local maxima
			if ( ( (*iM) > dForward ) && ( (*iM) > dBackward ) ){
				(*iT) = (*iM);
			}
		}
        ++iM; ++iT; ++iTx; ++iTy; ++iTz; i++;
	}
}

BZ_DECLARE_STENCIL5(veselness,I,M,Tx,Ty,Tz)
	T1::T_numtype A[3][3], w[3], V[3][3];
	T1::T_numtype *ATemp[3],*VTemp[3];

	A[0][0] = central22n(I,firstDim);				// Axx
	A[0][1] = mixed22n(I,firstDim,secondDim);		// Axy
	A[0][2] = mixed22n(I,firstDim,thirdDim);		// Axz
	A[1][0] = A[0][1];								// Axy
	A[1][1] = central22n(I,secondDim);				// Ayy
	A[1][2] = mixed22n(I,secondDim,thirdDim);		// Ayz
	A[2][0] = A[0][2];								// Axz
	A[2][1] = A[1][2];								// Ayz
	A[2][2] = central22n(I,thirdDim);				// Azz

	for ( int i = 0; i < 3; i++ ){
		ATemp[i] = A[i];
		VTemp[i] = V[i];
	}

	// diagonalize - eigenvalues 'w' in *descending* order (signed)
	vtkMath::Jacobi(ATemp,w,VTemp);

	// sort eigenvalues according to argument: |lambda[0]| > |lambda[1]| > |lambda[2]|
	T1::T_numtype lambda[3];

	if ( fabs( w[0] ) > fabs( w[2] ) )
	{
		// eigenvector (don't need to compute it fully, only vPy / vPx)
		// vPx = V[0][0], vPy = V[1][0], vPz = V[2][0];

		//T = atan( V[1][0] / V[0][0] );
		//Tx = V[0][0]; 
		//Ty = V[1][0]; 
		//Tz = V[2][0];

		Tz = atan( V[1][0] / V[0][0] ); 
		Ty = atan( V[2][0] / V[0][0] ); 
		Tx = atan( V[2][0] / V[1][0] );

		lambda[0] = w[0];
		if ( fabs( w[1] ) > fabs( w[2] ) ){
			lambda[1] = w[1];
			lambda[2] = w[2];
		}
		else{
			lambda[1] = w[2];
			lambda[2] = w[1];
		}
	}
	else{

		// eigenvector (don't need to compute it fully, only vPy / vPx)
		// vPx = V[0][2], vPy = V[1][2], vPz = V[2][2];

		//T = atan( V[1][2] / V[0][2] );
		//Tx = V[0][2]; 
		//Ty = V[1][2]; 
		//Tz = V[2][2];

		Tz = atan( V[1][2] / V[0][2] ); 
		Ty = atan( V[2][2] / V[0][2] ); 
		Tx = atan( V[2][2] / V[1][2] );

		lambda[0] = w[2];
		if ( fabs( w[0] ) > fabs( w[1] ) ){
			lambda[1] = w[0];
			lambda[2] = w[1];
		}
		else{
			lambda[1] = w[1];
			lambda[2] = w[0];
		}
	}

	// veselness measure
	//M = vtkMath::Norm(w) * ( lambda[0] > 0 ) * lambda[0];
	
	// ellipsoid orientation
	M = vtkMath::Norm(w) * ( lambda[0] > 0 ) * ( lambda[1] > 0 );

	// lambda[1] and lambda[2] don't really influence, so this is almost the same as a above:
	// M = vtkMath::Norm(w) * ( lambda[0] > 0 ) * lambda[0] * ( 1 - fabs(lambda[1]) ) * ( 1 - fabs(lambda[2]) );
BZ_STENCIL_END

BZ_DECLARE_STENCIL2(nonMaxSupress2D,M,T)
	// E
	if ( fabs(T) <= vtkMath::Pi() / 6.0 )
	{
		if ( !( ( M(0,0) > M(-1,0) ) & ( M > M(1,0) ) ) )
		{
			T = 0;
		}
		else{
			T = M(0,0);
		}
	} 
	else 
	{
		// NE 
		if ( ( T > vtkMath::Pi()/6.0 ) & ( T <= vtkMath::Pi()/3.0 ) ){
			if ( !( ( M(0,0) > M(-1,-1) ) & ( M(0,0) > M(1,1) ) ) ){
				T = 0;
			}
			else{
				T = M(0,0);
			}
		} 
		else 
		{
			// N
			if ( fabs(T) > vtkMath::Pi()/3.0 ){
				if ( !( ( M(0,0) > M(0,-1) ) & ( M(0,0) > M(0,1) ) ) ){
					T = 0;
				}
				else{
					T = M(0,0);
				}
			} 
			else {
				// NW
				if ( !( ( M(0,0) > M(-1,1) ) & ( M(0,0) > M(1,-1) ) ) ){
					T = 0;
				}
				else{
					T = M(0,0);
				}
			}
		}
	}
BZ_STENCIL_END

BZ_DECLARE_STENCIL2(nonMaxSupress3D,M,T)
	// E
	if ( fabs(T) <= vtkMath::Pi() / 6.0 )
	{
		if ( !( ( M(0,0) > M(-1,0) ) & ( M > M(1,0) ) ) )
		{
			T = 0;
		}
		else{
			T = M(0,0);
		}
	} 
	else 
	{
		// NE 
		if ( ( T > vtkMath::Pi()/6.0 ) & ( T <= vtkMath::Pi()/3.0 ) ){
			if ( !( ( M(0,0) > M(-1,-1) ) & ( M(0,0) > M(1,1) ) ) ){
				T = 0;
			}
			else{
				T = M(0,0);
			}
		} 
		else 
		{
			// N
			if ( fabs(T) > vtkMath::Pi()/3.0 ){
				if ( !( ( M(0,0) > M(0,-1) ) & ( M(0,0) > M(0,1) ) ) ){
					T = 0;
				}
				else{
					T = M(0,0);
				}
			} 
			else {
				// NW
				if ( !( ( M(0,0) > M(-1,1) ) & ( M(0,0) > M(1,-1) ) ) ){
					T = 0;
				}
				else{
					T = M(0,0);
				}
			}
		}
	}
BZ_STENCIL_END

#endif /* FILE_nbfVeselnessFilter */