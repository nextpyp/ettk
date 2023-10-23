#ifndef FILE_nbfGaussianFilter
#define FILE_nbfGaussianFilter

// TODO
// DO AS LAPLACIAN EVOLUTION. MUCH EASIER AND FAST?
//

#include <nbfArrayFilter.h>

// YOU MUST SPECIFY THE SIZE IN COMPILE TIME

#define FLUJOS_GAUSS_SIZE 3
#define FLUJOS_GAUSS_OFFSET (FLUJOS_GAUSS_SIZE-1)/2

// work-around for fast 3D neighbor filter
static float G[FLUJOS_GAUSS_SIZE][FLUJOS_GAUSS_SIZE][FLUJOS_GAUSS_SIZE];

template< class Pixel, int const Dim >
class nbfGaussianFilter : public nbfArrayFilter< Pixel, Dim >
{
public:

	// constructor takes source array as input, and sigma value
	nbfGaussianFilter( Array< Pixel, 3 > &, Pixel = 1/3.0 );

	~nbfGaussianFilter(){};

	// 2D //
	// TODO

	// 3D //
	// using array iterators (slow)
	void execute( Array< Pixel, 3 > & );

	// using laplacian pde, specify number of iterations
	void executePde( Array< Pixel, 3 > &, int = 10 );
	
	// G = exp( - sigma * (x^2+y^2+z^2) );
	void setSigma( Pixel s ){ this->sigma = s; }
	
protected:

	// build 3d gaussian kernel
	void buildKernel3D();

	Pixel sigma;
	Array< Pixel, Dim > kernel;

};

template< class Pixel, int const Dim >
nbfGaussianFilter< Pixel, Dim > :: nbfGaussianFilter( Array< Pixel, 3 > & input, Pixel sigma )
: nbfArrayFilter< Pixel, Dim >( input )
{
	this->sigma = sigma;
}

template< class Pixel, int const Dim >
void nbfGaussianFilter< Pixel, Dim > :: buildKernel3D()
{
	this->kernel.resize( FLUJOS_GAUSS_SIZE, FLUJOS_GAUSS_SIZE, FLUJOS_GAUSS_SIZE );

	this->kernel = exp( - this->sigma * ( sqr( tensor::i - FLUJOS_GAUSS_OFFSET )
		+ sqr( tensor::j - FLUJOS_GAUSS_OFFSET ) 
		+ sqr( tensor::k - FLUJOS_GAUSS_OFFSET ) ) );

	// normalize
	this->kernel = this->kernel / sum( this->kernel );

	// transfer to static array
	for ( int i = 0; i < FLUJOS_GAUSS_SIZE; i++ ){
		for ( int j = 0; j < FLUJOS_GAUSS_SIZE; j++ ){
			for ( int k = 0; k < FLUJOS_GAUSS_SIZE; k++ ){
				G[i][j][k] = kernel(i,j,k);
			}
		}
	}
}

template< class Pixel, int const Dim >
void nbfGaussianFilter< Pixel, Dim > :: execute( Array< Pixel, 3 > & output )
{
	// build 3d kernel
	this->buildKernel3D();

	output.resize( this->input->shape() );

#if 0 // slow
	output.resize( this->input->shape() );

	Array< Pixel, Dim > :: iterator _output = output.begin();

	int offset = ( this->size - 1 ) / 2;
	int x,y,z;
	Range I(-offset,offset);

	while ( _output != output.end() ){

		x = _output.position()(firstDim);
		y = _output.position()(secondDim);
		z = _output.position()(thirdDim);
		Array< Pixel, 3 > B( (*this->input)( I + x, I + y, I + z ) );
		(*_output) = sum( this->kernel * B );

		++_output;
	}
#else // fast
	output = gauss3D( *this->input );
#endif
}

#define FLUJOS_DECLARE_GAUSS_OPERATOR(name,size,offset)				\
template<class T>													\
inline _bz_typename T::T_numtype									\
name(T& A )															\
{																	\
	T::T_numtype res = 0;											\
	for ( int i = 0; i < size; i++ ){								\
		for ( int j = 0; j < size; j++ ){							\
			for ( int k = 0; k < size; k++ ){						\
				res += G[i][j][k] * A(i-offset,j-offset,k-offset);	\
			}														\
		}															\
	}																\
	return res;														\
}																	\

FLUJOS_DECLARE_GAUSS_OPERATOR(gauss,FLUJOS_GAUSS_SIZE,FLUJOS_GAUSS_OFFSET)

// declare stencil operator
BZ_DECLARE_STENCIL_OPERATOR1(gauss3D,A)
return gauss(A);
BZ_END_STENCIL_OPERATOR

// allow direct call of stencil operator
BZ_ET_STENCIL(gauss3D,P_numtype)


template< class Pixel, int const Dim >
void nbfGaussianFilter< Pixel, Dim > :: executePde( Array< Pixel, 3 > & output, int iters )
{
	for ( int i = 0; i < iters; i++ ){
		output = Laplacian3D(*this->input);
		cout << "[" << min(output) << "," << max(output) << "]" << endl;
		(*this->input) = (*this->input) + .1 * output;
	}
	output = (*this->input);
}

#endif /* FILE_nbfGaussianFilter */
