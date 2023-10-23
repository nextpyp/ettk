#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

//#include "fftw3.h"
#define BZ_GENERATE_GLOBAL_INSTANCES

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

#include <blitz/bzdebug.h>        // Test suite globals

#include <vtkMath.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkStructuredPointsWriter.h>

using namespace blitz;

//#include <complex>

#include <nbfTimer.h>
#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>

int main( int argc, char *argv[] )
{
  nbfTimer t;

  firstIndex i;
  secondIndex j;
  thirdIndex k;
  fourthIndex l;

  Array< complex< float >, 3 > A( 100, 100, 100 );
  real( A ) = sqrt( pow2( i - A.rows() / 2.0 ) + pow2( j - A.cols() / 2.0 ) + pow2( k - A.depth() / 2.0 ) );
  imag( A ) = 0;

  nbfMatlabWriter w;
  w.setFileName("p.matlab");
  w.write( real(A) );

  int m = 10;

  Array< complex< float >, 2 > E( A.depth(), m );
  real(E) = cos( - 2 * vtkMath::Pi() * ( i - A.rows() / 2 ) * j / A.rows() );
  imag(E) = sin( - 2 * vtkMath::Pi() * ( i - A.rows() / 2 ) * j / A.rows() );
  
  w.write( real(E) );
  w.write( imag(E) );
  
  Array< complex< float >, 3 > Fw( A.rows(), A.cols(), m );
  Fw = sum( A(i,j,l) * E(l,k), l );
	
  Array< complex< float >, 3 > Fv( A.rows(), m, m );
  Fv = sum( Fw(i,l,j) * E(l,k), l );
	
  Array< complex< float >, 3 > F( m, m, m );
  F = sum( Fv(l,i,j) * E(l,k), l );

  w.write( real(F) );
  w.write( imag(F) );

  return 0;
}

