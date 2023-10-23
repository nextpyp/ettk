#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "fftw3.h"

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

#include <vtkImageData.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkStructuredPointsWriter.h>

using namespace blitz;

#include <complex>

#include <nbfTimer.h>
#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>

int main( int argc, char *argv[] )
{
  nbfTimer t;

  Array< double, 3 > A( 100, 100, 100 );
  firstIndex i;
  secondIndex j;
  thirdIndex k;
  A = sqrt( pow2( i - A.rows() / 2.0 ) + pow2( j - A.cols() / 2.0 ) + pow2( k - A.depth() / 2.0 ) );

  //A = A * pow( -1, i + j + k );

  Array< complex< double >, 3 > F( A.shape() );
  Array< complex< double >, 3 > Fi( A.shape() );

  real(Fi) = A;
  imag(Fi) = 0;

  fftw_plan p;

  double * inD = reinterpret_cast<double*>(A.data());

  fftw_complex * in = reinterpret_cast<fftw_complex*>(Fi.data());
  fftw_complex * out = reinterpret_cast<fftw_complex*>(F.data());

  t.start();
  p = fftw_plan_dft_r2c_3d( A.rows(), A.cols(), A.depth(), inD, out, FFTW_PRESERVE_INPUT );    
  //p = fftw_plan_dft_r2c_3d( A.rows(), A.cols(), A.depth(), inD, out, FFTW_FORWARD, FFTW_ESTIMATE );    
  for ( int i = 0; i < 1; i++ ){
	fftw_execute(p); /* repeat as needed */
  }
  t.stop();
  //fftw_destroy_plan(p);
  cout << "fftw = " << t.elapsedSeconds() << endl;

  in = reinterpret_cast<fftw_complex*>(F.data());
  out = reinterpret_cast<fftw_complex*>(Fi.data());
    
  fftw_plan M_plan;
  t.start();
  M_plan = fftw_plan_dft_3d( A.rows(), A.cols(), A.depth(), in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
  fftw_execute(M_plan);
  //fftw_destroy_plan(M_plan);
  t.stop();
  cout << "ifftw = " << t.elapsedSeconds() << endl;

  vtkImageData * data = vtkImageData::New();
 // nbfVTKInterface::blitzToVtk( A, data );

 // vtkImageFFT * fft = vtkImageFFT::New();
 // fft->SetDimensionality(3);
 // fft->SetInput( data );
 // t.start();
 // for ( int i = 0; i < 1; i++ ){
	//fft->Modified();
	//fft->Update();
 // }
 // t.stop();
 // cout << " vtk fft = " << t.elapsedSeconds() << endl;

 // vtkImageRFFT * ifft = vtkImageRFFT::New();
 // ifft->SetDimensionality(3);
 // ifft->SetInput( fft->GetOutput() );
 // t.start();
 // ifft->Update();
 // t.stop();
 // cout << " vtk ifft = " << t.elapsedSeconds() << endl;

 // vtkImageFourierCenter * center = vtkImageFourierCenter::New();
 // center->SetInput( fft->GetOutput() );
 // center->Update();

 // Array< double, 3 > realP, imagP;
 // nbfVTKInterface::vtkToBlitz( ifft->GetOutput(), realP, imagP );

 // cout << "max fftw real = " << max( real(Fi) ) << endl;
 // cout << "max vtk real = " << max( realP ) << endl;

 // cout << "max fftw imag = " << max( imag(Fi) ) << endl;
 // cout << "max vtk imag = " << max( imagP ) << endl;

 // cout << "difference = " << sum( abs( real(Fi) - realP ) ) << endl;

 // cout << "difference = " << sum( abs( imag(Fi) - imagP ) ) << endl;

  vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
  //writer->SetFileName("vtkfft.vtk");
  //writer->SetInput( ifft->GetOutput() );
  //writer->Write();

  writer->SetFileName("fftw.vtk");
  Array< double, 3 > K( F.shape() );
  K = real(Fi) / K.size();

  //cout << K << endl;

  //K = K - min(K);
  //cout << min(K) << ", " << max(K) << endl;
  nbfVTKInterface::blitzToVtk( K, data );
  writer->SetInput( data );
  writer->Write();
  writer->Delete();

  //fft->Delete();
  data->Delete();

  return 0;
}

