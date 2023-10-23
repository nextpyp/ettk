#include "mpi.h"
// #define NBF_VERBOSE

#define BZ_GENERATE_GLOBAL_INSTANCES

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkImageReader.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageReslice.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkImageNonMaximumSuppression.h>
#include <vtkImageMathematics.h>
#include <vtkImageContinuousDilate3D.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkXMLPImageDataWriter.h>
#include <vtkPoints.h>
#include <vtkPolyVertex.h>
#include <vtkProbeFilter.h>
#include <vtkContourFilter.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMrcWriter.h>

#include <em/nbfFourierImageMetricCore.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfExtractPointsAndNormals3D.h>
#include <em/nbfTemplateSearchEM.h>
#include <em/nbfCorrelationImageMetric.h>

#include <em/nbfCutSubVolumes.h>


#define PIXEL double

int main( int argc, char ** argv )
{
	cout << " " << argv[0] << " compiled on " << __DATE__ << " at " << __TIME__ << ".\n\n";

	if ( argc != 3 ){
		std::cerr << "Usage: " << argv[0] << " input.mrc output.mrc " << endl;
		return 0;
	}

	char * input = argv[1];
	char * output = argv[2];

	nbfMrcReader reader;
	reader.setFileName( input );
	
	vtkImageData * data = vtkImageData :: New();
	reader.read( data );

	Array< PIXEL, 2 > A;
	nbfVTKInterface::vtkToBlitz( data, A );
	// A = A - mean(A);

	 nbfImageFilter< float, 2 > ifilter;
	 ifilter.setMaskSize( 900, 900, 900, 10 );
	 // ifilter.execute( A );
	
	nbfFourierFilter< PIXEL, 2 > ff;
	ff.initializeFFTW( TinyVector< int, 2 >( A.shape() ) );
	
	real( ff.blitzFFT1 ) = A * ff.shift;
	imag( ff.blitzFFT1 ) = 0;
	
	fftw_execute( ff.fftplan );
	
	Array< PIXEL, 2 > ctf( A.shape() );
	
	PIXEL box = 2048;
	// [x,y] = meshgrid(1:box,1:box);
	
	PIXEL pixel = 1.136;
	PIXEL voltage = 300 * 1e3; 
	PIXEL wgh = 0.07;
	PIXEL Cs = 2.7 * 1e7;
	PIXEL lambda = 12.26 / sqrt( voltage + 0.9785 * pow2(voltage) * 1e-6 );
	
	PIXEL z1 = 26522.04;
	PIXEL z2 = 20178.04;
	PIXEL angle = 58.90 * vtkMath::Pi() / 180;
	
	PIXEL centerX = A.cols()/2 + 1;
	PIXEL centerY = A.rows()/2 + 1;
		
	firstIndex i; secondIndex j;
	
	Array< PIXEL, 2 > theta( A.shape() );
	theta = atan2( i - centerX, j - centerY );
	
	Array< PIXEL, 2 > radius( A.shape() );
	radius = ( pow2( centerX - i ) + pow2( centerY - j ) ) / pow2(pixel) / pow2(box);
	
	Array< PIXEL, 2 > defocus( A.shape() );
	defocus = .5 * ( z1 + z2 + ( z1 - z2 ) * cos( 2 * ( theta - angle ) ) );
	
	Array< PIXEL, 2 > phase( A.shape() );
	phase = vtkMath::Pi() * lambda * defocus * radius - vtkMath::Pi()/2 * Cs * lambda * lambda * lambda * pow2(radius);
	
	ctf = - sqrt( 1 - pow2(wgh) ) * sin( phase ) - wgh * cos( phase );
	
	ff.blitzFFT2 *= ctf;
	
	//TinyVector<PIXEL,2> t(.4,.1);
	//ff.lowPassOn(t);
	//ff.execute( ff.blitzFFT2 );	
	
	fftw_execute(ff.ifftplan);	
	A = real(ff.blitzFFT1) * ff.shift / A.size();
	
	Array< float, 3 > B( A.rows(), A.cols(), 1 );
	// cout << blitz::minmax::min(A) << ", " << blitz::minmax::min(A) << endl;
	// A = A - blitz::minmax::min(A);
	// A = A / blitz::minmax::min(A);
	B( Range :: all(), Range :: all(), 0 ) = cast<float>( A );
	
	nbfMrcWriter mrcw;
	mrcw.setFileName( output );
	mrcw.write( B );

	B( Range :: all(), Range :: all(), 0 ) = cast<float>( ctf );
	mrcw.setFileName( "ctf.mrc" );
	mrcw.write( B );

	
	data->Delete();
}
