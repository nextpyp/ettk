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
	A = A - mean(A);

	PIXEL mask = ( A.rows() / 2 ) * .95;
	PIXEL apod = ( A.rows() / 2 ) * .05;
	nbfImageFilter< PIXEL, 2 > ifilter;
	ifilter.setMaskSize( mask, mask, mask, apod );
	ifilter.execute( A );
	
	Array< complex< PIXEL >, 2 > FFTA( A.shape() );

	nbfFourierFilter< PIXEL, 2 > ff;
	ff.initializeFFTW( TinyVector< int, 2 >( A.shape() ) );
	real( ff.blitzFFT1 ) = A * ff.shift;
	imag( ff.blitzFFT1 ) = 0;

	fftw_execute( ff.fftplan );
	
	TinyVector<PIXEL,2> t(.1,.1);
	ff.highPassOn(t);
	ff.execute( ff.blitzFFT2 );

	A = sqrt( real( ff.blitzFFT2 * conj( ff.blitzFFT2 ) ) );

	Array< float, 3 > B( A.rows(), A.cols(), 1 );
	A = A - min(A);
	A = A / max(A);
	A = log( 1 + 10 * A );
	B( Range :: all(), Range :: all(), 0 ) = cast<float>( A );
	nbfMrcWriter mrcw;
	mrcw.setFileName( output );
	mrcw.write( B );

	data->Delete();
}
