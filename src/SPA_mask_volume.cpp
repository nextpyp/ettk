#include "mpi.h"
#define NBF_VERBOSE

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


#define PIXEL float

int main( int argc, char ** argv )
{
	if ( argc != 8 ){
		cout << "USAGE: " << argv[0] << " [input mrc] [centerX centerY centerZ] [radius] [mode] [output mrc] \n" << endl;
		cout << "\t[mode] = 1: keep the spike only." << endl;
		cout << "\t       = 2: keep the background only." << endl;
		return 0;
	}

	int center_x = atoi( argv[2] );
	int center_y = atoi( argv[3] );
	int center_z = atoi( argv[4] );
	int radius = atoi( argv[5] );
	int mode = atoi( argv[6] );
	
	vtkImageData * data = vtkImageData :: New();
	nbfMrcReader reader;
	reader.setFileName( argv[1] );
	reader.read( data );
	
	Array< PIXEL, 3 > A;
	nbfVTKInterface::vtkToBlitzReference( data, A );

	Array< int, 3 > M( A.shape() );
	Range I( max( 0, center_x - radius / 2 ), min( center_x + radius / 2 - 1, M.ubound(firstDim) ) );
	Range J( max( 0, center_y - radius / 2 ), min( center_y + radius / 2 - 1, M.ubound(secondDim) ) );
	Range K( max( 0, center_z - radius / 2 ), min( center_z + radius / 2 - 1, M.ubound(thirdDim) ) );
	
	M = 0;
	M( I, J, K ) = 1;

	double m = mean( A( I, J, K ) );

	// firstIndex i;
	// secondIndex j;
	// thirdIndex k;
	// M = pow2( i - center_x ) + pow2( j - center_y ) + pow2( k - center_z );
	// M = where( M < radius * radius, 1, 0 );
	
	if ( mode == 1 ){	
		A = A * ( 1 - M ) + m * M;
	} else {
		A = A * M + m * ( 1 - M );
	}
	
	nbfMrcWriter writer;
	writer.setFileName( argv[7] );
	writer.write(data);
		
	data->Delete();
	
	return 0;
}
