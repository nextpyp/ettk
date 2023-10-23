#define NBF_AVERAGE_IN_RECIPROCAL_SPACE

#ifdef WIN32
	 #define NBF_VERBOSE
#endif

#include "mpi.h"

#define BZ_GENERATE_GLOBAL_INSTANCES

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

#include <blitz/bzdebug.h>        // Test suite globals

using namespace blitz;

#include <string.h>

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
#include <vtkImageCast.h>
#include <vtkImageExtractComponents.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMrcWriter.h>
#include <bs/nbfBordStrategyMirror.h>

#include <em/nbfImageMetric.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfWedgedAverageImage3D.h>

#include <em/nbfUnsupervisedLoopClustering.h>


#include <nbfCylindricalDomain3.h>

#include <mxml.h>

//extern "C" {
//#include <../svdlibc/svdlib.h>
//}

#include <nbfTimer.h>

#define PIXEL float

int main( int argc, char ** argv )
{	
	nbfMrcReader reader;
	reader.setFileName( argv[1] );
	vtkImageData * data = vtkImageData :: New();
	reader.read( data );

	Array< float, 3 > A;
	nbfVTKInterface :: vtkToBlitzReference( data, A );

	Array< float, 3 > cylinder( A.shape() );
	firstIndex i;
	secondIndex j;
	cylinder = atan2( i - 50.0, j - 50.0 ) * vtkMath::RadiansToDegrees();

	Array< float, 3 > cylinderradius( A.shape() );
	cylinderradius = sqrt( pow2(i - A.rows()/2.0) + pow2(j - A.cols()/2.0 ) );
	cylinderradius = where( ( cylinderradius < 40 ), 1, 0 );
	
	// Erin's VRC03
	cylinder = where( ( ( cylinder > -80 ) && ( cylinder < 10 ) ), 1, 0 );
	cylinder( Range :: all(), Range :: all(), Range(fromStart,36) ) = 0;
	cylinder( Range :: all(), Range :: all(), Range(80,toEnd) ) = 0;

	cylinder *= cylinderradius;

	// SIVmne_sCD4
	//cylinder = where( ( ( cylinder > -60 ) && ( cylinder < 60 ) ), 1, 0 );
	//cylinder( Range :: all(), Range :: all(), Range(0,35) ) = 0;

	//// Bal
	//cylinder = where( ( ( cylinder > -60 ) && ( cylinder < 60 ) ), 1, 0 );
	//cylinder( Range :: all(), Range :: all(), Range(0,35) ) = 0;

	// Fab
	//cylinder = where( ( ( cylinder > -20 ) && ( cylinder < 100 ) ), 1, 0 );
	//cylinder( Range :: all(), Range :: all(), Range(0,41) ) = 0;
	
	// CD4
	//	cylinder = where( ( ( cylinder > -35 ) && ( cylinder < 85 ) ), 1, 0 );
	//  cylinder( Range :: all(), Range :: all(), Range(0,41) ) = 0;

	vtkImageGaussianSmooth * smooth = vtkImageGaussianSmooth :: New();
	vtkImageData * cropVtk = vtkImageData :: New();
	nbfVTKInterface::blitzToVtk( cylinder, cropVtk );
	smooth->SetInput( cropVtk );
	smooth->SetDimensionality(3);
	smooth->SetStandardDeviations(1.0,1.0);
	smooth->SetRadiusFactors(1.0,1.0);
	smooth->Update();
	nbfVTKInterface::vtkToBlitz( smooth->GetOutput(), cylinder );
	smooth->Delete();

	A = cylinder;
	// A *= cylinder;
	// A *= -1;

	vtkTransform * t = vtkTransform :: New();
	
	// Erin's VRC03
	t->Translate(0,12,5);
		
	//// Bal
	//t->Translate(0,11,1);

	// Fab
	// t->Translate(9,7,1);

	// CD4
	// t->Translate(5,10,4);

	//vtkImageReslice * reslice = vtkImageReslice :: New();
	//reslice->SetInput( data );
	//reslice->SetResliceTransform(t);
	//reslice->SetBackgroundLevel(0.0);
	//reslice->SetInterpolationModeToCubic();
	//reslice->Update();

	nbfMrcWriter writer1;
	writer1.setFileName( argv[2] );
	writer1.write( data );
	//writer1.write( reslice->GetOutput() );
}