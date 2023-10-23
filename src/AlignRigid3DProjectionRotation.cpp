#define NBF_VERBOSE 1

#include "mpi.h"

#define BZ_GENERATE_GLOBAL_INSTANCES

#define NBF_VERBOSE

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

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
#include <em/nbfFourierImageMetric.h>
#include <em/nbfLoopClustering.h>
#include <em/nbfWedgedSubImage3D.h>

#define PIXEL float

int main( int argc, char ** argv )
{
	cout << argv[0] << "\n";
	cout << " Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";

	if ( argc != 6 ){
		std::cerr << "Rigid registration of two volumes by local maximixzation of their cross correlation.\n" << endl;
		std::cerr << "Usage: " << argv[0] << " input1 input2 shift_tolerance rot_tolerance output_vol\n";
		std::cerr << "  1. input1: mrc volume to be used as reference." << endl;
		std::cerr << "  2. input2: mrc volume to be aligned to the reference." << endl;
		std::cerr << "  3. shift_tolerance: maximum shift allowable." << endl;
		std::cerr << "  4. rot_tolerance: maximum rotation allowable." << endl;
		//std::cerr << "  4. binning factor: binning to compute the cross-correlation." << endl;
		std::cerr << "  5. output_vol: input2 aligned to input1." << endl;
		return 0;
	}

	nbfImageFilter< PIXEL, 3 > imfilter;
	imfilter.paddingOn(1);
	imfilter.squareMaskSize = 2;

	nbfFourierFilter< PIXEL, 3 > fffilter;
	fffilter.bandPassOn(.01,.5,.01,.05); // groel

	PIXEL restrictTranslationSearch = atof( argv[3] );
	PIXEL restrictRotationSearch = atof( argv[4] );

	// Metric to be used for minimization
	nbfProjectionRotationMetric3D< PIXEL > metric( &imfilter, &fffilter );
	metric.setTranslationSearchRestriction( restrictTranslationSearch );
	metric.setRotationSearchRestriction( restrictRotationSearch );
	metric.setToComputeOverlapNormalizedDistances( false );
	metric.setToUseMutualCorrelation( true );
	metric.setMissingWedgeCompensation( false );

	char * input1FileName = argv[1];
	char * input2FileName = argv[2];

	vtkImageData * data = vtkImageData :: New();
	
	nbfMrcReader reader;
	reader.setFileName( input1FileName );
	reader.read( data );

	//PIXEL mag = atof( argv[5] );

	//vtkImageResample * resample = vtkImageResample :: New();
	//resample->SetInput( data );
	//resample->SetDimensionality( 3 );
	//resample->SetAxisMagnificationFactor( 0, 1 / mag );
	//resample->SetAxisMagnificationFactor( 1, 1 / mag );
	//resample->SetAxisMagnificationFactor( 2, 1 / mag );
	//resample->Update();

	nbfWedgedSubImage3D< PIXEL > im1;
	im1.setFixedImage( data );

	TinyVector< int, 3 > original_size = im1.getDimensions();

	reader.setFileName( input2FileName );
	reader.read( data );

	//resample->SetInput( data );
	//resample->Update();

	nbfWedgedSubImage3D< PIXEL > im2;
	im2.setFixedImage( data );

	int volumeSizeX = min( im1.getDimensions()[0], im2.getDimensions()[0] );
	if ( ( volumeSizeX % 2 ) != 0 ) volumeSizeX--;
	int volumeSizeY = min( im1.getDimensions()[1], im2.getDimensions()[1] );
	if ( ( volumeSizeY % 2 ) != 0 ) volumeSizeY--;
	int volumeSizeZ = min( im1.getDimensions()[2], im2.getDimensions()[2] );
	if ( ( volumeSizeZ % 2 ) != 0 ) volumeSizeZ--;

	imfilter.setMaskSize( volumeSizeX / 2 - 5, volumeSizeY / 2 - 5, volumeSizeZ / 2 - 5, 2 );

	TinyVector< int, 3 > tsize( volumeSizeX, volumeSizeY, volumeSizeZ );
	im1.setDimensions( tsize );
	im2.setDimensions( tsize );
	//im1.setCutOffset( tsize[2] / 2 );
	//im2.setCutOffset( tsize[2] / 2 );

	metric.setInput1( reinterpret_cast< nbfWedgedImage3D<PIXEL>* >( &im1 ) );
	metric.setInput2( reinterpret_cast< nbfWedgedImage3D<PIXEL>* >( &im2 ) );

	// vtkTransform * t = vtkTransform::New();

	metric.getDistance();
	
	// // descent refinement
	// if ( restrictRotationSearch < 0 ){
		// metric.refine(t,0);
	// } else {
		// metric.setRotationSearchRestriction( restrictRotationSearch );
		// metric.refine(t);
	// }

	 vtkTransform * t = metric.getTransform();
	 
	 cout << "Distance = " << metric.getCorrelationPeak() << endl;
	 cout << "Transformation = " << *t->GetMatrix() << endl;
	
	//double matrix[16];
	//vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );

	//matrix[3] *= mag;
	//matrix[7] *= mag;
	//matrix[11] *= mag;

	//vtkTransform * tmag = vtkTransform :: New();
	//tmag->GetMatrix()->DeepCopy( matrix );

	nbfMrcWriter mrcw;
	mrcw.setFileName( argv[5] );

	//cout << "Mag transformation = " << *tmag->GetMatrix() << endl;

	//reader.read( data );
	//im2.setFixedImage( data );
	im2.setDimensions( original_size );
	im2.getImage( data, t );

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToFloat();
	cast->SetInput( data );
	cast->Update();

	mrcw.write( cast->GetOutput() );

	return 0;
}