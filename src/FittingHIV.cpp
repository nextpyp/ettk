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

	//Array< float, 3 > A1;
	//nbfVTKInterface::vtkToBlitz( data, A1 );
	//reader.setFileName( argv[2] );
	//reader.read( data );
	//Array< float, 3 > B;
	//nbfVTKInterface::vtkToBlitz( data, B );
	//Array< float, 3 > M( A1.shape() );
	//M = where( A1 > .1, 1, 0 ); // CD4
	//// extend mask downwards
	//for ( int i = 26; i < 30; i++ ){
	//	M( Range :: all(), Range :: all(), i ) = M( Range :: all(), Range :: all(), 31 );
	//}

	//vtkImageDilateErode3D * erode = vtkImageDilateErode3D :: New();
	//erode->SetKernelSize( 10, 10, 10 );
	//erode->SetDilateValue( 1 );
	//erode->SetErodeValue( 1 );
	//nbfVTKInterface::blitzToVtk( M, data );
	//erode->SetInput( data );

	//M( Range :: all(), Range :: all(), Range(31,toEnd) ) = 1;

	//M = where( A1 > 0.072, 1, 0 ); // BaL
	
	//TinyVector< int, 3 > center(51.5,49.5,26);
	//M( Range(47,55), Range(45,53), Range(26,29) ) = 1;

	//// smooth
	//vtkImageGaussianSmooth * smooth = vtkImageGaussianSmooth :: New();
	//vtkImageData * cropVtk = vtkImageData :: New();
	//nbfVTKInterface::blitzToVtk( M, cropVtk );
	//smooth->SetInput( cropVtk );
	//smooth->SetDimensionality(3);
	//smooth->SetStandardDeviations(1.0,1.0,1.0);
	//smooth->SetRadiusFactors(1.1,1.1,1.1);
	//smooth->Update();
	//nbfVTKInterface::vtkToBlitz( smooth->GetOutput(), M );
	//smooth->Delete();
	//B = B * M;
	//nbfVTKInterface::blitzToVtk( B, data );

	vtkTransform * t = vtkTransform :: New();
	t->Translate(1,0,0);

	//vtkImageChangeInformation * change = vtkImageChangeInformation::New();
	//change->SetInput( data );
	//change->CenterImageOn();

	//// Bal-b12
	////t->RotateWXYZ( -(162.73455498773788-131.18108821140643), 0.0024221896540988579, 0.031779885837984065, 0.99949195687279258 );
	////t->Translate( -(-24.726823787524442 +23.446317246110464)/4.1, -(-0.34529322109741867-3.7646108157708937)/4.1, -(47.071834100450076-58.755668624358563)/4.1 );

	//// BaL
	////t->RotateWXYZ( 167.57886388152332-162.73455498773788, 0.01039495292332381, 0.032844901470697796, 0.99940640252157753);
	////t->Translate( -(-27.382705411123414+23.446317246110464+7)/4.1, -(-0.28386004214465749-3.7646108157708937)/4.1, -(44.032420615221142-58.755668624358563) /4.1 );

	//vtkImageReslice * reslice = vtkImageReslice :: New();
	//reslice->SetInput( data );
	//reslice->SetResliceTransform(t);
	//reslice->SetBackgroundLevel(0.0);
	//reslice->SetInterpolationModeToCubic();
	//reslice->Update();
	//
	//nbfMrcWriter writer1;
	//writer1.setFileName( argv[3] );
	//writer1.write( reslice->GetOutput() );
	//return 0;


	//Array< float, 3 > I;
	//nbfVTKInterface::vtkToBlitz( data, I );
	//float minI = min(I);
	////I = I - 1.28;
	////I = I - minI;
	//float maxI = max(I);
	////I = I / maxI;
	////I = I - ( 1.28 - minI ) / ( maxI - minI );
	////I = I - 1.28;

	//reader.setFileName( argv[2] );
	//reader.read( data );
	//Array< float, 3 > J;
	//nbfVTKInterface::vtkToBlitz( data, J );
	////J = where( J > 0.0624, J, 0 ); // B12
	////J = where( J > 0.0244, J, 0 ); // CD4
	////J = where( J > 0.0796, J, 0 ); // BaL


	////// smooth
	////vtkImageGaussianSmooth * smooth = vtkImageGaussianSmooth :: New();
	////vtkImageData * cropVtk = vtkImageData :: New();
	////nbfVTKInterface::blitzToVtk( J, cropVtk );
	////smooth->SetInput( cropVtk );
	////smooth->SetDimensionality(3);
	////smooth->SetStandardDeviations(1.0,1.0,1.0);
	////smooth->SetRadiusFactors(1.5,1.5,1.5);
	////smooth->Update();
	////nbfVTKInterface::vtkToBlitz( smooth->GetOutput(), J );
	////smooth->Delete();

	//float minJ = min(J);
	//J = J - minJ;
	//float maxJ = max(J);
	//J = J / max(J);
	//J = J * ( maxI - minI ) + minI;
	//// J = ( J - ( 0.0591 - minJ ) / ( maxJ - minJ ) ) * ( maxI - minI ) + 1.28;
	////J = J - 0.0591;

	////I = I + 1.35 * J;
	//I = I + .65 * J;
	//nbfVTKInterface::blitzToVtk( I, data );

	////cout << min(I) << endl;
	////cout << max(I) << endl;

	////vtkImageResample * resample = vtkImageResample :: New();
	////resample->SetInput( data );
	////resample->SetDimensionality(3);
	////resample->SetAxisMagnificationFactor( 0, 3 );
	////resample->SetAxisMagnificationFactor( 1, 3 );
	////resample->SetAxisMagnificationFactor( 2, 3 );
	////resample->SetInterpolationModeToCubic();
	////resample->Update();

	//////Array< float, 3 > I;
	////nbfVTKInterface::vtkToBlitz( resample->GetOutput(), I );
	//////I( 0, 18, Range :: all() ) = -1;
	//////I( I.ubound(firstDim), 18, Range :: all() ) = -1;
	//////I=-I;
	////nbfVTKInterface::blitzToVtk( I, data );
	////cout << min(I) << endl;
	////cout << max(I) << endl;

	//nbfMrcWriter writer1;
	//writer1.setFileName( argv[3] );
	//writer1.write( data );
	//return 0;

//#if 0
//	Array< float, 3 > cylinder( 100, 100, 100 );
//	firstIndex i;
//	secondIndex j;
//	cylinder = atan2( i - 50.0, j - 50.0 ) * vtkMath::RadiansToDegrees();
//
//	// Bal
//	cylinder = where( ( ( cylinder > -60 ) && ( cylinder < 60 ) ), 1, 0 );
//	cylinder( Range :: all(), Range :: all(), Range(0,35) ) = 0;
//
//	// Fab
//	//cylinder = where( ( ( cylinder > -20 ) && ( cylinder < 100 ) ), 1, 0 );
//	//cylinder( Range :: all(), Range :: all(), Range(0,41) ) = 0;
//	
//	// CD4
//	//	cylinder = where( ( ( cylinder > -35 ) && ( cylinder < 85 ) ), 1, 0 );
//	//  cylinder( Range :: all(), Range :: all(), Range(0,41) ) = 0;
//
//	vtkImageGaussianSmooth * smooth = vtkImageGaussianSmooth :: New();
//	vtkImageData * cropVtk = vtkImageData :: New();
//	nbfVTKInterface::blitzToVtk( cylinder, cropVtk );
//	smooth->SetInput( cropVtk );
//	smooth->SetDimensionality(3);
//	smooth->SetStandardDeviations(1.0,1.0);
//	smooth->SetRadiusFactors(1.0,1.0);
//	smooth->Update();
//	nbfVTKInterface::vtkToBlitz( smooth->GetOutput(), cylinder );
//	smooth->Delete();
//
//	Array< float, 3 > image;
//	nbfVTKInterface::vtkToBlitz( data, image );
//	image *= cylinder;
//	nbfVTKInterface::blitzToVtk( image, data );
//
//	vtkTransform * t = vtkTransform :: New();
//	// Bal
//	t->Translate(0,11,1);
//
//	// Fab
//	// t->Translate(9,7,1);
//
//	// CD4
//	// t->Translate(5,10,4);
//
//	vtkImageReslice * reslice = vtkImageReslice :: New();
//	reslice->SetInput( data );
//	reslice->SetResliceTransform(t);
//	reslice->SetBackgroundLevel(0.0);
//	reslice->SetInterpolationModeToCubic();
//	reslice->Update();
//
//	nbfMrcWriter writer1;
//	writer1.setFileName( argv[2] );
//	writer1.write( reslice->GetOutput()  );
//
//#else
//
//	vtkImageResample * resample = vtkImageResample::New();
//	double factor = 6.667 / 4.1;
//	factor = 1;
//
//	resample->SetDimensionality(3);
//	resample->SetInterpolationModeToCubic();
//	resample->SetAxisMagnificationFactor( 0, factor );
//	resample->SetAxisMagnificationFactor( 1, factor );
//	resample->SetAxisMagnificationFactor( 2, factor );
//	resample->SetInput( data );
//	resample->SetBackgroundLevel(0);
//	resample->Update();
//
//	vtkImageChangeInformation * center = vtkImageChangeInformation :: New();
//	center->CenterImageOn();
//	center->SetInput( resample->GetOutput() );
//	center->Update();
//
//	int extent[6];
//	resample->GetOutput()->GetExtent(extent);
//
//	vtkImageConstantPad * pad = vtkImageConstantPad :: New();
//	pad->SetInput( center->GetOutput() );
//	int size[3];
//	size[0] = ( 64 - extent[1] - extent[0] ) / 2.0;
//	size[1] = ( 64 - extent[3] - extent[2] ) / 2.0;
//	size[2] = ( 64 - extent[5] - extent[4] ) / 2.0;
//	pad->SetOutputWholeExtent( extent[0] - size[0], extent[1] + size[0] - 2, extent[2] - size[1], extent[3] + size[1] - 2, extent[4] - size[2], extent[5] + size[2] - 2 );
//	pad->SetConstant(0);
//	pad->Update();
//
//	nbfMrcWriter writer;
//	writer.setFileName( argv[2] );
//	writer.write( pad->GetOutput() );
//	
//	data->Delete();
//	resample->Delete();
//	center->Delete();
//#endif

	// Setup image filter
	nbfImageFilter< PIXEL, 3 > imfilter;
	imfilter.setMaskSize( 32, 32, 32, 4, false );
	stringstream maskFile;
	maskFile << argv[1];
	float th = atof( argv[3] );
	imfilter.setMaskFile( maskFile, th );

	// Setup Fourier filter
	nbfFourierFilter< PIXEL, 3 > fffilter;
	fffilter.bandPassOn( 0, .5, .005, .01 );

	// Setup metric configuration
	nbfProjectionRotationMetric3D< PIXEL > metric( &imfilter, &fffilter );

	metric.setNumberOfCandidatePeaksToSearch( 50 );

	metric.setMissingWedgeCompensation( false );
	metric.setRotationSearchRestriction( 0.0 );
	metric.setTranslationSearchRestriction( 10.0 );

	metric.setNumberOfCandidates( 1 );
	metric.setToUseMutualCorrelation( false );
	metric.setToComputeOverlapNormalizedDistances( false );

	nbfWedgedSubImage3D< float > vol1, vol2;
	Array< float, 3 > A;
	nbfVTKInterface :: vtkToBlitz( data, A );

	vol1.setFixedImage( A );

	reader.setFileName( argv[2] );
	reader.read( data );
	nbfVTKInterface :: vtkToBlitz( data, A );
	vol2.setFixedImage( A );

	metric.setInput1( &vol1 );
	metric.setInput2( &vol2 );
	metric.getDistance();
}