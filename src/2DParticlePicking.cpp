#include "mpi.h"
// #define NBF_VERBOSE

#define BZ_GENERATE_GLOBAL_INSTANCES

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>
#include <blitz/array/convolve.h>

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

#include "itkImage.h"
#include "itkNormalizedCorrelationImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageKernelOperator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"

typedef double PIXEL;

typedef itk::Image<unsigned char, 2> UnsignedCharImageType;

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

	PIXEL m = mean(A);
	PIXEL v = sqrt( mean( pow2( A - m ) ) );

	Array< PIXEL, 2 > C( A.shape() );
	C = where( A < m - 2.5*v || A > m + 2.5*v, m, A );

	cout << "Mean, Variance = " << m << ", " << v << endl;
	
	Array< PIXEL, 2 > T( A(Range(1600,1700),Range(1930,2030)) );

	typedef itk::Image< PIXEL, 2 > ImageType;
	typedef itk::VTKImageToImageFilter<ImageType> vtk2itk;
	typedef itk::ImageToVTKImageFilter<ImageType> itk2vtk;
	
	vtk2itk::Pointer image = vtk2itk::New();
	image->SetInput( data );	

	vtk2itk::Pointer reference = vtk2itk::New();
	reference->SetInput( data );	

	// Perform normalized correlation
	// <input type, mask type (not used), output type>
	typedef itk::NormalizedCorrelationImageFilter<ImageType, ImageType, ImageType> CorrelationFilterType;
       
	itk::ImageKernelOperator< PIXEL > kernelOperator;
	kernelOperator.SetImageKernel( reference->GetOutput() );
       
	// The radius of the kernel must be the radius of the patch, NOT the size of the patch
	itk::Size<2> radius = reference->GetOutput()->GetLargestPossibleRegion().GetSize();
	radius[0] = (radius[0]-1) / 2;
	radius[1] = (radius[1]-1) / 2;
       
	kernelOperator.CreateToRadius(radius);
              
	CorrelationFilterType::Pointer correlationFilter = CorrelationFilterType::New();
	correlationFilter->SetInput( image->GetOutput() );
	correlationFilter->SetTemplate( kernelOperator );
	correlationFilter->Update();

	itk2vtk::Pointer result = itk2vtk::New();
	result->SetInput( correlationFilter->GetOutput() );
	
	nbfVTKInterface::vtkToBlitz( result->GetOutput(), A );
	
	typedef itk::MinimumMaximumImageCalculator <ImageType> MinimumMaximumImageCalculatorType;
	 
	MinimumMaximumImageCalculatorType::Pointer minimumMaximumImageCalculatorFilter = MinimumMaximumImageCalculatorType::New ();
	minimumMaximumImageCalculatorFilter->SetImage(correlationFilter->GetOutput());
	minimumMaximumImageCalculatorFilter->Compute();
       
	itk::Index<2> maximumCorrelationPatchCenter = minimumMaximumImageCalculatorFilter->GetIndexOfMaximum();
	std::cout << "Maximum: " << maximumCorrelationPatchCenter << std::endl;
       
	// Note that the best correlation is at the center of the patch we extracted (ie. (75, 75) rather than the corner (50,50)
       
	typedef itk::RescaleIntensityImageFilter< ImageType, UnsignedCharImageType > RescaleFilterType;
	typedef itk::ImageFileWriter<UnsignedCharImageType> WriterType;
	{
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(correlationFilter->GetOutput());
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->Update();
       
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(rescaleFilter->GetOutput());
	writer->SetFileName("correlation.png");
	writer->Update();
	}
       
	{
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(extractFilter->GetOutput());
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->Update();
       
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(rescaleFilter->GetOutput());
	writer->SetFileName("patch.png");
	writer->Update();
	}
       
	// Extract the best matching patch
	FloatImageType::IndexType bestPatchStart;
	bestPatchStart[0] = maximumCorrelationPatchCenter[0] - radius[0];
	bestPatchStart[1] = maximumCorrelationPatchCenter[1] - radius[1];
       
	FloatImageType::RegionType bestPatchRegion(bestPatchStart,patchSize);
       
	ExtractFilterType::Pointer bestPatchExtractFilter = ExtractFilterType::New();
	bestPatchExtractFilter->SetRegionOfInterest(bestPatchRegion);
	bestPatchExtractFilter->SetInput(reader->GetOutput());
	bestPatchExtractFilter->Update();
	
	//PIXEL mask = ( A.rows() / 2 ) * .95;
	//PIXEL apod = ( A.rows() / 2 ) * .05;
	//nbfImageFilter< PIXEL, 2 > ifilter;
	//ifilter.setMaskSize( mask, mask, mask, apod );
	//ifilter.execute( A );
	//
	//Array< complex< PIXEL >, 2 > FFTA( A.shape() );
	//
	//nbfFourierFilter< PIXEL, 2 > ff;
	//ff.initializeFFTW( TinyVector< int, 2 >( A.shape() ) );
	//real( ff.blitzFFT1 ) = A * ff.shift;
	//imag( ff.blitzFFT1 ) = 0;
	//
	//fftw_execute( ff.fftplan );
	//
	//TinyVector<PIXEL,2> t(.1,.1);
	//ff.highPassOn(t);
	//ff.execute( ff.blitzFFT2 );
	//
	//A = sqrt( real( ff.blitzFFT2 * conj( ff.blitzFFT2 ) ) );

	Array< float, 3 > B( A.rows(), A.cols(), 1 );
	//A = A - min(A);
	//A = A / max(A);
	//A = log( 1 + 10 * A );
	B( Range :: all(), Range :: all(), 0 ) = cast<float>( A );
	nbfMrcWriter mrcw;
	mrcw.setFileName( output );
	mrcw.write( B );

	data->Delete();
}
