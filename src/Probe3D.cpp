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
#include <vtkImageMathematics.h>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageCast.h>
#include <vtkPolyDataWriter.h>
#include <vtkProbeFilter.h>
#include <vtkContourFilter.h>
#include <vtkPolyData.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <nbfVeselnessFilter.h>
#include <bs/nbfBordStrategyMirror.h>

#include <nbf3Dalignment.h>

void main( int argv, char ** argc )
{
	vtkStructuredPointsReader * spreader = vtkStructuredPointsReader::New();
	spreader->SetFileName( argc[1] );
	spreader->Update();

	vtkImageData * myImage = vtkImageData::New();
	myImage->DeepCopy( spreader->GetOutput() );

	spreader->SetFileName( argc[2] );
	spreader->Update();

	vtkContourFilter * contour = vtkContourFilter::New();
	contour->SetInput( myImage );
	contour->SetValue(0,0);

	vtkProbeFilter * probe = vtkProbeFilter::New();
	probe->SetInput( vtkDataSet::SafeDownCast( (vtkObject*)contour->GetOutput() ) );
	probe->SetSource( spreader->GetOutput() );


	vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
	writer->SetFileName("surface.vtk");
	writer->SetInput( vtkPolyData::SafeDownCast( probe->GetOutput() ) );
	writer->Write();

	writer->Delete();
	probe->Delete();
	contour->Delete();
	myImage->Delete();
	spreader->Delete();
}