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
#include <vtkImageShrink3D.h>
#include <vtkImageEllipsoidSource.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <nbfVeselnessFilter.h>
#include <bs/nbfBordStrategyMirror.h>

// #include <nbf3Dalignment.h>

void main( int argv, char ** argc )
{
	vtkImageEllipsoidSource * ellipsoid = vtkImageEllipsoidSource::New();
	ellipsoid->SetWholeExtent(0,63,0,63,0,63);
	ellipsoid->SetCenter(32,32,32);
	ellipsoid->SetRadius(4,4,15);
	ellipsoid->SetInValue(1);
	ellipsoid->SetOutValue(0);
	ellipsoid->SetOutputScalarTypeToDouble();
	ellipsoid->Update();
	
	vtkImageChangeInformation * change1 = vtkImageChangeInformation::New();
	change1->SetOriginTranslation(-32,-32,15-32);
	change1->SetInput( ellipsoid->GetOutput() );

	vtkImageReslice * r1 = vtkImageReslice::New();
	r1->SetInput( change1->GetOutput() );
	vtkTransform * t1 = vtkTransform::New();
	t1->RotateX(90);
	r1->SetResliceTransform( t1 );
	r1->SetBackgroundLevel(0);

	vtkImageMathematics * add = vtkImageMathematics::New();
	add->SetOperationToMax();
	add->SetInput1( r1->GetOutput() );
	add->SetInput2( ellipsoid->GetOutput() );
	add->Update();

	vtkImageData * tmpo = vtkImageData::New();
	tmpo->DeepCopy( add->GetOutput() );

	vtkImageReslice * r2 = vtkImageReslice::New();
	r2->SetInput( change1->GetOutput() );
	vtkTransform * t2 = vtkTransform::New();
	t2->RotateY(90);
	t2->RotateZ(-20);
	r2->SetResliceTransform( t2 );
	r2->SetBackgroundLevel(0);
	r2->Update();
	
	add->SetInput2( tmpo );
	add->SetInput1( r2->GetOutput() );
	add->Update();

	tmpo->DeepCopy( add->GetOutput() );

	vtkImageReslice * r3 = vtkImageReslice::New();
	r3->SetInput( tmpo );
	vtkTransform * t3 = vtkTransform::New();
	t3->RotateY(-45);
	t3->RotateX(-45);
	r3->SetResliceTransform( t3 );
	r3->SetBackgroundLevel(0);
	r3->Update();

	tmpo->DeepCopy( r3->GetOutput() );

	vtkImageGaussianSmooth * smooth = vtkImageGaussianSmooth :: New();
	smooth->SetInput( tmpo );
	smooth->SetRadiusFactor(10);
	smooth->Update();

	double range[2];
	smooth->GetOutput()->GetPointData()->GetScalars()->GetRange( range );

	vtkImageMathematics * prod = vtkImageMathematics::New();
	prod->SetOperationToMultiplyByK();
	prod->SetConstantK(10.0 / range[1]);
	prod->SetInput1( smooth->GetOutput() );
	prod->Update();

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName("trimer2.vtk");
	writer->SetInput( prod->GetOutput() );
	writer->Write();

}