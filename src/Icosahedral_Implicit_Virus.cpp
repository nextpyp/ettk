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

#include <random/normal.h>
#include <random/uniform.h>
#include <time.h>

void main( int argc, char ** argv )
{
	vtkImageEllipsoidSource * ellipsoid = vtkImageEllipsoidSource::New();
	ellipsoid->SetWholeExtent(0,127,0,127,0,127);
	ellipsoid->SetCenter(64,64,64);
	ellipsoid->SetRadius(20,8,8);
	ellipsoid->SetInValue(1);
	ellipsoid->SetOutValue(0);
	ellipsoid->SetOutputScalarTypeToDouble();
	ellipsoid->Update();
	
	vtkImageChangeInformation * change1 = vtkImageChangeInformation::New();
	change1->CenterImageOn();
	// change1->SetOriginTranslation(-64,-64,-64);
	change1->SetInput( ellipsoid->GetOutput() );

	float phi = ( 1.0 + sqrt( 5.0 ) ) / 2.0;

	vtkImageReslice * r1 = vtkImageReslice::New();
	r1->SetInput( change1->GetOutput() );
	
	double axis1[3], axis2[3], axis[3];

	axis1[0] = 1; axis1[1] = axis1[2] = 0;

	vtkImageEllipsoidSource * sphere = vtkImageEllipsoidSource::New();
	sphere->SetWholeExtent(0,127,0,127,0,127);
	sphere->SetCenter(64,64,64);
	sphere->SetRadius(30,30,30);
	sphere->SetInValue(1);
	sphere->SetOutValue(0);
	sphere->SetOutputScalarTypeToDouble();
	sphere->Update();
	
	vtkImageData * tmpo = vtkImageData::New();
	r1->Update();
	tmpo->DeepCopy( sphere->GetOutput() );

	ranlib::Uniform<float> noise;
	noise.seed((unsigned int)time(0));
	Array< float, 1 > indexes(12);
	for ( int i = 0; i < indexes.size(); i++ ){
		indexes(i) = noise.random() > .75;
		indexes(i) = 1;
	}
	cout << indexes << endl;
	int count = 0;

	float factor = 30;
	for ( int coord = 0; coord < 3; coord++ ){
		int d1, d2;
		switch ( coord ){
			case 0:
				d1 = 1;
				d2 = 2;
				break;
			case 1:
				d1 = 2;
				d2 = 0;
				break;
			case 2:
				d1 = 0;
				d2 = 1;
				break;
		}

		for ( int j = -1; j < 2; j+=2 ){
			for ( int k = -1; k < 2; k+=2 ){

				if ( indexes(count) == 0 ){
					count++;
					continue;
				} else {
					count++;
				}

				axis2[coord] = 0;
				axis2[d1] = j;
				axis2[d2] = k * phi;

				cout << "t=[" << axis2[0] << "," << axis2[1] << "," << axis2[2] << "]" << endl;

				vtkMath :: Normalize( axis2 );

				vtkMath :: Cross( axis1, axis2, axis );

				vtkTransform * t1 = vtkTransform::New();
				
				float theta = vtkMath :: RadiansToDegrees() * acos( vtkMath :: Dot( axis1, axis2 ) );
				t1->RotateWXYZ( -theta, axis[0], axis[1], axis[2] );

				axis2[0] *= factor;
				axis2[1] *= factor;
				axis2[2] *= factor;
				t1->Translate( axis2 );

				r1->SetResliceTransform( t1 );
				r1->SetBackgroundLevel(0);
				r1->Update();

				vtkImageMathematics * add = vtkImageMathematics::New();
				add->SetOperationToMax();
				add->SetInput1( r1->GetOutput() );
				add->SetInput2( tmpo );
				add->Update();

				tmpo->DeepCopy( add->GetOutput() );
				add->Delete();
				t1->Delete();
			}
		}
	}

	//vtkImageReslice * r2 = vtkImageReslice::New();
	//r2->SetInput( change1->GetOutput() );
	//vtkTransform * t2 = vtkTransform::New();
	//t2->Translate( 0, 1, -phi);
	//r2->SetResliceTransform( t2 );
	//r2->SetBackgroundLevel(0);
	//r2->Update();
	//
	//add->SetInput2( tmpo );
	//add->SetInput1( r2->GetOutput() );
	//add->Update();

	//tmpo->DeepCopy( add->GetOutput() );

	//vtkImageReslice * r3 = vtkImageReslice::New();
	//r3->SetInput( tmpo );
	//vtkTransform * t3 = vtkTransform::New();
	//t3->Translate( 0, -1, phi);
	//r3->SetResliceTransform( t3 );
	//r3->SetBackgroundLevel(0);
	//r3->Update();

	//tmpo->DeepCopy( r3->GetOutput() );

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
	writer->SetFileName( argv[1] );
	writer->SetInput( prod->GetOutput() );
	writer->Write();

}