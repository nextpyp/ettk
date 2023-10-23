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
#include <vtkCylindricalTransform.h>
#include <vtkImageMedian3D.h>
#include <vtkImageButterworthHighPass.h>
#include <vtkImageConstantPad.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMrcWriter.h>

#include <nbfTimer.h>

#include <vtkPNGReader.h>

#include <em/nbfReconstruction3D.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <random/normal.h>
#include <random/uniform.h>
#include <time.h>

#define PIXEL float

void main( int argc, char ** argv  )
{
	if ( argc != 10 ){
		cout << "Usage: input angleX angleY angleZ wedgeLow wedgeHigh noiseAmplitude shiftAmplitude output" << endl;
		exit(0);
	}
	
	char * fileIn = argv[1];

	double angleX = atof(argv[2]);
	double angleY = atof(argv[3]);
	double angleZ = atof(argv[4]);

	double wedgeLow = atof(argv[5]);
	double wedgeHigh = atof(argv[6]);

	double noiseAmplitude = atof(argv[7]);

	double shiftAmplitude = atof(argv[8]);

	char * fileOut = argv[9];

	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	vtkImageData * inputData = vtkImageData::New();

	vtkStructuredPointsReader * reader = NULL;

	// read input data (VTK format)
	//reader = vtkStructuredPointsReader::New();
	//reader->SetFileName( fileIn );
	//reader->Update();
	//inputData->DeepCopy( reader->GetOutput() );

	// read input data (MRC format)
	nbfMrcReader rmrc;
	rmrc.setFileName(fileIn);
	rmrc.read( inputData );

	int extent[6];
	inputData->GetExtent(extent);

	double size = 63;
	vtkImageResample * resample = vtkImageResample::New();
	resample->SetInput( inputData );
	resample->SetAxisMagnificationFactor( 0, size / extent[1] );
	resample->SetAxisMagnificationFactor( 1, size / extent[3] );
	resample->SetAxisMagnificationFactor( 2, size / extent[5] );
	resample->Update();
	inputData->DeepCopy( resample->GetOutput() );

	// pad output so final reconstruction does not have gaps coming from the circular reconstruction domain
	vtkImageConstantPad * pad = vtkImageConstantPad::New();
	pad->SetConstant( inputData->GetScalarComponentAsDouble(0,0,0,0) );
	pad->SetInput( inputData );

	inputData->GetExtent(extent);
	int padding = extent[1] * ( sqrt(2.0) - 1.0 );

	pad->SetOutputWholeExtent( extent[0] - padding, extent[1] + padding, extent[2], extent[3], extent[4] - padding, extent[5] + padding );
	pad->Update();

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToDouble();
	cast->SetInput( pad->GetOutput() );

	vtkImageChangeInformation * center = vtkImageChangeInformation::New();
	center->CenterImageOn();
	center->SetInput( cast->GetOutput() );
	center->Update();

	vtkImageReslice * reslice = vtkImageReslice::New();
	reslice->WrapOff();
	reslice->MirrorOff();
	reslice->SetBackgroundLevel(0);
	reslice->SetInput( center->GetOutput() );
	reslice->SetInterpolationModeToCubic();

	vtkTransform * transform = vtkTransform::New();
	transform->RotateZ( angleX );
	transform->RotateY( angleY );
	transform->RotateZ( angleZ );
	transform->Translate(0,0,0);
	reslice->SetResliceTransform( transform );
	reslice->Update();

	firstIndex i;
	Array< PIXEL, 1 > angles( ceil( ( wedgeHigh - wedgeLow ) / 3.0 ) );
	angles = wedgeHigh * i / (PIXEL)angles.ubound(0) + wedgeLow * ( 1.0 - i / (PIXEL)angles.ubound(0) );

	//cout << angles << endl;

	//w.write(angles);
	angles *= vtkMath::DegreesToRadians();

	Array< PIXEL, 3 > V;
	nbfVTKInterface::vtkToBlitz( reslice->GetOutput(), V );

	Array< short, 3 > Vs;
	nbfVTKInterface::vtkToBlitz( reslice->GetOutput(), Vs );
	//w.write(Vs);

	Array< PIXEL, 3 > projections( V.rows(), V.cols(), angles.size() );

	Array< PIXEL, 3 > P;

	// Instantiate Radon object to compute projections
	nbfReconstruction3D< PIXEL > rec3Dp;
	rec3Dp.setAngles( angles );
	rec3Dp.setProjections( projections );

	ranlib::Uniform<float> noise;
    noise.seed((unsigned int)time(0));
	
	// store 3D reconstruction
	Array< PIXEL, 3 > reconstruction( projections.rows(), projections.cols(), projections.rows() );
	reconstruction = V;
	// reconstruction = 0;
	rec3Dp.setImage( reconstruction );

	rec3Dp.getProjections( P );

	w.write(P);

	// simulate image shifts as tilt series alignment errors
	for ( int i = 0; i < projections.depth(); i++ ){
		vtkImageData * slice = vtkImageData :: New();
		nbfVTKInterface :: blitzToVtk( P( Range::all(), Range::all(), i ), slice );
		vtkTransform * t = vtkTransform :: New();
		t->Translate( noise.random() * shiftAmplitude, noise.random() * shiftAmplitude, 0 );
		vtkImageReslice * reslice = vtkImageReslice :: New();
		reslice->SetInput( slice );
		reslice->SetInterpolationModeToCubic();
		reslice->SetBackgroundLevel( 0.0 );
		reslice->SetResliceTransform( t );
		reslice->Update();
		Array< PIXEL, 2 > shiftedProjection;
		nbfVTKInterface :: vtkToBlitz( reslice->GetOutput(), shiftedProjection );
		P( Range::all(), Range::all(), i ) = shiftedProjection;
		slice->Delete();
		t->Delete();
		reslice->Delete();
	}

	w.write(P);

	// add noise to projections
	Array< PIXEL, 3 > :: iterator iter = P.begin();
	Array< PIXEL, 3 > :: iterator iterP = projections.begin();
	while ( iter != P.end() ){
		(*iter) = (*iter) + noiseAmplitude * ( noise.random() - .5 );
		(*iterP) = (*iter);
		++iter; ++iterP;
	}
	//w.write(projections);

	// Instantiate second Radon object to compute reconstruction
	nbfReconstruction3D< PIXEL > rec3D;
	rec3D.setAngles( angles );
	rec3D.setProjections( projections );
	rec3D.setImage( reconstruction );

	reconstruction = 0;

	//rec3D.sirt(50);
	rec3D.art(5);

	reconstruction -= min(reconstruction);
	reconstruction = reconstruction / max(reconstruction) * 2.0 * numeric_limits<short>::max() - numeric_limits<short>::max();

	Array< PIXEL, 3 > C( reconstruction.shape() );
	//C = reconstruction.transpose( thirdDim, secondDim, firstDim );
	//reconstruction = C.reverse(secondDim);

	// crop center part to keep only area inside reconstruction circle
	C.resize( extent[1] + 1, extent[3] + 1, extent[5] + 1 );
	Range I( padding, padding + extent[1] );
	C = reconstruction( I, Range::all(), I );

	vtkImageData * filtered = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( C, filtered );

	//w.write(C);

	//vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	//writer->SetInput( filtered );
	//writer->SetFileName( fileOut );
	//writer->Write();
	//writer->Delete();

	vtkImageCast * cast2 = vtkImageCast::New();
	cast2->SetOutputScalarTypeToShort();
	cast2->SetInput( filtered );
	cast2->Update();

	nbfMrcWriter mrc;
	mrc.setFileName( fileOut );
	mrc.write( cast2->GetOutput() );

	filtered->Delete();
	center->Delete();
	reslice->Delete();
	transform->Delete();
	if ( reader != NULL ){
		reader->Delete();
	}
	cast->Delete();
	cast2->Delete();
	resample->Delete();
	inputData->Delete();
}