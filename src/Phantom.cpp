/*
 * Phantom.cpp
 *
 *  Created on: Jan 31, 2010
 *      Author: fefo
 */

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
#include <vtkSphereSource.h>
#include <io/nbfVTKInterface.h>
#include <io/nbfMrcWriter.h>

int main( int argc, char *argv[] ) {

	Array< float, 3 > T( 128, 128, 2 );
	T = 0;
	T( 64, 63, Range :: all() ) = 1;
	//T( Range(127,128), Range(127,128), Range(127,128) ) = 1;
	nbfMrcWriter wr;
	wr.setFileName("tiltseries.mrc");
	wr.write(T);
	return 0;


	int membraneRadius;
	int membraneThickness;
	int volumeSize[3];
	double volumeCenter[3];
	int trimerSize[3];
	int numberOfTrimers;
	int monomerAngle = 35;
	double smoothRadiusFactor = 2.0;
	double smoothStandardDeviation = 1.0;

	// Default values
	volumeSize[0] = 256/2;
	volumeSize[1] = 256/2;
	volumeSize[2] = 128/2;
	membraneRadius = 40/2;
	membraneThickness = 7;
	trimerSize[0] = 3;
	trimerSize[1] = 3;
	trimerSize[2] = 9;

	// Read the number of spikes per virion.
	numberOfTrimers = atoi( argv[2] );

	// Process the input parameters.
	for (int i=1; i<argc; i++) {
		if ( !strcmp(argv[i],"-mR") ) membraneRadius = atoi(argv[++i]);
		if ( !strcmp(argv[i],"-mT") ) membraneThickness = atoi(argv[++i]);
		if ( !strcmp(argv[i],"-mA") ) monomerAngle = atoi(argv[++i]);
		if ( !strcmp(argv[i],"-sRF") ) smoothRadiusFactor = atof(argv[++i]);
		if ( !strcmp(argv[i],"-sSD") ) smoothStandardDeviation = atof(argv[++i]);
		if ( !strcmp(argv[i],"-vS") ) {
			volumeSize[0] = atoi(argv[++i]);
			volumeSize[1] = atoi(argv[++i]);
			volumeSize[2] = atoi(argv[++i]);
		}
		if ( !strcmp(argv[i],"-tS") ) {
			trimerSize[0] = atoi(argv[++i]);
			trimerSize[1] = atoi(argv[++i]);
			trimerSize[2] = atoi(argv[++i]);
		}
	}
	cout << "Membrane Radius: " << membraneRadius << endl;
	cout << "Membrane Thickness: " << membraneThickness << endl;
	cout << "Monomer Angle: " << monomerAngle << endl;
	cout << "Volume Size: (" << volumeSize[0] << "," << volumeSize[1] << "," << volumeSize[2] << ")" << endl;
	cout << "Trimer Size: (" << trimerSize[0] << "," << trimerSize[1] << "," << trimerSize[2] << ")" << endl;
	cout << "Smoothing Radius Factor " << smoothRadiusFactor << endl;
	cout << "Smoothing Standard Deviation: " << smoothStandardDeviation << endl;

	// Calculate the volume center.
	volumeCenter[0] = floor(volumeSize[0]/2);
	volumeCenter[1] = floor(volumeSize[1]/2);
	volumeCenter[2] = floor(volumeSize[2]/2);

	// Writer for vtk files.
	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();

	// Temporal volume for building the trimer.
//	vtkImageData * tmp = vtkImageData::New();
	vtkImageData * tmpTrimer = vtkImageData::New();
	vtkImageData * tmpVirus = vtkImageData::New();

	tmpTrimer->SetScalarTypeToFloat();
	tmpVirus->SetScalarTypeToFloat();

	// Operators.
	vtkImageMathematics * add = vtkImageMathematics::New();
	add->SetOperationToMax();
	vtkImageMathematics * substract = vtkImageMathematics::New();
	substract->SetOperationToAdd();

	// The virus (exterior) membrane.
	vtkImageEllipsoidSource * membraneEllipse = vtkImageEllipsoidSource::New();
	membraneEllipse->SetWholeExtent(0,volumeSize[0]-1,0,volumeSize[1]-1,0,volumeSize[2]-1);
	membraneEllipse->SetCenter(volumeCenter[0],volumeCenter[1],volumeCenter[2]);
	membraneEllipse->SetRadius(membraneRadius,membraneRadius,membraneRadius);
	membraneEllipse->SetInValue(1);
	membraneEllipse->SetOutValue(0);
	membraneEllipse->SetOutputScalarTypeToDouble();
	membraneEllipse->Update();


#if 1	
	
	// The virus (interior) membrane.
	vtkImageEllipsoidSource * membraneEllipse2 = vtkImageEllipsoidSource::New();
	membraneEllipse2->SetWholeExtent(0,volumeSize[0]-1,0,volumeSize[1]-1,0,volumeSize[2]-1);
	membraneEllipse2->SetCenter(volumeCenter[0],volumeCenter[1],volumeCenter[2]);
	membraneEllipse2->SetRadius(membraneRadius-membraneThickness,membraneRadius-membraneThickness,membraneRadius-membraneThickness);
	membraneEllipse2->SetInValue(-1);
	membraneEllipse2->SetOutValue(0);
	membraneEllipse2->SetOutputScalarTypeToDouble();
	membraneEllipse2->Update();

	// Compute the membrane.
	substract->SetInput1( membraneEllipse->GetOutput());
	substract->SetInput2( membraneEllipse2->GetOutput() );
	substract->Update();
	tmpVirus->DeepCopy( substract->GetOutput() );

	// First a monomer.
	vtkImageEllipsoidSource * monomer = vtkImageEllipsoidSource::New();
	monomer->SetWholeExtent(0,volumeSize[0]-1,0,volumeSize[1]-1,0,volumeSize[2]-1);
	monomer->SetCenter(volumeCenter[0],volumeCenter[1],volumeCenter[2]);
	monomer->SetRadius(trimerSize[0],trimerSize[1],trimerSize[2]);
	monomer->SetInValue(1);
	monomer->SetOutValue(0);
	monomer->SetOutputScalarTypeToDouble();
	monomer->Update();

	// Transformation to define the center of the volume.
	vtkImageChangeInformation * change1 = vtkImageChangeInformation::New();

	// Translate up the monomer.
	vtkImageReslice * r0 = vtkImageReslice::New();
	vtkTransform * t0 = vtkTransform::New();
	r0->SetInput( monomer->GetOutput() );
	t0->Translate(0,0,-trimerSize[2]);
	r0->SetResliceTransform( t0 );
	r0->SetBackgroundLevel( 0 );
	r0->Update();

	vtkImageReslice * r01 = vtkImageReslice::New();
	vtkImageReslice * r02 = vtkImageReslice::New();
	vtkImageReslice * r03 = vtkImageReslice::New();
	vtkTransform * t01 = vtkTransform::New();
	vtkTransform * t02 = vtkTransform::New();
	vtkTransform * t03 = vtkTransform::New();

	change1->SetInput( r0->GetOutput() );
	change1->SetOutputOrigin(-volumeCenter[0],-volumeCenter[1],-volumeCenter[2]);
	change1->Update();
	r01->SetInput( change1->GetOutput() );
	t01->RotateX( monomerAngle );
	r01->SetResliceTransform( t01 );
	r01->SetBackgroundLevel( 0 );
	r01->Update();

	r02->SetInput( change1->GetOutput() );
	t02->RotateX( monomerAngle );
	t02->RotateZ( 120 );
	r02->SetResliceTransform( t02 );
	r02->SetBackgroundLevel( 0 );
	r02->Update();

	r03->SetInput( change1->GetOutput() );
	t03->RotateX( monomerAngle );
	t03->RotateZ( 240 );
	r03->SetResliceTransform( t03 );
	r03->SetBackgroundLevel( 0 );
	r03->Update();

	add->SetInput1( r01->GetOutput() );
	add->SetInput2( r02->GetOutput() );
	add->Update();
	tmpTrimer->DeepCopy( add->GetOutput() );

	add->SetInput1( tmpTrimer );
	add->SetInput2( r03->GetOutput() );
	add->Update();
	tmpTrimer->DeepCopy( add->GetOutput() );

	// Translate the trimer to the top of the membrane.
	vtkImageReslice * rTrimer = vtkImageReslice::New();
	rTrimer->SetInput( tmpTrimer );
	vtkTransform * transformTrimer = vtkTransform::New();
	transformTrimer->Translate(0, 0, - membraneRadius );

	rTrimer->SetResliceTransform( transformTrimer );
	rTrimer->SetBackgroundLevel(0);
	rTrimer->Update();
	tmpTrimer->DeepCopy( rTrimer->GetOutput() );	// And ready to be used.

	// Initialize random seed
	srand ( time(NULL) );
	int angleX, angleZ;
	vtkImageChangeInformation * change2 = vtkImageChangeInformation::New();
	change2->SetOutputOrigin(-volumeCenter[0],-volumeCenter[1],-volumeCenter[2]);
	vtkImageChangeInformation * change3 = vtkImageChangeInformation::New();
	change3->SetOutputOrigin(-volumeCenter[0],-volumeCenter[1],-volumeCenter[2]);

	int arrayAngleX[3];
	int arrayAngleZ[12];
	arrayAngleX[0] = 90;
	arrayAngleX[1] = 90+20;
	arrayAngleX[2] = 90-30;
	arrayAngleZ[0] = 0;
	arrayAngleZ[1] = 30;
	arrayAngleZ[2] = 60;
	arrayAngleZ[3] = 90;
	arrayAngleZ[4] = 120;
	arrayAngleZ[5] = 150;
	arrayAngleZ[6] = 180;
	arrayAngleZ[7] = 210;
	arrayAngleZ[8] = 240;
	arrayAngleZ[9] = 270;
	arrayAngleZ[10] = 300;
	arrayAngleZ[11] = 330;

	vtkImageReslice * r4 = vtkImageReslice::New();
	vtkTransform * t4 = vtkTransform::New();
	int trimerCounter = 1;
	for (int indZ = 0; indZ < 12; indZ+=3 ) {
		for (int indX = 0; indX < 1; indX++ ) {

			if ( trimerCounter >= numberOfTrimers ){
				continue;
			}
			
			change2->SetInput( tmpTrimer );
			change2->Update();
			t4->Identity();	// Reset the transformation
			r4->SetInput( change2->GetOutput() );
			angleX =  arrayAngleX[ indX ];
			angleZ =  arrayAngleZ[ indZ ];
			t4->RotateX( angleX );
			t4->RotateWXYZ( angleZ, 0, 0, 1 );
			r4->SetResliceTransform( t4 );
			r4->SetBackgroundLevel(0);
			r4->Update();
			printf("[Spike %02d] Euler angles: (%3d,%3d)\n",trimerCounter++, angleX, angleZ);
//			cout << "Spike " << trimerCounter++ << ". Euler angles: (" << angleX << ", " << angleZ << ")" << endl;

			// Add the trimer to the membrane to make the virus.
			add->SetInput1( r4->GetOutput() );
			add->SetInput2( tmpVirus );
			add->Update();
			tmpVirus->DeepCopy( add->GetOutput() );

		}
	}
#endif

	vtkImageGaussianSmooth * smooth = vtkImageGaussianSmooth :: New();
	smooth->SetInput( tmpVirus );
	
	// if no trimers, assume we want to generate the segmentation for the membrane
	if ( numberOfTrimers == 0 ){
		smooth->SetInput( membraneEllipse->GetOutput() );
	}
	smooth->SetRadiusFactor( smoothRadiusFactor );
	smooth->SetStandardDeviation( smoothStandardDeviation);
	smooth->Update();

	double range[2];
	smooth->GetOutput()->GetPointData()->GetScalars()->GetRange( range );

	vtkImageMathematics * prod = vtkImageMathematics::New();
	prod->SetOperationToMultiplyByK();
	prod->SetConstantK( 10.0 / range[1] );
	prod->SetInput1( smooth->GetOutput() );
	prod->Update();

	// writer->SetFileName( argv[1] );
	// writer->SetInput( prod->GetOutput() );
	// writer->Write();

	Array< double, 3 > A;
	nbfVTKInterface::vtkToBlitz( prod->GetOutput(), A );
	Array< float, 3 > C( A.shape() );
	C = cast< float >(A);

	// set marker (single pixel) at fixed level
	//TinyVector< int, 3 > c = C.shape() / 2;
	//C( c ) = 100;
	
	nbfMrcWriter w;
	w.setFileName( argv[1] );
	w.write( C );
}
