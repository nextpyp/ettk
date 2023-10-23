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

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfBlitzWriter.h>
#include <io/nbfBlitzReader.h>
#include <io/nbfMrcWriter.h>
#include <nbfVeselnessFilter.h>
#include <bs/nbfBordStrategyMirror.h>

#include <em/nbfCutSubVolumes.h>

void main( int argv, char ** argc )
{
	// argc[1] - input volume
	// argc[2] - input volume (full resolution)
	// argc[3] - input implicit surface
	// argc[4] - surface offset

	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName(argc[1]);
	reader->Update();

	// read 3D image data
	vtkImageData * myImage = vtkImageData::New();
	myImage->DeepCopy( reader->GetOutput() );

	// read 3D image data (full resolution)
	reader->SetFileName(argc[2]);
	reader->Update();

	vtkImageData * myImageFullResolution = vtkImageData::New();
	myImageFullResolution->DeepCopy( reader->GetOutput() );

	// read 3D implicit surface
	reader->SetFileName(argc[3]);
	reader->Update();

	vtkImageData * mySurface = reader->GetOutput();

	Array< float, 3 > flow;
	nbfVTKInterface::vtkToBlitz( mySurface, flow );
	flow = flow - atof(argc[4]);

	cout << flow.shape() << endl;

	double range[2];
	mySurface->GetPointData()->GetScalars()->GetRange(range);
	cout << range[0] << ", " << range[1] << endl;

	// specify cropped region size
	int size[3];
	size[0] = 30; // x
	size[1] = 30; // y
	size[2] = 30; // z (height)

	// build spherical generic template
	Array< short, 3 > myCropBlitz( size[0], size[1], size[2] );
	
	firstIndex i; secondIndex j; thirdIndex k;
	myCropBlitz = sqrt( pow2( 1.0 * i - size[0] / 2.0 ) + pow2( 1.0 * j - size[1] / 2.0 )  + pow2( 1.0 * k - size[2] / 2.0 ) );

	// build spherical generic template
	Array< short, 3 > myCropFullResolutionBlitz( 2.5*size[0], 2.5*size[1], 2.5*size[2] );

	nbfCutSubVolumes< float > cutter;
	cutter.setTemplate( myCropBlitz );
	cutter.setSurface( mySurface );
	cutter.setImage( myImage );
	cutter.setFullResolutionImage( myImageFullResolution );
	cutter.setTemplateFullResolution( myCropFullResolutionBlitz );
	cutter.setCorrelationTh(-1.0);
	cutter.setOffsetUnderSurface(0.0);
	cutter.setOffsetUnderSurfaceFullResolution(5.0);
	//cutter.setSurfaceLbound(-.5);
	//cutter.setSurfaceUbound(0);
	cutter.setSurfaceLbound(-0.05);
	cutter.setSurfaceUbound(0.05);
	cutter.setMinimumSeparation(20);
	cutter.setHeightRange(thirdDim,55,125);
	cutter.setNotToIgnoreValuesBelowSurface();

	vector< Array< short, 3 > > volumes;
	vector< vtkTransform * > transforms;
	cutter.execute( volumes, transforms );

	//cout << volumes.size() << endl;
	
	//nbfMrcWriter mrwriter;
	//mrwriter.setFileName("test.mrc");
	//mrwriter.write( myImage );

	return;

	vtkImageData * output = vtkImageData::New();

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	for ( int i = 1; i <= volumes.size(); i++ ){
		stringstream fileName;
		if ( i < 10 ){
			fileName << argc[5] << "_0" << i << ".vtk";
		}
		else{
			fileName << argc[5] << "_" << i << ".vtk";
		}
		writer->SetFileName(fileName.str().c_str());
		nbfVTKInterface converter;
		//cout << volumes[i-1].shape() << endl;
		//transforms[i]->Print(cout);
		converter.blitzToVtk( volumes[i-1], output );
		writer->SetInput( output );
		writer->Write();
		cout << "File: " << fileName.str() << " written." << endl;
	}


	// Store wedge limits and orientation
	/////////////////////////////////////
	// Values are stored in a 2D blitz array.
	// There is one row for each volume.
	// Columns (8) are as follows:
	// wedgeL wedgeU wedgeAngleX wedgeAngleY wedgeAngleZ positionX positionY positionZ

	double lWedge = -59.99;
	double uWedge = 59.99;
	Array< double, 2 > save( transforms.size(), 8 );
	save = 0;
	for ( int i = 0; i < save.rows(); i++ ){
		save(i,0) = lWedge;
		save(i,1) = uWedge;
		double orient[3];
		transforms[i]->GetOrientation(orient);
		for ( int j = 0; j <= 2; j++ ){
			save(i,2+j) = orient[j];
		}
		double position[3];
		transforms[i]->GetPosition(position);
		for ( int j = 0; j <= 2; j++ ){
			save(i,5+j) = position[j];
		}
	}
	
	nbfMatlabWriter bwriter;
	bwriter.setFileName("geometry.matlab");
	bwriter.write(save);

	Array< float, 2 > K;
	nbfBlitzReader breader;
	breader.setFileName("geometry.blitz");
	breader.read(K);
	cout << K << endl;

	writer->Delete();
	myImage->Delete();
	reader->Delete();
}