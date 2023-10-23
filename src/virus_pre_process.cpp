#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMrcReader.h>
#include <io/nbfMrcWriter.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfVTKInterface.h>
#include <nbfMaximalFlow.h>
#include <nbfDifferentials.h>

#include <vtkImageCast.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageMathematics.h>
#include <vtkMath.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>

#include <vtkPoints.h>
#include <vtkImageData.h>

#define PIXEL float

void main( int argc, char ** argv )
{

	if ( argc != 5 ){
		cout << "USAGE: virus_pre_process.exe input binning smooth output" << endl;
	}

	// read image volume
	nbfMrcReader mrc;
	mrc.setFileName( argv[1] );
	vtkImageData * g = vtkImageData::New();
	mrc.read(g);

	float binning = 1 / atof( argv[2] );
	float smooth = atof( argv[3] );

	vtkImageResample * resample = vtkImageResample :: New();
	resample->SetAxisMagnificationFactor( 0, binning );
	resample->SetAxisMagnificationFactor( 1, binning );
	resample->SetAxisMagnificationFactor( 2, binning );
	resample->SetInput( g );

	vtkImageGaussianSmooth * smoothf = vtkImageGaussianSmooth :: New();
	smoothf->SetDimensionality(3);
	smoothf->SetRadiusFactors( smooth, smooth, smooth);
	smoothf->SetInput( resample->GetOutput() );
	smoothf->Update();

	// store in file

	// fix geometry for MRC writer
	//Array< float, 3 > A, B;
	//nbfVTKInterface::vtkToBlitz( smoothf->GetOutput(), A );
	//A.transposeSelf(thirdDim,secondDim,firstDim);
	//B.resize( A.shape() );
	//B = A.reverse(secondDim);
	//nbfVTKInterface::blitzToVtk( B, smoothf->GetOutput() );

	nbfMrcWriter writerm;
	writerm.setFileName( argv[4] );
	writerm.write( smoothf->GetOutput() );
	cout << "File: " << argv[4] << " written." << endl;

	smoothf->Delete();
	resample->Delete();
	g->Delete();
}