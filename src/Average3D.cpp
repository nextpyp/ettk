#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkImageData.h>
#include <vtkImageRFFT.h>
#include <vtkImageExtractComponents.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>

#include <io/nbfVTKInterface.h>

void main( int argc, char ** argv )
{
	if ( argc != 4 ){
		cout << "Usage: filepattern #volumes" << endl;
		exit(0);
	}

	char * pattern = argv[1];
	int numImages = atoi( argv[2] );
	int useWedge = atoi( argv[3] );
	
	// store average image
	Array< double, 3 > aReal;
	Array< double, 3 > aImag;
	Array< double, 3 > accum;
	
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	for ( int i = 1; i <= numImages; i++ ){
		stringstream fileName;
		fileName << pattern << i+1 << ".vtk";
		cout << fileName.str() << endl;
		reader->SetFileName( fileName.str().c_str() );
		reader->Update();
		Array< double, 3 > real;
		Array< double, 3 > imag;
		nbfVTKInterface::vtkToBlitz( reader->GetOutput(), real, imag );
		
		// initialize storage for average image
		if ( i == 1 ){
			aReal.resize( real.shape() );
			aImag.resize( real.shape() );
			accum.resize( real.shape() );
			aReal = 0;
			aImag = 0;
			accum = 0;
		}

		if ( useWedge == 1 ){
			accum = where( abs( real ) > 0, accum + 1, accum );
			aReal += real; 
			aImag += imag; 
		}
		else{
			aReal = ( ( i - 1.0 ) * aReal + real ) / ( i * 1.0 );
			aImag = ( ( i - 1.0 ) * aImag + imag ) / ( i * 1.0 );
		}
	}

	if ( useWedge == 1 ){
		aReal = where( accum > 0, aReal / accum, aReal );
		aImag = where( accum > 0, aImag / accum, aImag );
		cout << min(accum) << ", " << max(accum) << endl;
	}


	// convert result to VTK
	
	Array< double, 3 > real;
	Array< double, 3 > imag;

	// create new image data
	vtkImageData * average = vtkImageData::New();

	// grab geometry from reader
	average->DeepCopy( reader->GetOutput() );

	// take blitz references
	nbfVTKInterface::vtkToBlitz( average, real, imag );

	// operate on references
	real = aReal;
	imag = aImag;

	// now use the new changed values
	vtkImageRFFT * rfft = vtkImageRFFT::New();
	rfft->SetDimensionality(3);
	rfft->SetInput( average );
	rfft->Update();

	vtkImageExtractComponents * realp = vtkImageExtractComponents::New();
	realp->SetComponents(0);
	realp->SetInput( rfft->GetOutput() );
	realp->Update();

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	stringstream fileName;
	fileName << pattern << "average.vtk";
	writer->SetFileName( fileName.str().c_str() );
	writer->SetInput( realp->GetOutput() );
	writer->Write();

	rfft->Delete();
	realp->Delete();
	writer->Delete();
	average->Delete();
	reader->Delete();

}