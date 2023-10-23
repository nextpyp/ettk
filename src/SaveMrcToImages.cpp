#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

#include <vtkPNGWriter.h>
#include <vtkImageData.h>
#include <vtkImageCast.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>

using namespace blitz;

#include <complex>

#include <nbfTimer.h>
#include <io/nbfVTKInterface.h>
#include <io/nbfMrcReader.h>

int main( int argc, char *argv[] )
{
	nbfMrcReader reader;
	reader.setBigEndian(true);
	reader.setFileName( argv[1] );

	vtkImageData * data = vtkImageData::New();
	Array< float, 1 > angles;	
	Array< float, 1 > means;

	reader.read( data, angles, means );

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName(argv[2]);
	writer->SetFileTypeToBinary();
	writer->SetInput( data );
	writer->Write();
}

