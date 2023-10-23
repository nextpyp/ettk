#ifndef FILE_nbfDTIReader
#define FILE_nbfDTIReader

#include <vtkImageReader.h>

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>

#include <blitz/array.h>
using namespace blitz;

#include <io/nbfVTKInterface.h>
#include <io/nbfFile.h>

class nbfDTIReader : public nbfFile
{

public:

	// Read data into array. Return 1 if error, 0 otherwise.
	template< class Pixel >
	int read( Array< Pixel, 2 > & );

	template< class Pixel >
	int read( Array< Pixel, 3 > & );

	// read all tensor components
	template< class Pixel >
	int read( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< Pixel, 3 > &,
	          Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< Pixel, 3 > & );

	// read all tensor components
	template< class Pixel >
	int readRaw( Array< Pixel, 3 > &, int, int, int, float, float, float );

	// read all tensor components
	template< class Pixel >
	int readRaw( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< Pixel, 3 > &,
	             Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< Pixel, 3 > & );

	// read into 3D vtk volume
	template< class Pixel >
	int readRaw( vtkImageData * );

	// read into 3D vtk volume
	int read( vtkImageData * );

	nbfDTIReader(){};

};

template< class Pixel >
int nbfDTIReader :: read( Array< Pixel, 2 > & A ){
	return 0;
}

template< class Pixel >
int nbfDTIReader :: read( Array< Pixel, 3 > & A ){
	vtkImageData * output = vtkImageData::New();
	this->read( output );
	nbfVTKInterface::vtkToBlitz( output, A );
	output->Delete();
	return 0;
}

template< class Pixel >
int nbfDTIReader :: read( Array< Pixel, 3 > & Dxx, Array< Pixel, 3 > & Dxy, Array< Pixel, 3 > & Dxz,
						  Array< Pixel, 3 > & Dyy, Array< Pixel, 3 > & Dyz, Array< Pixel, 3 > & Dzz ){
	// initialize reader
	vtkImageReader * reader;
	reader = vtkImageReader::New();
	reader->SetDataByteOrderToBigEndian();
	reader->SetDataExtent(0, 125, 0, 159, 0, 46);
	reader->SetDataSpacing( -1.25, 1.25, 2.5);
	reader->SetDataScalarTypeToFloat();
	reader->SetFileDimensionality(3);
	
	// get all tensor components
	reader->SetFileName("Dxx.img");
	reader->Update();
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), Dxx );

	reader->SetFileName("Dxy.img");
	reader->Update();
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), Dxy );

	reader->SetFileName("Dxz.img");
	reader->Update();
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), Dxz );

	reader->SetFileName("Dyy.img");
	reader->Update();
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), Dyy );

	reader->SetFileName("Dyz.img");
	reader->Update();
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), Dyz );

	reader->SetFileName("Dzz.img");
	reader->Update();
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), Dzz );

	return 0;
}

template< class Pixel >
int nbfDTIReader :: readRaw( Array< Pixel, 3 > & D,
							 int rows, int cols, int depth,
							 float xspace = 1, float yspace = 1, float zspace = 1 )
{	
	// initialize reader
	vtkImageReader * reader;
	reader = vtkImageReader::New();
	reader->SetDataByteOrderToBigEndian();
	reader->SetDataExtent(0, rows-1, 0, cols-1, 0, depth-1);
	reader->SetDataSpacing( xspace, yspace, zspace );
	reader->SetDataScalarTypeToUnsignedShort();
	reader->SetFileDimensionality(3);
	reader->SetFileName(this->fileName);
	reader->Update();
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), D );
	return 0;
}

int nbfDTIReader :: read( vtkImageData * output ){

	// initialize reader
	vtkImageReader * reader;
	reader = vtkImageReader::New();
	reader->SetDataByteOrderToBigEndian();
	reader->SetDataExtent(0, 125, 0, 159, 0, 46);
	reader->SetDataSpacing( 1.25, 1.25, 2.5);
	reader->SetDataScalarTypeToFloat();
	reader->SetFileDimensionality(3);
	
	// first read coordinates + FA image
	reader->SetFileName("FA.img");
	reader->Update();

	// store full data set (3D coordinates + FA + tensor)
	output->DeepCopy( reader->GetOutput() );

	// store tensor values (3x3)
	vtkFloatArray * newTensors = vtkFloatArray::New();
	newTensors->SetNumberOfComponents(6);
	newTensors->SetNumberOfTuples( output->GetNumberOfPoints() );

	// get Dxx
	reader->SetFileName("Dxx.img");
	reader->Update();
	newTensors->CopyComponent( 0, reader->GetOutput()->GetPointData()->GetScalars(), 0 );

	// get Dxy
	reader->SetFileName("Dxy.img");
	reader->Update();
	newTensors->CopyComponent( 1, reader->GetOutput()->GetPointData()->GetScalars(), 0 );
	//newTensors->CopyComponent( 3, reader->GetOutput()->GetPointData()->GetScalars(), 0 );

	// get Dxz
	reader->SetFileName("Dxz.img");
	reader->Update();
	newTensors->CopyComponent( 2, reader->GetOutput()->GetPointData()->GetScalars(), 0 );
	//newTensors->CopyComponent( 6, reader->GetOutput()->GetPointData()->GetScalars(), 0 );

	// get Dyy
	reader->SetFileName("Dyy.img");
	reader->Update();
	newTensors->CopyComponent( 3, reader->GetOutput()->GetPointData()->GetScalars(), 0 );

	// get Dyz
	reader->SetFileName("Dyz.img");
	reader->Update();
	newTensors->CopyComponent( 4, reader->GetOutput()->GetPointData()->GetScalars(), 0 );
	//newTensors->CopyComponent( 7, reader->GetOutput()->GetPointData()->GetScalars(), 0 );

	// get Dzz
	reader->SetFileName("Dzz.img");
	reader->Update();
	newTensors->CopyComponent( 5, reader->GetOutput()->GetPointData()->GetScalars(), 0 );

	// assign tensors to output
	output->GetPointData()->SetScalars( newTensors );

	reader->Delete();
	newTensors->Delete();

	return 0;
}

#endif // FILE_nbfDTIReader