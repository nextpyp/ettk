#ifndef FILE_nbfDTIFieldReader
#define FILE_nbfDTIFieldReader

#include <vtkImageReader.h>

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkArrayCalculator.h>

#include <blitz/array.h>
using namespace blitz;

#include <io/nbfVTKInterface.h>

class nbfDTIFieldReader
{

public:

	// Read data into array. Return 1 if error, 0 otherwise.
	template< class Pixel >
	int read( Array< Pixel, 2 > & );

	template< class Pixel >
	int read( Array< Pixel, 3 > & );

	template< class Pixel >
	int read( Array< Pixel, 3 > &, Array< Pixel, 3 > &, Array< Pixel, 3 > & );

	// read into 3D vtk volume
	int read( vtkImageData * );

	nbfDTIFieldReader(){};

};

template< class Pixel >
int nbfDTIFieldReader :: read( Array< Pixel, 2 > & A ){
	return 0;
}

template< class Pixel >
int nbfDTIFieldReader :: read( Array< Pixel, 3 > & A ){
	vtkImageData * output = vtkImageData::New();
	this->read( output );
	nbfVTKInterface::vtkToBlitz( output, A );
	output->Delete();
	return 0;
}

template< class Pixel >
int nbfDTIFieldReader :: read( Array< Pixel, 3 > & Ax, Array< Pixel, 3 > & Ay, Array< Pixel, 3 > & Az )
{
	// initialize reader
	vtkImageReader * reader;
	reader = vtkImageReader::New();
	reader->SetDataByteOrderToBigEndian();
	reader->SetDataExtent(0, 125, 0, 159, 0, 46);
	reader->SetDataSpacing(1.25, 1.25, 2.5);
	reader->SetDataScalarTypeToFloat();
	reader->SetFileDimensionality(3);
	
	// first read coordinates + FA image
	reader->SetFileName("C:/home/DTI/data/ISMRM/FA.img");
	reader->Update();

	// store full data set (3D coordinates + FA + tensor)
	//output->DeepCopy( reader->GetOutput() );

	// get FA image
	Array< Pixel, 3 > A;

	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), A );

	// normalize to [0,1]
	A = A / max(A);

	// get Vx
	reader->SetFileName("C:/home/DTI/data/ISMRM/evec_max_compx.img");
	reader->Update();
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), Ax );

	// get Vy
	reader->SetFileName("C:/home/DTI/data/ISMRM/evec_max_compy.img");
	reader->Update();
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), Ay );

	// get Vz
	reader->SetFileName("C:/home/DTI/data/ISMRM/evec_max_compz.img");
	reader->Update();
	nbfVTKInterface::vtkToBlitz( reader->GetOutput(), Az );

	// scale all three components by FA image
	Ax = Ax * A;
	Ay = Ay * A;
	Az = Az * A;

	reader->Delete();
	return 0;
}

int nbfDTIFieldReader :: read( vtkImageData * output ){

	// initialize reader
	vtkImageReader * reader;
	reader = vtkImageReader::New();
	reader->SetDataByteOrderToBigEndian();
	reader->SetDataExtent(0, 125, 0, 159, 0, 46);
	reader->SetDataSpacing(-1.25, 1.25, 2.5);
	reader->SetDataScalarTypeToFloat();
	reader->SetFileDimensionality(3);
	
	// first read coordinates + FA image
	reader->SetFileName("C:/home/DTI/data/ISMRM/FA.img");
	reader->Update();

	// store full data set (3D coordinates + FA + tensor)
	output->DeepCopy( reader->GetOutput() );

	// store tensor values (3x3)
	vtkFloatArray * newTensors = vtkFloatArray::New();
	newTensors->SetNumberOfComponents(3);
	newTensors->SetNumberOfTuples( output->GetNumberOfPoints() );

	// get Vx
	reader->SetFileName("C:/home/DTI/data/ISMRM/evec_max_compx.img");
	reader->Update();
	newTensors->CopyComponent( 0, reader->GetOutput()->GetPointData()->GetScalars(), 0 );

	// get Vy
	reader->SetFileName("C:/home/DTI/data/ISMRM/evec_max_compy.img");
	reader->Update();
	newTensors->CopyComponent( 1, reader->GetOutput()->GetPointData()->GetScalars(), 0 );

	// get Vz
	reader->SetFileName("C:/home/DTI/data/ISMRM/evec_max_compz.img");
	reader->Update();
	newTensors->CopyComponent( 2, reader->GetOutput()->GetPointData()->GetScalars(), 0 );

	// assign tensors to output
	output->GetPointData()->SetVectors( newTensors );

	reader->Delete();
	newTensors->Delete();

	return 0;
}

#endif // FILE_nbfDTIFieldReader