/** @file nbfVTKInterface.cpp
	VTK format conversion. Part of IO suite.
*/

#include <blitz/array.h>
using namespace blitz;

#include <io/nbfVTKInterface.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkShortArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkImageCast.h>
#include <vtkImagePermute.h>

#define NBF_DECLARE_VTK_TO_BLITZ_2(type1,type2,type3)									\
int nbfVTKInterface :: vtkToBlitz( vtkImageData * input, Array< type1, 2 > & output ){	\
int dimensions[3];																		\
input->GetDimensions( dimensions );														\
nbfVTKInterface::setScalars( input->GetPointData()->GetScalars(), output, dimensions ); \
return 1;																				\
}

NBF_DECLARE_VTK_TO_BLITZ_2(short,Short,SHORT)
NBF_DECLARE_VTK_TO_BLITZ_2(float,Float,FLOAT)
NBF_DECLARE_VTK_TO_BLITZ_2(double,Double,DOUBLE)
NBF_DECLARE_VTK_TO_BLITZ_2(int,Int,INT)
NBF_DECLARE_VTK_TO_BLITZ_2(unsigned char,UnsignedChar,UNSIGNED_CHAR)

#define NBF_DECLARE_BLITZ_TO_VTK_2(type1,type2)											\
int nbfVTKInterface :: blitzToVtk( Array< type1, 2 > & input, vtkImageData * output ){	\
output->SetDimensions( input.rows(), input.cols(), 1 );									\
output->SetScalarTypeTo ## type2 ();													\
output->AllocateScalars();																\
nbfVTKInterface::setScalars( input, output->GetPointData()->GetScalars() ); 			\
return 0;																				\
}

//#define NBF_DECLARE_BLITZ_TO_VTK_2(type1,type2)											\
//int nbfVTKInterface :: blitzToVtk( Array< type1, 2 > & input, vtkImageData * output ){	\
//output->SetDimensions( input.rows(), input.cols(), 1 );									\
//output->SetScalarTypeTo ## type2 ();													\
//output->AllocateScalars();																\
//vtk ## type2 ## Array * scalars = vtk ## type2 ## Array::New();							\
//nbfVTKInterface::setScalars( input, scalars);											\
//output->GetPointData()->SetScalars( scalars );											\
//return 0;																				\
//}

NBF_DECLARE_BLITZ_TO_VTK_2(double,Double)
NBF_DECLARE_BLITZ_TO_VTK_2(float,Float)
NBF_DECLARE_BLITZ_TO_VTK_2(short,Short)
NBF_DECLARE_BLITZ_TO_VTK_2(int,Int)
NBF_DECLARE_BLITZ_TO_VTK_2(unsigned char,UnsignedChar)


#define NBF_DECLARE_VTK_TO_BLITZ_3(type1,type2,type3)																								\
int nbfVTKInterface :: vtkToBlitz( vtkImageData * input, Array< type1, 3 > & output ){																\
int dimensions[3];																																	\
input->GetDimensions( dimensions );																													\
if ( dimensions[0] + dimensions[1] + dimensions[2] > 0 ){																							\
vtk ## type2 ## Array * pointerToData;																												\
if ( input->GetNumberOfScalarComponents() == 1 ){																									\
	vtkImageCast * cast = vtkImageCast::New();																										\
	if ( input->GetPointData()->GetScalars()->GetDataType() == VTK_ ## type3 ){																		\
	pointerToData = vtk ## type2 ## Array::SafeDownCast( input->GetPointData()->GetScalars() );														\
	}																																				\
	else{																																			\
	cast->SetOutputScalarTypeTo ## type2 ();																										\
	cast->SetInput( input );																														\
	cast->Update();																																	\
	pointerToData = vtk ## type2 ## Array::SafeDownCast( cast->GetOutput()->GetPointData()->GetScalars() );											\
	}																																				\
	Array< type1, 3 > A( pointerToData->GetPointer(0), shape(dimensions[0],dimensions[1],dimensions[2]), neverDeleteData, ColumnMajorArray<3>() );	\
	output.resize( A.shape() );																														\
	output = A;																																		\
	cast->Delete();																																	\
	return 1;																																		\
}																																					\
else{																																				\
int dimensions[3];																																	\
input->GetDimensions(dimensions);																													\
nbfVTKInterface::setScalars( input->GetPointData()->GetScalars(), output, dimensions );																\
}																																					\
}																																					\
return 1;																																			\
}

NBF_DECLARE_VTK_TO_BLITZ_3(short,Short,SHORT)
NBF_DECLARE_VTK_TO_BLITZ_3(float,Float,FLOAT)
NBF_DECLARE_VTK_TO_BLITZ_3(double,Double,DOUBLE)
NBF_DECLARE_VTK_TO_BLITZ_3(int,Int,INT)
NBF_DECLARE_VTK_TO_BLITZ_3(unsigned char,UnsignedChar,UNSIGNED_CHAR)

#define NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_2(type1,type2,type3)\
int nbfVTKInterface :: vtkToBlitzReference( vtkImageData * input, Array< type1, 2 > & output ){\
if ( ( input->GetPointData()->GetScalars()->GetDataType() != VTK_ ## type3 ) || ( input->GetNumberOfScalarComponents() != 1 ) ){\
	assert("Incompatible data type.");\
	cout << "ERROR - cannot convert VTK data. Incompatible data type." << endl;\
	return 0;\
}\
else{\
int dimensions[3];\
input->GetDimensions( dimensions );\
TinyVector< int, 2 > strides( 1, dimensions[0] );\
vtk ## type2 ## Array * pointerToData = vtk ## type2 ## Array::SafeDownCast( input->GetPointData()->GetScalars() );\
Array< type1, 2 > A( pointerToData->GetPointer(0), shape(dimensions[0],dimensions[1]), strides, neverDeleteData, ColumnMajorArray<2>() );\
output.reference( A );\
return 0;\
}\
}

NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_2(short,Short,SHORT)
NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_2(float,Float,FLOAT)
NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_2(double,Double,DOUBLE)
NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_2(int,Int,INT)
NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_2(unsigned char,UnsignedChar,UNSIGNED_CHAR)


#define NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_3(type1,type2,type3)\
int nbfVTKInterface :: vtkToBlitzReference( vtkImageData * input, Array< type1, 3 > & output ){\
if ( ( input->GetScalarType() != VTK_ ## type3 ) || ( input->GetNumberOfScalarComponents() != 1 ) ){\
	cerr << "Cannot convert VTK to Blitz as reference. Input type: " << input->GetPointData()->GetScalars()->GetDataType() << ", output type: " << VTK_ ## type3 << ", Number of scalar components: " << input->GetNumberOfScalarComponents() << endl;\
	cerr << "ERROR" << __FILE__ << ", " << __LINE__ << endl;\
	assert(0);\
	return 0;\
}\
else{\
int dimensions[3];\
input->GetDimensions( dimensions );\
TinyVector< int, 3 > strides( 1, dimensions[0], dimensions[0] * dimensions[1] );\
Array< type1, 3 > A( (type1*)input->GetScalarPointer(), shape(dimensions[0],dimensions[1],dimensions[2]), strides, neverDeleteData, ColumnMajorArray<3>() );\
output.reference( A );\
return 0;\
}\
}

NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_3(short,Short,SHORT)
NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_3(float,Float,FLOAT)
NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_3(double,Double,DOUBLE)
NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_3(int,Int,INT)
NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_3(unsigned char,UnsignedChar,UNSIGNED_CHAR)


#define NBF_DECLARE_BLITZ_TO_VTK_3(type1,type2)														\
int nbfVTKInterface :: blitzToVtk( Array< type1, 3 > & input, vtkImageData * output ){				\
output->SetDimensions( input.rows(), input.cols(), input.depth() );									\
output->SetScalarTypeTo ## type2 ();																\
vtk ## type2 ## Array * scalars = vtk ## type2 ## Array::New();										\
if ( false && input.isStorageContiguous() ){																	\
	cout << "continuous storage" << endl;												\
	scalars->SetArray( input.data(), input.size(), 1 );												\
	output->GetPointData()->SetScalars( scalars );													\
	if ( input.isMinorRank(firstDim) ){																\
		output->SetDimensions( input.depth(), input.cols(), input.rows() );							\
		output->GetPointData()->SetScalars( scalars );												\
		vtkImagePermute * permute = vtkImagePermute::New();											\
		permute->SetFilteredAxes(2,1,0);															\
		permute->SetInput( output );																\
		permute->Update();																			\
		output->DeepCopy( permute->GetOutput() );													\
		output->SetUpdateExtent( permute->GetOutput()->GetUpdateExtent() );							\
		permute->Delete();																			\
	}																								\
}																									\
else{																								\
nbfVTKInterface::setScalars( input, scalars );														\
output->GetPointData()->SetScalars( scalars );														\
}																									\
scalars->Delete();																					\
return 0;																							\
}

NBF_DECLARE_BLITZ_TO_VTK_3(double,Double)
NBF_DECLARE_BLITZ_TO_VTK_3(float,Float)
NBF_DECLARE_BLITZ_TO_VTK_3(short,Short)
NBF_DECLARE_BLITZ_TO_VTK_3(int,Int)
NBF_DECLARE_BLITZ_TO_VTK_3(unsigned char,UnsignedChar)


template < class Pixel, const int Dim >
void nbfVTKInterface :: setScalars( Array< Pixel, Dim > & A, vtkDataArray * scalars )
{
	if ( Dim == 3 ){
		for ( int k = 0; k < A.depth(); k++ )
		{
			int kOffset = k *  A.cols() * A.rows();
			for ( int j = 0; j <  A.cols(); j++ )
			{
				int jOffset = j * A.rows();
				for ( int i = 0; i < A.rows() ; i++ )
				{
					int offset = kOffset + jOffset + i;
					scalars->InsertTuple1( offset, A( i, j, k ) );
				}
			}
		}
	}
	else if ( Dim == 2 ){
		for ( int j = 0; j <  A.cols(); j++ )
		{
			int jOffset = j * A.rows();
			for ( int i = 0; i < A.rows() ; i++ )
			{
				int offset = jOffset + i;
				scalars->InsertTuple1( offset, A( i, j ) );
			}
		}
	}
}

template < class Pixel, const int Dim >
void nbfVTKInterface :: setScalars( Array< complex< Pixel >, Dim > & A, vtkDataArray * scalars )
{
	scalars->SetNumberOfComponents(2);
	if ( Dim == 3 ){
		for ( int k = 0; k < A.depth(); k++ )
		{
			int kOffset = k *  A.cols() * A.rows();
			for ( int j = 0; j <  A.cols(); j++ )
			{
				int jOffset = j * A.rows();
				for ( int i = 0; i < A.rows() ; i++ )
				{
					int offset = kOffset + jOffset + i;
					complex< Pixel > c = A(i,j,k);
					scalars->InsertTuple2( offset, real(c), imag(c) );
				}
			}
		}
	}
	else if ( Dim == 2 ){
		for ( int j = 0; j <  A.cols(); j++ )
		{
			int jOffset = j * A.rows();
			for ( int i = 0; i < A.rows() ; i++ )
			{
				int offset = jOffset + i;
				complex< Pixel > c = A(i,j);
				scalars->InsertTuple2( offset, real(c), imag(c) );
			}
		}
	}
}

template < class Pixel, const int Dim >
void nbfVTKInterface :: setScalars( vtkDataArray * scalars, Array< Pixel, Dim > & A, int * dimensions )
{
	if ( Dim == 3 ){
		A.resize( dimensions[0], dimensions[1], dimensions[2] );
		for ( int k = 0; k < dimensions[2]; k++ )
		{
			int kOffset = k *  dimensions[1] * dimensions[0];
			for ( int j = 0; j <  dimensions[1]; j++ )
			{
				int jOffset = j * dimensions[0];
				for ( int i = 0; i < dimensions[0] ; i++ )
				{
					int offset = kOffset + jOffset + i;
					A( i, j, k ) = (Pixel)scalars->GetComponent( offset, 0 );
				}
			}
		}
	}
	else if ( Dim == 2 ){
		A.resize( dimensions[0], dimensions[1] );
		for ( int j = 0; j <  dimensions[1]; j++ )
		{
			int jOffset = j * dimensions[0];
			for ( int i = 0; i < dimensions[0] ; i++ )
			{
				int offset = jOffset + i;
				A( i, j ) = (Pixel)scalars->GetComponent( offset, 0 );
			}
		}
	}
}

#define NBF_DECLARE_BLITZ_TO_VTK_COMPLEX_2(type1,type2,type3)\
int nbfVTKInterface :: vtkToBlitz( vtkImageData * input, Array< type1, 2 > & real, Array< type1, 2 > & imag ){\
int dimensions[3];\
input->GetDimensions( dimensions );\
vtk ## type2 ## Array * pointerToData;\
int comp = input->GetNumberOfScalarComponents();\
if ( comp == 2 ){\
	vtkImageCast * cast = vtkImageCast::New();																										\
	if ( input->GetPointData()->GetScalars()->GetDataType() == VTK_ ## type3 ){																		\
	pointerToData = vtk ## type2 ## Array::SafeDownCast( input->GetPointData()->GetScalars() );														\
	}																																				\
	else{																																			\
	cast->SetOutputScalarTypeTo ## type2 ();																										\
	cast->SetInput( input );																														\
	cast->Update();																																	\
	pointerToData = vtk ## type2 ## Array::SafeDownCast( cast->GetOutput()->GetPointData()->GetScalars() );											\
	}																																				\
	TinyVector< int, 2 > strides( comp, comp * dimensions[0] );\
	Array< type1, 2 > A( pointerToData->GetPointer(0), shape(dimensions[0],dimensions[1]), strides, neverDeleteData, ColumnMajorArray<2>() );\
	real.resize( A.shape() );\
	real = A;\
	Array< type1, 2 > B( pointerToData->GetPointer(0) + 1, shape(dimensions[0],dimensions[1]), strides, neverDeleteData, ColumnMajorArray<2>() );\
	imag.resize( B.shape() );\
	imag = B;\
	cast->Delete();\
}\
else{\
	cerr << "ERROR converting to vtk. Data does not have 2 components." << endl;\
	return 0;\
}\
return 1;\
}

NBF_DECLARE_BLITZ_TO_VTK_COMPLEX_2(float,Float,FLOAT)
NBF_DECLARE_BLITZ_TO_VTK_COMPLEX_2(double,Double,DOUBLE)

#define NBF_DECLARE_BLITZ_TO_VTK_REFERENCE_COMPLEX_2(type1,type2,type3)\
int nbfVTKInterface :: vtkToBlitzReference( vtkImageData * input, Array< type1, 2 > & real, Array< type1, 2 > & imag ){\
int dimensions[3];\
input->GetDimensions( dimensions );\
vtk ## type2 ## Array * pointerToData;\
int comp = input->GetNumberOfScalarComponents();\
if ( comp == 2 ){\
	if ( input->GetPointData()->GetScalars()->GetDataType() == VTK_ ## type3 ){\
		TinyVector< int, 2 > strides( comp, comp * dimensions[0] );\
		pointerToData = vtk ## type2 ## Array::SafeDownCast( input->GetPointData()->GetScalars() );\
		Array< type1, 2 > A( pointerToData->GetPointer(0), shape(dimensions[0],dimensions[1]), strides, neverDeleteData, ColumnMajorArray<2>() );\
		real.reference( A );\
		Array< type1, 2 > B( pointerToData->GetPointer(0) + 1, shape(dimensions[0],dimensions[1]), strides, neverDeleteData, ColumnMajorArray<2>() );\
		imag.reference( B );\
	}\
	else{\
	    cerr << "ERROR converting to vtk. Data type does not match." << endl;\
		return 0;\
	}\
}\
else{\
	cerr << "ERROR converting to vtk. Data does not have 2 components." << endl;\
	return 0;\
}\
return 1;\
}

NBF_DECLARE_BLITZ_TO_VTK_REFERENCE_COMPLEX_2(float,Float,FLOAT)
NBF_DECLARE_BLITZ_TO_VTK_REFERENCE_COMPLEX_2(double,Double,DOUBLE)

#define NBF_DECLARE_VTK_TO_BLITZ_COMPLEX_3(type1,type2,type3)\
int nbfVTKInterface :: vtkToBlitz( vtkImageData * input, Array< type1, 3 > & real, Array< type1, 3 > & imag ){\
int dimensions[3];\
input->GetDimensions( dimensions );\
vtk ## type2 ## Array * pointerToData;\
int comp = input->GetNumberOfScalarComponents();\
if ( comp == 2 ){\
	if ( input->GetPointData()->GetScalars()->GetDataType() == VTK_ ## type3 ){\
		TinyVector< int, 3 > strides( comp, comp * dimensions[0], comp * dimensions[0] * dimensions[1] );\
		pointerToData = vtk ## type2 ## Array::SafeDownCast( input->GetPointData()->GetScalars() );\
		Array< type1, 3 > A( pointerToData->GetPointer(0), shape(dimensions[0],dimensions[1],dimensions[2]), strides, neverDeleteData, ColumnMajorArray<3>() );\
		real.reference( A );\
		Array< type1, 3 > B( pointerToData->GetPointer(0) + 1, shape(dimensions[0],dimensions[1],dimensions[2]), strides, neverDeleteData, ColumnMajorArray<3>() );\
		imag.reference( B );\
	}\
	else{\
		if ( input->GetPointData()->GetScalars()->GetDataType() == VTK_DOUBLE ){\
			TinyVector< int, 3 > strides( comp, comp * dimensions[0], comp * dimensions[0] * dimensions[1] );\
			vtkDoubleArray * pointer = vtkDoubleArray::SafeDownCast( input->GetPointData()->GetScalars() );\
			Array< double, 3 > A( pointer->GetPointer(0), shape(dimensions[0],dimensions[1],dimensions[2]), strides, neverDeleteData, ColumnMajorArray<3>() );\
			real.resize( A.shape() );\
			real = cast< type1 >(A);\
			Array< double, 3 > B( pointer->GetPointer(0) + 1, shape(dimensions[0],dimensions[1],dimensions[2]), strides, neverDeleteData, ColumnMajorArray<3>() );\
			imag.resize( B.shape() );\
			imag = cast< double >(B);\
		}\
		else{\
		    cerr << "ERROR converting to vtk." << endl;\
			assert(false);\
		}\
	}\
}\
else{\
	cerr << "ERROR converting to vtk. Data does not have 2 components." << endl;\
	assert(false);\
}\
return 1;\
}

NBF_DECLARE_VTK_TO_BLITZ_COMPLEX_3(float,Float,FLOAT)
NBF_DECLARE_VTK_TO_BLITZ_COMPLEX_3(double,Double,DOUBLE)

#define NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_COMPLEX_3(type1,type2,type3)\
int nbfVTKInterface :: vtkToBlitzReference( vtkImageData * input, Array< type1, 3 > & real, Array< type1, 3 > & imag ){\
int dimensions[3];\
input->GetDimensions( dimensions );\
vtk ## type2 ## Array * pointerToData;\
int comp = input->GetNumberOfScalarComponents();\
if ( comp == 2 ){\
	if ( input->GetPointData()->GetScalars()->GetDataType() == VTK_ ## type3 ){\
		TinyVector< int, 3 > strides( comp, comp * dimensions[0], comp * dimensions[0] * dimensions[1] );\
		pointerToData = vtk ## type2 ## Array::SafeDownCast( input->GetPointData()->GetScalars() );\
		Array< type1, 3 > A( pointerToData->GetPointer(0), shape(dimensions[0],dimensions[1],dimensions[2]), strides, neverDeleteData, ColumnMajorArray<3>() );\
		real.reference( A );\
		Array< type1, 3 > B( pointerToData->GetPointer(0) + 1, shape(dimensions[0],dimensions[1],dimensions[2]), strides, neverDeleteData, ColumnMajorArray<3>() );\
		imag.reference( B );\
	}\
	else{\
	    cerr << "ERROR converting to vtk. Data type does not match." << __FILE__ << endl;\
		return 0;\
	}\
}\
else{\
	cerr << "ERROR converting to vtk. Data does not have 2 components." << __FILE__ << endl;\
	assert(false);\
}\
return 1;\
}

NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_COMPLEX_3(float,Float,FLOAT)
NBF_DECLARE_VTK_TO_BLITZ_REFERENCE_COMPLEX_3(double,Double,DOUBLE)

#define NBF_DECLARE_BLITZ_TO_VTK_COMPLEX_3(type1,type2,type3)\
int nbfVTKInterface :: blitzToVtk( Array< complex< type1 >, 3 > & input, vtkImageData * output  ){	\
output->SetDimensions( input.rows(), input.cols(), input.depth() );									\
output->SetScalarTypeTo ## type2 ();																\
output->SetNumberOfScalarComponents(2);																\
vtk ## type2 ## Array * scalars = vtk ## type2 ## Array::New();										\
nbfVTKInterface::setScalars( input, scalars );														\
output->GetPointData()->SetScalars( scalars );														\
scalars->Delete();																					\
return 0;																							\
}

NBF_DECLARE_BLITZ_TO_VTK_COMPLEX_3(float,Float,FLOAT)
NBF_DECLARE_BLITZ_TO_VTK_COMPLEX_3(double,Double,DOUBLE)