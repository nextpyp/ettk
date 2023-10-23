#pragma once

/** @file nbfImageWriter.h
	Write images and .vtk files. Part of IO.
*/

#include <io/nbfFileWriter.h>
#include <io/nbfVTKInterface.h>

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCharArray.h>
#include <vtkShortArray.h>
#include <vtkImageCast.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageWriter.h>
#include <vtkBMPWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkPNMWriter.h>
#include <vtkPNGWriter.h>

#include <vtkImagePermute.h>

#include <vtkCellData.h>

/** Write Blitz arrays into images (2D) and vtkStructuredPoints (3D).

	@see nbfImageReader
*/
class nbfImageWriter : public nbfFileWriter
{

public:

	/** Save 2D Blitz array data into image file.
		As image data must be cast to _unsigned char_ before writing, all image 
		values are assumed to be in the [0,255] range to avoid scaling problems.
		File format is deducted from the filename extension. 
		Supported formats include: bmp, tif, png, jpg and pnm.

		@caveats Only gray scale images are currently supported.
		
		@todo Extend to color images.
	*/
	template< class Pixel >
	int write( Array< Pixel, 2 > & );

	/** Save 3D Blitz array data into vtk format (vtkStructuredPoints).
	*/
	template< class Pixel >
	int write( Array< Pixel, 3 > & );
};

template< class Pixel >
int nbfImageWriter :: write( Array< Pixel, 2 > & A )
{
  vtkImageData * image = vtkImageData::New();

  nbfVTKInterface::blitzToVtk( A, image );

  vtkImageCast * cast = vtkImageCast::New();
  cast->SetInput( image );
  cast->SetOutputScalarTypeToUnsignedChar();

  vtkImageWriter * writer;

  int lenght = strlen( this->fileName );

  if ( ( this->fileName[ lenght - 4 ] == '.' ) &&
       ( this->fileName[ lenght - 3 ] == 'b' ) &&
       ( this->fileName[ lenght - 2 ] == 'm' ) &&
       ( this->fileName[ lenght - 1 ] == 'p' ) )
    {
      writer = vtkBMPWriter::New();
    }
  else if ( ( this->fileName[ lenght - 4 ] == '.' ) &&
	        ( this->fileName[ lenght - 3 ] == 't' ) &&
	        ( this->fileName[ lenght - 2 ] == 'i' ) &&
	        ( this->fileName[ lenght - 1 ] == 'f' ) )
    {
      writer = vtkTIFFWriter::New();
    }
  else if ( ( this->fileName[ lenght - 4 ] == '.' ) &&
	        ( this->fileName[ lenght - 3 ] == 'p' ) &&
	        ( this->fileName[ lenght - 2 ] == 'n' ) &&
	        ( this->fileName[ lenght - 1 ] == 'm' ) )
    {
      writer = vtkPNMWriter::New();
    }
  else if ( ( this->fileName[ lenght - 4 ] == '.' ) &&
	        ( this->fileName[ lenght - 3 ] == 'j' ) &&
	        ( this->fileName[ lenght - 2 ] == 'p' ) &&
	        ( this->fileName[ lenght - 1 ] == 'g' ) )
    {
      writer = vtkJPEGWriter::New();
    }
  else if ( ( this->fileName[ lenght - 4 ] == '.' ) &&
	        ( this->fileName[ lenght - 3 ] == 'p' ) &&
	        ( this->fileName[ lenght - 2 ] == 'n' ) &&
	        ( this->fileName[ lenght - 1 ] == 'g' ) )
    {
      writer = vtkJPEGWriter::New();
    }
  else{
    return 1;
  }

  writer->SetFileName( this->fileName );
  writer->SetInput( cast->GetOutput() );
  writer->Write();
  
  cast->Delete();
  writer->Delete();
  image->Delete();
  
  return 0;
}

template< class Pixel >
int nbfImageWriter :: write( Array< Pixel, 3 > & A )
{
  vtkImageData * image = vtkImageData::New();

  nbfVTKInterface::blitzToVtk( A, image );

  vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();

  writer->SetFileName( this->fileName );
  writer->SetInput( image );
  writer->SetFileTypeToBinary();
  writer->Write();
  
  writer->Delete();
  image->Delete();
  
  return 0;
}