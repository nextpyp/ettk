#pragma once

/** @file nbfImageWriter.h
*	Image file writer. Part of IO.
*/

#include <io/nbfFileReader.h>

#include <vtkImageReader2.h>
#include <vtkBMPReader.h>
#include <vtkPNMReader.h>
#include <vtkTIFFReader.h>
#include <vtkJPEGReader.h>
#include <vtkPNGReader.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>

/** Read images (2D) and vtkStructuredPoints (3D) into Blitz arrays.

	@see nbfImageReader
*/
class nbfImageReader : public nbfFileReader
{

public:

	/** Read image into 2D Blitz array.
		Supported formats include: bmp, tif, png, jpg and pnm.

		@caveats Only gray scale images are currently supported. Color images are converted to b&w.
		
		@todo Extend to color images.
	*/
	template< class Pixel >
	int read( Array< Pixel, 2 > & );

	/** Read vtkStructuredPoints file into 3D Blitz array.
	*/
	template< class Pixel >
	int read( Array< Pixel, 3 > & );

};

template< class Pixel >
int nbfImageReader :: read( Array< Pixel, 2 > & A )
{
  vtkImageReader2 * reader;
  int lenght = strlen( this->fileName );
  
  if ( ( this->fileName[ lenght - 4 ] == '.' ) &&
       ( this->fileName[ lenght - 3 ] == 'b' ) &&
       ( this->fileName[ lenght - 2 ] == 'm' ) &&
       ( this->fileName[ lenght - 1 ] == 'p' ) )
    {
      reader = vtkBMPReader::New();
    }
  else if ( ( this->fileName[ lenght - 4 ] == '.' ) &&
	        ( this->fileName[ lenght - 3 ] == 't' ) &&
	        ( this->fileName[ lenght - 2 ] == 'i' ) &&
	        ( this->fileName[ lenght - 1 ] == 'f' ) )
    {
      reader = vtkTIFFReader::New();
    }
  else if ( ( this->fileName[ lenght - 4 ] == '.' ) &&
	        ( this->fileName[ lenght - 3 ] == 'p' ) &&
	        ( this->fileName[ lenght - 2 ] == 'n' ) &&
	        ( this->fileName[ lenght - 1 ] == 'm' ) )
    {
      reader = vtkPNMReader::New();
    }
  else if ( ( this->fileName[ lenght - 4 ] == '.' ) &&
	        ( this->fileName[ lenght - 3 ] == 'j' ) &&
	        ( this->fileName[ lenght - 2 ] == 'p' ) &&
	        ( this->fileName[ lenght - 1 ] == 'g' ) )
    {
      reader = vtkJPEGReader::New();
    }
  else if ( ( this->fileName[ lenght - 4 ] == '.' ) &&
	        ( this->fileName[ lenght - 3 ] == 'p' ) &&
	        ( this->fileName[ lenght - 2 ] == 'n' ) &&
	        ( this->fileName[ lenght - 1 ] == 'g' ) )
    {
      reader = vtkPNGReader::New();
    }
  else{
    return 1;
  }
  
  reader->SetFileName( this->fileName );
  reader->Update();

  if ( !reader->GetOutput() || !reader->GetOutput()->GetNumberOfPoints() ){
	  return 1;
  }
  else{

	  int rows = reader->GetOutput()->GetDimensions()[0];
	  int cols = reader->GetOutput()->GetDimensions()[1];
	  A.resize( rows, cols );

	 typename  Array< Pixel, 2 > :: iterator iter = A.begin();

	  for ( int i = 0; i < rows * cols; i++ ){
		  Pixel value = 0;
		  if ( reader->GetOutput()->GetPointData()->GetScalars()->GetNumberOfComponents() == 3 ){
			  double readed[3];
			  reader->GetOutput()->GetPointData()->GetScalars()->GetTuple(i,readed);
			  value = 0.299 * readed[0] + 0.587 * readed[1] + 0.114 * readed[2];
		  }
		  else{
			  value = reader->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i);
		  }
		  A( (int)reader->GetOutput()->GetPoint(i)[0], (int)reader->GetOutput()->GetPoint(i)[1] ) = value;
	  }
  }
  return 0;			  // todo ok.!
}

template< class Pixel >
int nbfImageReader :: read( Array< Pixel, 3 > & A )
{
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName( this->fileName );
	if ( reader->IsFileStructuredPoints() ){
		reader->Update();
		nbfVTKInterface::vtkToBlitz( reader->GetOutput(), A );
		return 1;
	}
	return 0;
}