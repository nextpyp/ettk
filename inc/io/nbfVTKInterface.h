#pragma once

/** @file nbfVTKInterface.h
	VTK format conversion. Part of IO suite.
*/

#include <vtkImageData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkShortArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkImageCast.h>
#include <vtkImagePermute.h>

/** Image format conversion between Blitz++ and VTK libraries.
    
	Converts between Blitz++ Arrays< Pixel, Dim > and VTK's vtkImageData.
	Conversion is done *without copying* whenever possible. 
	
	If coverting from Blitz to VTK and Array storage is not contiguous
	then data has to be copied over.
		
	Also, when types are different (this only applies when converting VTK to Blitz) data has
	to be cast and therefore *copied* into the Blitz array.

	@warning VTK and Blitz use different storage order. Be warned that Blitz arrays
	produced with this class are stored in ColumnMajor order (opposed to the default RowMajor order). 
	Although this is done to avoid copying data, it prevents these arrays 
	to be used in methods that handle pointers to Array's data directly.
*/

class nbfVTKInterface
{

public:

	nbfVTKInterface(){
	}

	/** Convert vtkImageData to Array< Pixel, 2 >.
		If Array type is same as in vtkImageData, data is converted *without* copying.
		If types are different, VTK data is cast and copied into Blitz.
	*/
	static int vtkToBlitz( vtkImageData *, Array< double, 2 > & );
	static int vtkToBlitz( vtkImageData *, Array< float, 2 > & );
	static int vtkToBlitz( vtkImageData *, Array< short, 2 > & );
	static int vtkToBlitz( vtkImageData *, Array< int, 2 > & );
	static int vtkToBlitz( vtkImageData *, Array< unsigned char, 2 > & );


	/** Convert Array< Pixel, 2 > to vtkImageData.
		Output vtkImageData type is set to Blitz pixel type.
		Conversion can be done without copying *only* when Blitz data is contiguous and column-major.
		Otherwise, data has to be *copied* into vtkImageData.
	
		@caveat Before calling this method, you may generate a copy of the Array with contiguous storage
		so conversion can be done without copying.
	*/
	static int blitzToVtk( Array< double, 2 > &, vtkImageData * );
	static int blitzToVtk( Array< float, 2 > &, vtkImageData * );
	static int blitzToVtk( Array< short, 2 > &, vtkImageData * );
	static int blitzToVtk( Array< int, 2 > &, vtkImageData * );
	static int blitzToVtk( Array< unsigned char, 2 > &, vtkImageData * );

	/** Convert vtkImageData to Array< Pixel, 3 >.
		If Array type is same as in vtkImageData, data is converted *without* copying.
		If types are different, VTK data is cast and copied into Blitz.
	*/
	static int vtkToBlitz( vtkImageData *, Array< double, 3 > & );
	static int vtkToBlitz( vtkImageData *, Array< float, 3 > & );
	static int vtkToBlitz( vtkImageData *, Array< short, 3 > & );
	static int vtkToBlitz( vtkImageData *, Array< int, 3 > & );
	static int vtkToBlitz( vtkImageData *, Array< unsigned char, 3 > & );

	static int vtkToBlitz( vtkImageData *, Array< float, 3 > &, Array< float, 3 > & );
	static int vtkToBlitz( vtkImageData *, Array< double, 3 > &, Array< double, 3 > & );

	static int vtkToBlitzReference( vtkImageData *, Array< float, 3 > &, Array< float, 3 > & );
	static int vtkToBlitzReference( vtkImageData *, Array< double, 3 > &, Array< double, 3 > & );

	static int vtkToBlitz( vtkImageData *, Array< float, 2 > &, Array< float, 2 > & );
	static int vtkToBlitz( vtkImageData *, Array< double, 2 > &, Array< double, 2 > & );

	static int vtkToBlitzReference( vtkImageData *, Array< float, 2 > &, Array< float, 2 > & );
	static int vtkToBlitzReference( vtkImageData *, Array< double, 2 > &, Array< double, 2 > & );

	/** Create reference of vtkImageData to Array< Pixel, 2 >.
		Only if Array type is same as in vtkImageData, operation is succesful.
		If types are different, an assert message is generated.
	*/
	static int vtkToBlitzReference( vtkImageData *, Array< double, 2 > & );
	static int vtkToBlitzReference( vtkImageData *, Array< float, 2 > & );
	static int vtkToBlitzReference( vtkImageData *, Array< short, 2 > & );
	static int vtkToBlitzReference( vtkImageData *, Array< int, 2 > & );
	static int vtkToBlitzReference( vtkImageData *, Array< unsigned char, 2 > & );

	/** Create reference of vtkImageData to Array< Pixel, 3 >.
		Only if Array type is same as in vtkImageData, operation is succesful.
		If types are different, an assert message is generated.
	*/
	static int vtkToBlitzReference( vtkImageData *, Array< double, 3 > & );
	static int vtkToBlitzReference( vtkImageData *, Array< float, 3 > & );
	static int vtkToBlitzReference( vtkImageData *, Array< short, 3 > & );
	static int vtkToBlitzReference( vtkImageData *, Array< int, 3 > & );
	static int vtkToBlitzReference( vtkImageData *, Array< unsigned char, 3 > & );


	/** Convert Array< Pixel, 3 > to vtkImageData.
		Output vtkImageData type is set to Blitz pixel type.
		Conversion can be done without copying *only* when Blitz data is contiguous and column-major.
		Otherwise, data has to be *copied* into vtkImageData.
	
		@caveat Before calling this method, you may generate a copy of the Array with contiguous storage
		so conversion can be done without copying.
	*/
	static int blitzToVtk( Array< double, 3 > &, vtkImageData * );
	static int blitzToVtk( Array< float, 3 > &, vtkImageData * );
	static int blitzToVtk( Array< short, 3 > &, vtkImageData * );
	static int blitzToVtk( Array< int, 3 > &, vtkImageData * );
	static int blitzToVtk( Array< unsigned char, 3 > &, vtkImageData * );

	static int blitzToVtk( Array< complex< float >, 3 > &, vtkImageData * );
	static int blitzToVtk( Array< complex< double >, 3 > &, vtkImageData * );

protected:

	/** Copy data from Blitz array to VTK format. 
		Uses brute-force (3 for loops) for handling indexes correctly.
	*/
	template < class Pixel, const int Dim >
	static void setScalars( Array< Pixel, Dim > &, vtkDataArray * );

	/** Copy data from Blitz array to VTK format. 
		Uses brute-force (3 for loops) for handling indexes correctly.
	*/
	template < class Pixel, const int Dim >
	static void setScalars( Array< complex< Pixel >, Dim > &, vtkDataArray * );

	/** Copy data from VTK format to Blitz array. 
		Uses brute-force (3 for loops) for handling indexes correctly.
		If VTK data is multicomponent, only first component is considered.
	*/
	template < class Pixel, const int Dim >
	static void setScalars( vtkDataArray *, Array< Pixel, Dim > &, int * );

};