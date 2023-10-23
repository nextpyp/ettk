#pragma once

/** @file nbfBlitzWriter.h
*	Blitz file writer. Part of IO.
*/

#include <io/nbfFileWriter.h>

/** BLITZ file writer.
	Dump Blitz arrays into file.

	@see nbfBlitzReader
*/
class nbfBlitzWriter : public nbfFileWriter
{

public:

	/// Write 1D blitz data to file. Return 1 if error ocurred.
	template< class Pixel >
	int write( Array< Pixel, 1 > & );

	/// Write 2D blitz data to file. Return 1 if error ocurred.
	template< class Pixel >
	int write( Array< Pixel, 2 > & );
	
	/// Write 2D blitz data to file. Return 1 if error ocurred.
	template< class Pixel >
	static int write( char *, Array< Pixel, 2 > & );

	/// Write 3D blitz data to file. Return 1 if error ocurred.
	template< class Pixel >
	int write( Array< Pixel, 3 > & );
};


#define BZ_DECLARE_BLITZ_FILE_WRITER(dimension)				\
template< class Pixel >										\
int nbfBlitzWriter :: write( Array< Pixel, dimension > & A )\
{															\
  ofstream outFile( this->fileName );						\
  if ( outFile.bad() ){										\
    return 1;												\
  }															\
  outFile << A << endl;										\
  outFile.close();											\
  return 0;													\
}

BZ_DECLARE_BLITZ_FILE_WRITER(1)
BZ_DECLARE_BLITZ_FILE_WRITER(2)
BZ_DECLARE_BLITZ_FILE_WRITER(3)

template< class Pixel >
int nbfBlitzWriter :: write( char * file, Array< Pixel, 2 > & A ){
	ofstream outFile( file );
	if ( outFile.bad() ){
		return 1;
	}
	for ( int i = 0; i < A.rows(); i++ ){
		for ( int j = 0; j < A.cols(); j++ ){
            outFile << A(i,j);
			if ( j < A.cols() - 1 ){
				outFile << "\t";
			}
		}
		outFile << endl;
	}
	outFile.close();
	return 0;
}