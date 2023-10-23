#pragma once

/** @file nbfBlitzReader.h
*	Blitz file reader. Part of IO.
*/

#include <io/nbfFileReader.h>

/** BLITZ file reader.
	Read dumped Blitz files (nbfBlitzWriter).

	@see nbfBlitzWriter
*/
class nbfBlitzReader : public nbfFileReader
{

public:

	/// Read data into 1D array. Return 1 if error, 0 otherwise.
	template< class Pixel >
	int read( Array< Pixel, 1 > & );

	/// Read data into 2D array. Return 1 if error, 0 otherwise.
	template< class Pixel >
	int read( Array< Pixel, 2 > & );

	/// Read data into 3D array. Return 1 if error, 0 otherwise.
	template< class Pixel >
	int read( Array< Pixel, 3 > & );

};


// read blitz data from file. Return 1 if error ocurred.

#define BZ_DECLARE_BLITZ_FILE_READER(dimension)				\
template< class Pixel >										\
int nbfBlitzReader :: read( Array< Pixel, dimension > & A )	\
{															\
  ifstream inFile( this->fileName );						\
  if ( inFile.bad() ){										\
    return 1;												\
  }															\
  inFile >> A;												\
  inFile.close();											\
  return 0;													\
}															


BZ_DECLARE_BLITZ_FILE_READER(1)
BZ_DECLARE_BLITZ_FILE_READER(2)
BZ_DECLARE_BLITZ_FILE_READER(3)