#ifndef FILE_nbfFileReader
#define FILE_nbfFileReader

#include <io/nbfFile.h>

class nbfFileReader : public nbfFile
{

public:

	// Read data into array. Return 1 if error, 0 otherwise.
	template< class Pixel >
	int read( Array< Pixel, 2 > & ){};

	template< class Pixel >
	int read( Array< Pixel, 3 > & ){};
};


#endif // FILE_nbfReadFile