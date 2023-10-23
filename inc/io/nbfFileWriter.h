#ifndef FILE_nbfFileWriter
#define FILE_nbfFileWriter

#include <io/nbfFile.h>

class nbfFileWriter : public nbfFile
{

public:

	// Read data into array. Return 1 if error, 0 otherwise.
	template< class Pixel >
	int write( Array< Pixel, 2 > & ){};
	
	template< class Pixel >
	int write( Array< Pixel, 3 > & ){};
};


#endif // FILE_nbfFileWriter