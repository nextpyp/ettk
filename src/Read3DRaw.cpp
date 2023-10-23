#include <blitz/array.h>
using namespace blitz;

#include <io/nbf3DReader.h>
#include <io/nbfImageWriter.h>

void main( int argv, char ** argc ){

	Array< short, 3 > A;
	nbf3DReader reader;

	// set mode to big endian if needed
	// reader.setBigEndian(true);

	// specify volume size (normally from file name)
	reader.setDimensions(1024,235,1024);

	// set type to short (signed) 16 bits
	reader.setDataType(1);

	reader.setFileName( argc[1] );

	reader.setSubSample( 1, 1, 1 );

	// read the full volume
	if ( reader.read( A ) ){
		cout << "Error reading file.\n";
	}

	cout << A.shape() << endl;

	cout << "[" << max(A) << "," << min(A) << "]" << endl;

	nbfImageWriter writer;

	writer.setFileName( argc[2] );
    writer.writeFastShort( A );
}