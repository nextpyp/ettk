#include <blitz/array.h>

using namespace blitz;

#include <vector>

#include <Flujos.hh>
#include <NarrowBand/nbfMatlabReader.h>
#include <NarrowBand/nbfMatlabWriter.h>
#include <IO/FlujosIO.hh>

void main( int argv, char ** argc ){

	Array< float, 2 > A;
	nbfMatlabReader reader;
	reader.setFileName( argc[1] );
	if ( reader.read( A ) ){
		cout << "Error reading file.\n";
	}

	nbfMatlabWriter writer;
	writer.setFileName( argc[2] );
	writer.write( A );

	cout << A(Range(10,14),Range(15,18)) << endl;

	SaveToFile(A,"full.blitz");
}
