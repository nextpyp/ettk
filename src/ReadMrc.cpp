#include <blitz/array.h>

using namespace blitz;

#include <io/nbfMrcReader.h>
#include <io/nbfMatlabWriter.h>

void main( int argv, char ** argc ){

	Array< short, 3 > A;
	nbfMrcReader reader;
//	reader.setBigEndian(true);

	reader.setFileName( argc[1] );

	reader.setSubSample( 2, 1, 2 );

#if 0
	// read a single slice ( slice number, slicing dimension )
	if ( reader.readSlice( A, atoi( argc[2] ), atoi( argc[3] ) ) ){
		cout << "Error reading file.\n";
	}
#elif 0
	// read the full volume
	if ( reader.read( A ) ){
		cout << "Error reading file.\n";
	}
#else
	// re-writes the file only sub-sampled
	if ( reader.rewrite( argc[2] ) ){
		cout << "Error reading file.\n";
	}
#endif

	cout << A.shape() << endl;
	cout << "[" << max(A) << "," << min(A) << "]" << endl;

	nbfMatlabWriter writer;

	//writer.setFileName( argc[4] );
	//writer.write( A );

#if 0 // write individual slices
	for ( int i = 0; i < A.cols(); i++ ){
		// build output file name: suffix + number + .array
		string outFile;
		outFile = argc[2];
		if ( i < 10 )
			outFile += "00";
		else if ( i < 100 )
			outFile += '0';
		char nada[20];
		sprintf(nada,"%d",(int)i);
		outFile += nada;
		outFile += ".array";
#if 0 // array or bmp
		writer.setFileName( outFile.c_str() );
		cout << "Writing " << outFile.c_str() << " ..." << endl;
	    writer.write( A(Range::all(),i,Range::all()) );
#else
		char tmp[256];
		strcpy(tmp,outFile.c_str());
		WriteArray( A(Range::all(),i,Range::all()), tmp );
#endif // array or bmp
		break;
	}
#endif // write individual slices
}
