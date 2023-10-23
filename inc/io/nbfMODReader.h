#pragma once

/** @file nbfMODReader.h
*	MRC file reader. Part of IO suite.
*/

#include <io/nbf3DReader.h>
#include <vtkImageReader.h>

#include <vector>

/** MOD file reader.
	Read MOD files into Blitz Array. Only suports a single contour with multiple points in it. 
*/
class nbfMODReader : public nbf3DReader
{

public:

	/** Read 3D coordinates from IMOD contour data.
	*/
	void read( Array< float, 2 > & );
	
	void readRaw( Array< float, 2 > & );
	void readRawEul( Array< float, 2 > & );

protected:

};

void nbfMODReader :: read( Array< float, 2 > & A ){

	// open file for reading (assume BINARY and always BIG ENDIAN)
	ifstream inFile( this->fileName, ios::in | ios::binary );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		return;		
	}

	// Get file format description
	char format[4];
	inFile.read((char *) &format, sizeof (format)) ;

	// Get version format description
	char version[4];
	inFile.read((char *) &version, sizeof (version)) ;

	// Get object description
	inFile.seekg(8+26*4+128, ios::beg);
	char object[4];
	inFile.read((char *) &object, sizeof (object));

	inFile.seekg(176,ios::cur);
	char cont[4];
	inFile.read((char *) &cont, sizeof (cont));

	long psize, flags;
	int type, surf;

	// Get number of points in contour
	inFile.read( (char *) &psize, sizeof (psize) );
	inFile.read( (char *) &flags, sizeof (flags) );
	inFile.read( (char *) &type, sizeof (type) );
	inFile.read( (char *) &surf, sizeof (surf) );

	// Assume file is always BIG ENDIAN
	vtkByteSwap::Swap4BE(&psize);

	A.resize( psize, 3 );
	float dummy;
	for ( int i = 0; i < psize; i++ ){
		for ( int comp = 0; comp < 3; comp++ ){
			inFile.read( (char *) &dummy, sizeof(dummy));
			vtkByteSwap::Swap4BE(&dummy);
			A( i, comp ) = dummy;
		}
	}

	inFile.close();
}

void nbfMODReader :: readRaw( Array< float, 2 > & A ){

	// open file for reading (assume BINARY and always BIG ENDIAN)
	//ifstream inFile( this->fileName, ios::in | ios::binary );
	ifstream inFile( this->fileName, ios::in );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		return;		
	}

	vector< float > points;
	float dummy;

	while ( true ){
		inFile >> dummy;
		points.push_back( dummy );
		if ( inFile.fail() ){
			break;
		}
	}

	//do {
	//	inFile >> dummy;
	//	points.push_back( dummy );
	//} while ( !inFile.eof() );

	inFile.close();

	int csize = ( points.size() - 1 ) / 3;
	A.resize( csize, 3 );
	int count = 0;
	for ( int i = 0; i < csize; i++ ){
		for ( int comp = 0; comp < 3; comp++ ){
			A(i,comp) = *( points.begin() + count );
			count++;
		}
	}
}

void nbfMODReader :: readRawEul( Array< float, 2 > & A ){

	// open file for reading (assume BINARY and always BIG ENDIAN)
	//ifstream inFile( this->fileName, ios::in | ios::binary );
	ifstream inFile( this->fileName, ios::in );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		return;		
	}

	vector< float > points;
	float dummy;

	while ( true ){
		inFile >> dummy;
		points.push_back( dummy );
		if ( inFile.fail() ){
			break;
		}
	}

	//do {
	//	inFile >> dummy;
	//	points.push_back( dummy );
	//} while ( !inFile.eof() );

	inFile.close();

	int csize = ( points.size() - 1 ) / 6;
	A.resize( csize, 6 );
	int count = 0;
	for ( int i = 0; i < csize; i++ ){
		for ( int comp = 0; comp < 6; comp++ ){
			A(i,comp) = *( points.begin() + count );
			count++;
		}
	}
}