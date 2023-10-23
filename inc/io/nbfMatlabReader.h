#pragma once

/** @file nbfMatlabReader.h
	Matlab file reader. Part of IO suite.
*/

#include <io/nbfFileReader.h>

/** Matlab file reader.
	Read Matlab files into Blitz. 

	@see nbfMatlabWriter
*/
class nbfMatlabReader : public nbfFileReader
{

public:

	template< class Pixel >
	int read( Array< Pixel, 2 > & );

	template< class Pixel >
	int read( Array< Pixel, 3 > & );

	template< class Pixel >
	int read( Array< Pixel, 4 > & );

	nbfMatlabReader(){
		this->fileName = NULL;
	}

};

template< class Pixel >
int nbfMatlabReader :: read( Array< Pixel, 2 > & A ){

	ifstream inFile( this->fileName, ios::in | ios::binary );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		return 1;		
	}

	// Get data dimensions
	int Nx, Ny;
	inFile.read((char *) &Nx, sizeof (Nx)) ;
	inFile.read((char *) &Ny, sizeof (Ny)) ;
	
	A.resize( Nx, Ny );

	// get pixel data
	for ( int j = 0; j < Ny; j++ ){
		for ( int i = 0; i < Nx; i++ ){
			Pixel dFloat;
			if ( !inFile.eof() ){
				inFile.read((char *) &dFloat, sizeof(dFloat) );
				A( i , j ) = (Pixel)dFloat;
			}
			else{
				// bad file size
				inFile.close();
				return 1;
			}
		}
	}

	inFile.close();
	return 0;
}

template< class Pixel >
int nbfMatlabReader :: read( Array< Pixel, 3 > & A ){

	ifstream inFile( this->fileName, ios::in | ios::binary );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		return 1;		
	}

	// Get data dimensions
	int Nx, Ny, Nz;
	inFile.read((char *) &Nx, sizeof (Nx)) ;
	inFile.read((char *) &Ny, sizeof (Ny)) ;
	inFile.read((char *) &Nz, sizeof (Nz)) ;
	
	A.resize( Nx, Ny, Nz );

	// get pixel data
	for ( int k = 0; k < Nz; k++ ){
		for ( int j = 0; j < Ny; j++ ){
			for ( int i = 0; i < Nx; i++ ){
				Pixel dFloat;
				if ( !inFile.eof() ){
					inFile.read((char *) &dFloat, sizeof(dFloat) );
					A( i , j, k ) = (Pixel)dFloat;
				}
				else{
					// bad file size
					inFile.close();
					return 1;
				}
			}
		}
	}

	inFile.close();
	return 0;
}

template< class Pixel >
int nbfMatlabReader :: read( Array< Pixel, 4 > & A ){

	ifstream inFile( this->fileName, ios::in | ios::binary );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		return 1;		
	}

	// Get data dimensions
	int Nx, Ny, Nz, Nw;
	inFile.read((char *) &Nx, sizeof (Nx)) ;
	inFile.read((char *) &Ny, sizeof (Ny)) ;
	inFile.read((char *) &Nz, sizeof (Nz)) ;
	inFile.read((char *) &Nw, sizeof (Nw)) ;
	
	A.resize( Nx, Ny, Nz, Nw );

	// get pixel data
	for ( int p = 0; p < Nw; p++ ){
		for ( int k = 0; k < Nz; k++ ){
			for ( int j = 0; j < Ny; j++ ){
				for ( int i = 0; i < Nx; i++ ){
					Pixel dFloat;
					if ( !inFile.eof() ){
						inFile.read((char *) &dFloat, sizeof(dFloat) );
						A( i , j, k, p ) = (Pixel)dFloat;
					}
					else{
						// bad file size
						inFile.close();
						return 1;
					}
				}
			}
		}
	}

	inFile.close();
	return 0;
}