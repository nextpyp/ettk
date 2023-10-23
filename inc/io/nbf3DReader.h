#pragma once

/** @file nbf3DReader.h
*	3D file reader. Part of IO.
*/

#include <fstream>
#include <vtkByteSwap.h>
#include <vtkImageData.h>

#include <io/nbfFileReader.h>

/** Interface for reading 3D datasets.
	Specialized classes provide readers for various formats.
*/
class nbf3DReader : public nbfFileReader
{

public:

	nbf3DReader(){
		this->subSample = 1;
		this->fileName = NULL;
		this->subX = 1;
		this->subY = 1;
		this->subZ = 1;
		this->Nx = 0;
		this->Ny = 0;
		this->Nz = 0;
		this->bigEndian = false;
	}

	/// Set to read big/little endian data.
	void setBigEndian(bool b){ this->bigEndian = b;}

	/// Set sub-sampling factor. If only one parameter, apply that factor to every dimension.
	void setSubSample( int _subX, int _subY = 0, int _subZ = 0 ){
		if ( _subY == 0 ){
			this->subX = _subX;
			this->subY = _subX;
			this->subZ = _subX;
		}
		else{
            this->subX = _subX;
            this->subY = _subY;
            this->subZ = _subZ;
		}
	}

	/// Set data type (for raw data only)
	void setDataType( long type ){
		this->type = type;
	}

	/// Read 3D data into array. Return 1 if error, 0 otherwise.
	template< class Pixel, int const Dim >
	int read( Array< Pixel, Dim > & );

	/// Read 3D data into VTK. Return 1 if error, 0 otherwise.
	int read( vtkImageData * );

	/// Read slice into 2D array. Return 1 if error, 0 otherwise.
	template< class Pixel >
	int readSlice( Array< Pixel, 2 > &, int, int );

	/// Read and re-write 3D data (normally subsampled). Return 1 if error, 0 otherwise.
	int rewrite( const char *);

	TinyVector< int, 3 > getDims(){ this->readHeader(); return TinyVector< int, 3 >( this->Nx, this->Ny, this->Nz ); }

	int getDimY(){ this->readHeader(); return this->Ny; }

protected:

	/// Set data dimension (for raw data only)
	void setDimensions( int _Nx, int _Ny, int _Nz ){
		this->Nx = _Nx;
		this->Ny = _Ny;
		this->Nz = _Nz;
	}
    
	/// Read file header to get data dimensions etc (format dependent).
	virtual int readHeader(){ 
		return 0;
	}

	/// Advance file pointer to data section
	virtual void advancePointerToData( ifstream & inFile ){};

	int subSample;

	/// Subsampling factors
	int subX, subY, subZ;

	/// Data dimensions
	long Nx, Ny, Nz;

	/// Pixel data type (only MRC format currently)
	int type;

	/// Data ordering
	bool bigEndian;

};

template< class Pixel, int const Dim >
int nbf3DReader :: read( Array< Pixel, Dim > & A ){

	ifstream inFile( this->fileName, ios::in | ios::binary );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		return 1;		
	}

	// Get data dimensions
	this->readHeader();

	A.resize( ceil( (float)Nx / this->subX ), 
			  ceil( (float)Ny / this->subY ),
			  ceil( (float)Nz / this->subZ ) );

	// advance to data section
	this->advancePointerToData( inFile );

	// get pixel data
	for ( int k = 0; k < Nz; k++ ){
		for ( int j = 0; j < Ny; j++ ){
			for ( int i = 0; i < Nx; i++ ){
				char dChar;
				short dShort;
                if ( ~inFile.eof() ){
					switch ( type ){
						case 0:
							inFile.read((char *) &dChar, sizeof(dChar) );
							if ( this->bigEndian == true ){
								vtkByteSwap::Swap4BE(&dChar);
							}
							if ( !(k % this->subZ) && !(j % this->subY) && !(i % this->subX) ){
								A( i / this->subX,
								   j / this->subY,
								   k / this->subZ ) = (Pixel)dChar;
							}
							break;
						case 1:
							inFile.read((char *) &dShort, sizeof(dShort) );
							if ( this->bigEndian == true ){
                                vtkByteSwap::Swap2BE(&dShort);
							}
							if ( !(k % this->subZ) && !(j % this->subY) && !(i % this->subX) ){
								A( i / this->subX,
								   j / this->subY,
								   k / this->subZ ) = (Pixel)dShort;
							}
							break;
						default:
							inFile.close();
							return 1;
					}
				}
				else{
					// bad file size
					inFile.close();
					return 1;
				}
			}
		}
		cout << k << " of " << Nz << endl << flush;
		fflush(stdout);
	}

	inFile.close();
	return 0;
}

template< class Pixel >
int nbf3DReader :: readSlice( Array< Pixel, 2 > & A, int slice, int dim ){

	ifstream inFile( this->fileName, ios::in | ios::binary );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		return 1;		
	}

	// Get data dimensions
	this->readHeader();

	cout << "Data size = [" << Nx << "," << Ny << ","  << Nz << "]." << endl;
	cout << "Extracting slice " << slice << " from the " << dim << " dimension." << endl;

	switch ( dim ){
		case firstDim:
			A.resize( ceil( (float)Ny / this->subY ), 
				      ceil( (float)Nz / this->subZ ) );
			break;
		case secondDim:
            A.resize( ceil( (float)Nx / this->subX ), 
			     	  ceil( (float)Nz / this->subZ ) );
			break;
		case thirdDim:
			A.resize( ceil( (float)Nx / this->subX ), 
					  ceil( (float)Ny / this->subY ) );
			break;
	}

	cout << A.shape() << endl;

	if ( slice < 0 || slice >= Nz ){
		inFile.close();
		return 1;		
	}

	// advance to data section
	this->advancePointerToData( inFile );

	// get pixel data
	for ( int k = 0; k < Nz; k++ ){
		for ( int j = 0; j < Ny; j++ ){
			for ( int i = 0; i < Nx; i++ ){
				char dChar;
				short dShort;
                if ( ~inFile.eof() ){
					switch ( type ){
						case 0:
							inFile.read((char *) &dChar, sizeof(dChar) );
							if ( this->bigEndian == true ){
								vtkByteSwap::Swap4BE(&dChar);
							}
							if ( !(k % this->subZ) && !(j % this->subY) && !(i % this->subX) )
							{
								switch (dim){
									case firstDim:
										if ( slice == i ){
											A( j / this->subY, k / this->subZ ) = (Pixel)dChar;
										}
										break;
									case secondDim:
										if ( slice == j ){
											A( i / this->subX, k / this->subZ ) = (Pixel)dChar;
										}
										break;
									case thirdDim:
										if ( slice == k ){
											A( i / this->subX, j / this->subY ) = (Pixel)dChar;
										}
								}
							}
							break;
						case 1:
							inFile.read((char *) &dShort, sizeof(dShort) );
							if ( this->bigEndian == true ){
								vtkByteSwap::Swap2BE(&dShort);
							}
							if ( !(k % this->subZ) && !(j % this->subY) && !(i % this->subX) )
							{
								switch (dim){
									case firstDim:
										if ( slice == i ){
											A( j / this->subY, k / this->subZ ) = (Pixel)dShort;
										}
										break;
									case secondDim:
										if ( slice == j ){
											A( i / this->subX, k / this->subZ ) = (Pixel)dShort;
										}
										break;
									case thirdDim:
										if ( slice == k ){
											A( i / this->subX, j / this->subY ) = (Pixel)dShort;
										}
								}
							}
							break;
						default:
							inFile.close();
							return 1;
					}
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

int nbf3DReader :: rewrite( const char * _filename ){

	ifstream inFile( this->fileName, ios::in | ios::binary );
	ofstream outFile( _filename, ios::out | ios::binary );
	if ( inFile.is_open() != 1 || outFile.is_open() != 1 ){
		inFile.close();
		outFile.close();
		return 1;		
	}

	// Data dimensions
	this->readHeader();

	cout << "[" << ceil( (float)Nx / this->subX ) << "," <<
		ceil( (float)Ny / this->subY ) << "," <<
		ceil( (float)Nz / this->subZ ) << "]" << endl;

	this->advancePointerToData( inFile );

	// get pixel data
	for ( int k = 0; k < Nz; k++ ){
		for ( int j = 0; j < Ny; j++ ){
			for ( int i = 0; i < Nx; i++ ){
				char dChar;
				short dShort;
                if ( inFile.eof() == false ){
					switch ( type ){
						case 0:
							inFile.read((char *) &dChar, sizeof(dChar) );
							if ( this->bigEndian == true ){
								vtkByteSwap::Swap4BE(&dChar);
							}
							if ( !(k % this->subZ) && !(j % this->subY) && !(i % this->subX) ){
								outFile.write((char *) &(dChar), sizeof(dChar) );
							}
							break;
						case 1:
							inFile.read((char *) &dShort, sizeof(dShort) );
							if ( this->bigEndian == true ){
								vtkByteSwap::Swap2BE(&dShort);
							}
							if ( !(k % this->subZ) && !(j % this->subY) && !(i % this->subX) ){
								outFile.write((char *) &(dShort), sizeof(dShort) );
							}
							break;
						default:
							inFile.close();
							outFile.close();
							return 1;
					}
				}
				else{
					// bad file size
					inFile.close();
					outFile.close();
					return 1;
				}
			}
		}
	}

	inFile.close();
	outFile.close();
	return 0;
}
