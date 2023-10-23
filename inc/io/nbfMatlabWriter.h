#pragma once

/** @file nbfMatlabWriter.h
*	Matlab file writer. Part of IO suite.
*/

#include <io/nbfFileWriter.h>

/** Matlab file writer.
	Write Matlab files from Blitz. 

	@see nbfMatlabReader
*/
class nbfMatlabWriter : public nbfFileWriter
{

public:

	template< class Pixel >
	bool write( Array< Pixel, 1 > & );

	template< class Pixel >
	bool write( Array< Pixel, 2 > & );
	
	template< class Pixel >
	bool write( Array< Pixel, 3 > & );

	template< class Pixel >
	bool write( Array< Pixel, 4 > & );

	template< class Pixel, int const Dim >
	bool write( Array< Pixel, Dim > & );

	nbfMatlabWriter(){
		this->fileName = NULL;
	}
};


//template< class Pixel >
//int nbfMatlabWriter :: write( Array< Pixel, 1 > & A ){
//
//	ofstream outFile( this->fileName, ios::out | ios::binary );
//	if ( outFile.is_open() != 1 ){
//		outFile.close();
//		return 1;		
//	}
//
//	// Get data dimensions
//	int D = 1;
//	int Nx = A.rows();
//
//	outFile << D << "\t" << Nx << endl;
//	
//	// save pixel data
//	for ( int i = 0; i < Nx; i++ ){
//		outFile << A(i);
//		if ( i < Nx - 1 ){
//			outFile << "\t";
//		}
//	}
//
//	outFile.close();
//	return 0;
//}
//
//template< class Pixel >
//int nbfMatlabWriter :: write( Array< Pixel, 2 > & A ){
//
//	ofstream outFile( this->fileName, ios::out );
//	if ( outFile.is_open() != 1 ){
//		outFile.close();
//		return 1;		
//	}
//
//	// Get data dimensions
//	int Nx = A.rows();
//	int Ny = A.cols();
//
//	outFile << Nx << "\t" << Ny << endl;
//	
//	// save pixel data
//	for ( int j = 0; j < Ny; j++ ){
//		for ( int i = 0; i < Nx; i++ ){
//			outFile << A(i,j);
//			if ( i < Nx - 1 ){
//				outFile << "\t";
//			}
//		}
//		outFile << endl;
//	}
//
//	outFile.close();
//	return 0;
//}
//
//template< class Pixel >
//int nbfMatlabWriter :: write( Array< Pixel, 3 > & A ){
//
//	ofstream outFile( this->fileName, ios::out | ios::binary );
//	if ( outFile.is_open() != 1 ){
//		outFile.close();
//		return 1;		
//	}
//
//	// Get data dimensions
//	int Nx = A.rows();
//	int Ny = A.cols();
//	int Nz = A.depth();
//
//	outFile.write((char *) &Nx, sizeof (Nx)) ;
//	outFile.write((char *) &Ny, sizeof (Ny)) ;
//	outFile.write((char *) &Nz, sizeof (Nz)) ;
//	
//	// save pixel data
//	for ( int k = 0; k < Nz; k++ ){
//		for ( int j = 0; j < Ny; j++ ){
//			for ( int i = 0; i < Nx; i++ ){
//				float dFloat = A( i , j, k );
//					outFile.write((char *) &(A( i , j, k )), sizeof(A( i , j, k )) );
//			}
//		}
//	}
//
//	outFile.close();
//	return 0;
//}

template< class Pixel >
bool nbfMatlabWriter :: write( Array< Pixel, 1 > & A ){

	ofstream outFile( this->fileName, ios::out | ios::binary );
	if ( outFile.is_open() != 1 ){
		outFile.close();
		cerr << "ERROR: File " << this->fileName << " is already open." << endl;
		return false;		
	}

	// Get data dimensions
	int Nx = A.rows();

	outFile.write((char *) &Nx, sizeof (Nx)) ;

	if ( outFile.bad() == true ){
		outFile.close();
		return false;
	}

	// save pixel data
	for ( int i = 0; i < Nx; i++ ){
		float dFloat = A( i );
		outFile.write((char *) &(A( i )), sizeof(A( i )) );
		if ( outFile.bad() == true ){
			outFile.close();
			return false;
		}
	}

	outFile.close();
	return true;
}

template< class Pixel >
bool nbfMatlabWriter :: write( Array< Pixel, 2 > & A ){

	ofstream outFile( this->fileName, ios::out | ios::binary );
	if ( outFile.is_open() != 1 ){
		outFile.close();
		cerr << "ERROR: File " << this->fileName << " is already open." << endl;
		return false;		
	}

	// Get data dimensions
	int Nx = A.rows();
	int Ny = A.cols();

	outFile.write((char *) &Nx, sizeof (Nx)) ;
	if ( outFile.bad() == true ){
		outFile.close();
		return false;
	}
	outFile.write((char *) &Ny, sizeof (Ny)) ;
	if ( outFile.bad() == true ){
		outFile.close();
		return false;
	}
	
	// save pixel data
	for ( int j = 0; j < Ny; j++ ){
		for ( int i = 0; i < Nx; i++ ){
			float dFloat = A( i , j );
			outFile.write((char *) &(A( i , j )), sizeof(A( i , j )) );
			if ( outFile.bad() == true ){
				outFile.close();
				return false;
			}
		}
	}

	outFile.close();
	return true;
}

template< class Pixel >
bool nbfMatlabWriter :: write( Array< Pixel, 3 > & A ){

	ofstream outFile( this->fileName, ios::out | ios::binary );
	if ( outFile.is_open() != 1 ){
		outFile.close();
		cerr << "ERROR: File " << this->fileName << " is already open." << endl;
		return false;
	}

	// Get data dimensions
	int Nx = A.rows();
	int Ny = A.cols();
	int Nz = A.depth();

	outFile.write((char *) &Nx, sizeof (Nx)) ;
	if ( outFile.bad() == true ){
		outFile.close();
		return false;
	}
	outFile.write((char *) &Ny, sizeof (Ny)) ;
	if ( outFile.bad() == true ){
		outFile.close();
		return false;
	}
	outFile.write((char *) &Nz, sizeof (Nz)) ;
	if ( outFile.bad() == true ){
		outFile.close();
		return false;
	}
	
	// save pixel data
	for ( int k = 0; k < Nz; k++ ){
		for ( int j = 0; j < Ny; j++ ){
			for ( int i = 0; i < Nx; i++ ){
				float dFloat = A( i , j, k );
				outFile.write((char *) &(A( i , j, k )), sizeof(A( i , j, k )) );
				if ( outFile.bad() == true ){
					outFile.close();
					return false;
				}
			}
		}
	}

	outFile.close();
	return true;
}

template< class Pixel >
bool nbfMatlabWriter :: write( Array< Pixel, 4 > & A ){

	ofstream outFile( this->fileName, ios::out | ios::binary );
	if ( outFile.is_open() != 1 ){
		outFile.close();
		cerr << "ERROR: File " << this->fileName << " is already open." << endl;
		return false;
	}

	// Get data dimensions
	int Nx = A.rows();
	int Ny = A.cols();
	int Nz = A.depth();
	int Nw = A.extent(fourthDim);

	outFile.write((char *) &Nx, sizeof (Nx)) ;
	if ( outFile.bad() == true ){
		outFile.close();
		return false;
	}

	outFile.write((char *) &Ny, sizeof (Ny)) ;
	if ( outFile.bad() == true ){
		outFile.close();
		return false;
	}

	outFile.write((char *) &Nz, sizeof (Nz)) ;
	if ( outFile.bad() == true ){
		outFile.close();
		return false;
	}

	outFile.write((char *) &Nw, sizeof (Nw)) ;
	if ( outFile.bad() == true ){
		outFile.close();
		return false;
	}
	
	// save pixel data
	for ( int w = 0; w < Nw; w++ ){
		for ( int k = 0; k < Nz; k++ ){
			for ( int j = 0; j < Ny; j++ ){
				for ( int i = 0; i < Nx; i++ ){
					float dFloat = A( i, j, k, w );
					outFile.write((char *) &(A( i, j, k, w )), sizeof(A( i, j, k, w )) );
					if ( outFile.bad() == true ){
						outFile.close();
						return false;
					}
				}
			}
		}
	}

	outFile.close();
	return true;
}

//template< class Pixel, int const Dim >
//int nbfMatlabWriter :: write( Array< Pixel, Dim > & A ){
//
//	Array< Pixel, Dim > B();
//	B.resize( A.shape() );
//
//	ofstream outFile( this->fileName, ios::out | ios::binary );
//	if ( outFile.is_open() != 1 ){
//		outFile.close();
//		return 1;		
//	}
//
//	// save data dimensions
//	outFile.write((char *) &Dim, sizeof (Dim)) ;
//
//	for ( int i = 0; i < Dim; i++ ){
//		int c = A.extent(i);
//		outFile.write((char *) &c, sizeof (c)) ;
//	}
//
//	// data type: int 0, float 1, double 2
//	int type;
//	if ( sizeOf(Pixel) == sizeOf(int) ){
//		type = 0;
//	} else if ( sizeOf(Pixel) == sizeOf(float) ){
//		type = 1;
//	} else if ( sizeOf(Pixel) == sizeOf(double) ){
//		type = 2;
//	}
//	outFile.write((char *) &type, sizeof (type)) ;
//
//	Array< Pixel, Dim > :: iterator iter = A.begin();
//	while ( iter != A.end() ){
//		outFile.write((char *) &(*iter), sizeof(Pixel) );
//	}
//
//	outFile.close();
//	return 0;
//}