#pragma once

/** @file nbfMrcWriter.h
*	MRC file writer. Part of IO suite.
*/

#include <io/nbfVTKInterface.h>
#include <io/nbfFileWriter.h>
#include <vtkImageReader.h>

/** MRC file reader.
	Reads MRC files into Blitz or VTK. 
*/
class nbfMrcWriter : public nbfFileWriter
{

public:

	/** Read directly into VTK image using vtkImageReader.
	*/
	bool write( vtkImageData *, bool reverse = true );

	/// Write data as shorts directly to file (no copying of data, done specially for big files)
	bool write( Array< float, 3 > &, bool reverse = true );
	bool write( Array< short, 3 > &, bool reverse = true );

	template< class Pixel >
	bool header( ofstream &, Array< Pixel, 3 > &, int );

};

bool nbfMrcWriter :: write( vtkImageData * image, bool reverse )
{
	if ( image->GetPointData()->GetScalars()->GetDataType() == VTK_SHORT ){
		Array< short, 3 > A;
		nbfVTKInterface::vtkToBlitzReference(image,A);

		return this->write(A,reverse);
		//ofstream outFile( this->fileName, ios::out | ios::binary );
		//if ( outFile.is_open() != 1 ){
		//	outFile.close();
		//	cerr << "ERROR: File " << this->fileName << " cannot be accessed." << endl;
		//	return false;
		//}

		//long type = 1;
		//
		//if ( this->header( outFile, A, type ) == false ){
		//	outFile.close();
		//	return false;
		//}

		////// Get data dimensions
		////long Nx = A.rows();
		////long Ny = A.cols();
		////long Nz = A.depth();

		////outFile.write((char *) &Nx, 4 );
		////outFile.write((char *) &Ny, 4 );
		////outFile.write((char *) &Nz, 4 ) ;

		////long type = 1;
		////outFile.write((char *) &type, 4 ) ;

		//// advance pointer
		//
		////Array< short, 1 > header( ( 1024 - 22 * 4 ) / sizeof(short) );
		////// ignore everything exept for data dimensions
		////header = 0;

		////short * data = header.data();
		////outFile.write( (char *) data, sizeof(short) * header.numElements() );

		//// write image data
		//short * data = A.data();
		//outFile.write( (char *) data, sizeof(short) * A.numElements() );
		//if ( outFile.bad() == true ){
		//	outFile.close();
		//	return false;
		//}

		//outFile.close();
		//return true;
	}
	else if ( image->GetPointData()->GetScalars()->GetDataType() == VTK_FLOAT ){
		Array< float, 3 > A;
		nbfVTKInterface::vtkToBlitzReference(image,A);

		return this->write(A,reverse);
		//ofstream outFile( this->fileName, ios::out | ios::binary );
		//if ( outFile.is_open() != 1 ){
		//	outFile.close();
		//	cerr << "ERROR: File " << this->fileName << " cannot be accessed." << endl;
		//	return false;
		//}

		//long type = 2;
		//if ( this->header( outFile, A, type ) == false ){
		//	outFile.close();
		//	return false;
		//}

		////// Get data dimensions
		////long Nx = A.rows();
		////long Ny = A.cols();
		////long Nz = A.depth();

		////outFile.write((char *) &Nx, 4 );
		////outFile.write((char *) &Ny, 4 );
		////outFile.write((char *) &Nz, 4 ) ;

		////long type = 2;
		////outFile.write((char *) &type, 4 ) ;

		////// advance pointer
		////
		////Array< float, 1 > header( ( 1024 - 22 * 4 ) / sizeof(float) );
		////// ignore everything exept for data dimensions
		////header = 0;

		////float * data = header.data();
		////outFile.write( (char *) data, sizeof(float) * header.numElements() );

		//// write image data
		//float * data = A.data();
		//outFile.write( (char *) data, sizeof(float) * A.numElements() );
		//if ( outFile.bad() == true ){
		//	outFile.close();
		//	return false;
		//}

		//outFile.close();
		//return true;
	} else{
		cout << "WARNING - unsuported data type\n";
		return false;
	}
}

bool nbfMrcWriter :: write( Array< float, 3 > & A, bool reverse )
{
	ofstream outFile( this->fileName, ios::out | ios::binary );
	if ( outFile.is_open() != 1 ){
		outFile.close();
		cerr << "ERROR: File " << this->fileName << " cannot be accessed." << endl;
		return false;
	}

    // in 32-bits sizeof(long)=4 but in 64-bits sizeof(long)=8!
    // long type = 2;
	int type = 2;
	if ( this->header( outFile, A, type ) == false ){
		outFile.close();
		return false;
	}

	//// Get data dimensions
	//long Nx = A.rows();
	//long Ny = A.cols();
	//long Nz = A.depth();

	//outFile.write((char *) &Nx, 4 );
	//outFile.write((char *) &Ny, 4 );
	//outFile.write((char *) &Nz, 4 ) ;

	//outFile.write((char *) &type, 4 ) ;

	//// advance pointer

	//Array< short, 1 > header( ( 1024 - 22 * 4 ) / sizeof(short) );
	//// ignore everything exept for data dimensions
	//header = 0;

	//short * data = header.data();
	//outFile.write( (char *) data, sizeof(short) * header.numElements() );

	//// write image data
	//float * imagedata = A.data();
	//outFile.write( (char *) imagedata, sizeof(float) * A.numElements() );

	if ( reverse == true ){
		for ( int i = 0; i < A.depth(); i++ ){
			for ( int j = A.ubound(secondDim); j >= 0; j-- ){
				for ( int k = 0; k < A.rows(); k++ ){
					float pixel = A(k,j,i);
					outFile.write( (char *) &pixel, sizeof(float) );
					if ( outFile.bad() == true ){
						outFile.close();
						return false;
					}
				}
			}
		}
	} else {
		for ( int i = 0; i < A.depth(); i++ ){
			for ( int j = 0; j < A.cols(); j++ ){
				for ( int k = 0; k < A.rows(); k++ ){
					float pixel = A(k,j,i);
					outFile.write( (char *) &pixel, sizeof(float) );
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

template< class Pixel >
bool nbfMrcWriter :: header( ofstream & outFile, Array< Pixel, 3 > & A, int type )
{
	// Get data dimensions
	long Nx = A.rows();
	long Ny = A.cols();
	long Nz = A.depth();

	outFile.write((char *) &Nx, 4 ); if ( outFile.bad() == true ){ return false; }
	outFile.write((char *) &Ny, 4 ); if ( outFile.bad() == true ){ return false; }
	outFile.write((char *) &Nz, 4 ); if ( outFile.bad() == true ){ return false; }

	outFile.write((char *) &type, 4 ); if ( outFile.bad() == true ){ return false; }

	// Starting point of sub image.
	long tmp = 0;
	outFile.write((char *) &tmp, 4 ); if ( outFile.bad() == true ){ return false; }
	outFile.write((char *) &tmp, 4 ); if ( outFile.bad() == true ){ return false; }
	outFile.write((char *) &tmp, 4 ); if ( outFile.bad() == true ){ return false; }

	// Grid size in X, Y, and Z
	outFile.write((char *) &Nx, 4 ); if ( outFile.bad() == true ){ return false; }
	outFile.write((char *) &Ny, 4 ); if ( outFile.bad() == true ){ return false; }
	outFile.write((char *) &Nz, 4 ); if ( outFile.bad() == true ){ return false; }

	// Cell size; pixel spacing = xlen/mx
	float ftmp = 410;
	outFile.write((char *) &ftmp, 4 ); if ( outFile.bad() == true ){ return false; }
	//ftmp = Ny;
	outFile.write((char *) &ftmp, 4 ); if ( outFile.bad() == true ){ return false; }
	//ftmp = Nz;
	outFile.write((char *) &ftmp, 4 ); if ( outFile.bad() == true ){ return false; }

	// cell angles
	ftmp = 90.0;
	outFile.write((char *) &ftmp, 4 ); if ( outFile.bad() == true ){ return false; }
	outFile.write((char *) &ftmp, 4 ); if ( outFile.bad() == true ){ return false; }
	outFile.write((char *) &ftmp, 4 ); if ( outFile.bad() == true ){ return false; }

	// map coloumn 1=x,2=y,3=z. 
	tmp = 1;
	outFile.write((char *) &tmp, 4 ); if ( outFile.bad() == true ){ return false; }
	// map row 1=x,2=y,3=z. 
	tmp = 2;
	outFile.write((char *) &tmp, 4 ); if ( outFile.bad() == true ){ return false; }
	// map section 1=x,2=y,3=z. 
	tmp = 3;
	outFile.write((char *) &tmp, 4 ); if ( outFile.bad() == true ){ return false; }

	// Minimum pixel value.
	float minima = min(A( Range::all(), Range::all(), Nz / 2 ));
	outFile.write((char *) &minima, 4 ); if ( outFile.bad() == true ){ return false; }

	// Maximum pixel value.
	float maxima = max(A( Range::all(), Range::all(), Nz / 2 ));
	outFile.write((char *) &maxima, 4 ); if ( outFile.bad() == true ){ return false; }

	// Mean pixel value.
	float meanv = mean(A( Range::all(), Range::all(), Nz / 2 ));
	outFile.write((char *) &meanv, 4 ); if ( outFile.bad() == true ){ return false; }

	Array< short, 1 > header( ( 1024 - 22 * 4 ) / sizeof(short) );
	// ignore everything exept for data dimensions
	header = 0;

	short * data = header.data();
	outFile.write( (char *) data, sizeof(short) * header.numElements() );

	return !outFile.bad();
}

bool nbfMrcWriter :: write( Array< short, 3 > & A, bool reverse )
{
	ofstream outFile( this->fileName, ios::out | ios::binary );
	if ( outFile.is_open() != 1 ){
		outFile.close();
		return false;
	}

	int type = 1;

	if ( this->header( outFile, A, type ) == false ){
		outFile.close();
		return false;
	}

	//// Get data dimensions
	//long Nx = A.rows();
	//long Ny = A.cols();
	//long Nz = A.depth();

	//outFile.write((char *) &Nx, 4 );
	//outFile.write((char *) &Ny, 4 );
	//outFile.write((char *) &Nz, 4 ) ;

	//long type = 1;
	//outFile.write((char *) &type, 4 ) ;

	//// advance pointer

	//Array< short, 1 > header( ( 1024 - 22 * 4 ) / sizeof(short) );
	//// ignore everything exept for data dimensions
	//header = 0;

	//short * data = header.data();
	//outFile.write( (char *) data, sizeof(short) * header.numElements() );

	//// write image data
	//float * imagedata = A.data();
	//outFile.write( (char *) imagedata, sizeof(float) * A.numElements() );

	if ( reverse == true ){
		for ( int i = 0; i < A.depth(); i++ ){
			for ( int j = A.ubound(secondDim); j >= 0; j-- ){
				for ( int k = 0; k < A.rows(); k++ ){
					short pixel = A(k,j,i);
					outFile.write( (char *) &pixel, sizeof(short) );
					if ( outFile.bad() == true ){ outFile.close();	return false; }
				}
			}
		}
	} else {
		for ( int i = 0; i < A.depth(); i++ ){
			for ( int j = 0; j <= A.cols(); j++ ){
				for ( int k = 0; k < A.rows(); k++ ){
					short pixel = A(k,j,i);
					outFile.write( (char *) &pixel, sizeof(short) );
					if ( outFile.bad() == true ){ outFile.close();	return false; }
				}
			}
		}
	}

	outFile.close();
	return true;
}
