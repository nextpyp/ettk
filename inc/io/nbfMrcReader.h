#pragma once

/** @file nbfMrcReader.h
*	MRC file reader. Part of IO suite.
*/

#include <io/nbf3DReader.h>
#include <vtkImageReader.h>

/** MRC file reader.
	Reads MRC files into Blitz or VTK. 
*/
class nbfMrcReader : public nbf3DReader
{

public:

	/** Read directly into VTK image using vtkImageReader.
	*/
	void read( vtkImageData *, Array< float, 1 > &, Array< float, 1 > & );
	void read( Array< float, 1 > &, Array< float, 1 > & );
	void read( vtkImageData * );

protected:

	/// Redefine from father.
	int readHeader();

	/// Redefine from father.
	void advancePointerToData( ifstream & );

};

int nbfMrcReader :: readHeader(){

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

	if ( this->bigEndian == true ){
		vtkByteSwap::Swap4BE(&Nx);
		vtkByteSwap::Swap4BE(&Ny);
		vtkByteSwap::Swap4BE(&Nz);
	}

	this->setDimensions(Nx,Ny,Nz);

	// Get data type:
	//	0 - image : signed 8-bit bytes range -128 to 127
	//  1 - image : 16-bit halfwords
	//  2 - image : 32-bit reals
	//  3 - transform : complex 16-bit integers
	//  4 - transform : complex 32-bit reals  
	inFile.read((char *) &((this->type)), sizeof ((this->type))) ;
	if ( this->bigEndian == true ){
		vtkByteSwap::Swap4BE(&((this->type)));
	}

	inFile.close();
	return 0;
}

void nbfMrcReader :: advancePointerToData( ifstream & inFile ){
	
	long fileSize;
	inFile.seekg(0, ios::beg);
	fileSize = inFile.tellg();
	inFile.seekg(0, ios::end);
	fileSize = (int)inFile.tellg() - fileSize;
	inFile.close();

	int headerSize = fileSize - this->Nx * this->Ny * this->Nz * sizeof( this->type );

	// skip the header
	inFile.seekg( headerSize, ios_base::beg );
}

void nbfMrcReader :: read( vtkImageData * image, Array< float, 1 > & angles, Array< float, 1 > & means )
{
	vtkImageReader * mrcReader = vtkImageReader::New();
	mrcReader->SetFileName( this->fileName );
	if ( this->bigEndian == true ){
		mrcReader->SetDataByteOrderToBigEndian();
	}
	else{
		mrcReader->SetDataByteOrderToLittleEndian();
	}
	mrcReader->SetFileDimensionality(3);

	// read header manually
	ifstream inFile( this->fileName, ios::in | ios::binary );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		cout << "ERROR: Unable to open file " << this->fileName << endl;
		return;		
	}

	// get dimensions
	inFile.read((char *)&this->Nx, 4);
	inFile.read((char *)&this->Ny, 4);
	inFile.read((char *)&this->Nz, 4);
	inFile.read((char *)&this->type, 4);

	if ( this->bigEndian == true ){
		vtkByteSwap::Swap4BE(&this->Nx);
		vtkByteSwap::Swap4BE(&this->Ny);
		vtkByteSwap::Swap4BE(&this->Nz);
		vtkByteSwap::Swap4BE(&this->type);
	}

	// After the first 1024 bytes comes following data (blocks of 128 bytes, one for each tilt):
	//
	// SIZE  	DATA  	NAME  	DESCRIPTION
	// 4 	float 	a_tilt 	Alpha tilt (deg)
	// 4 	float 	b_tilt 	Beta tilt (deg)
	// 4 	float 	x_stage 	Stage x position (Unit=m. But if value>1, unit=µm)
	// 4 	float 	y_stage 	Stage y position (Unit=m. But if value>1, unit=µm)
	// 4 	float 	z_stage 	Stage z position (Unit=m. But if value>1, unit=µm)
	// 4 	float 	x_shift 	Image shift x (Unit=m. But if value>1, unit=µm)
	// 4 	float 	y_shift 	Image shift y (Unit=m. But if value>1, unit=µm)
	// 4 	float 	z_shift 	Image shift z (Unit=m. But if value>1, unit=µm)
	// 4 	float 	defocus 	Defocus Unit=m. But if value>1, unit=µm)
	// 4 	float 	exp_time 	Exposure time (s)
	// 4 	float 	mean_int 	Mean value of image
	// 4 	float 	tilt_axis 	Tilt axis (deg)
	// 4 	float 	pixel_size 	Pixel size of image (m)
	// 4 	float 	magnification 	Magnification used
	// 4 	float 	remainder 	Not used (filling up to 128 bytes)

	// Get mean values from file
	means.resize( this->Nz );
	for ( int i = 0; i < this->Nz; i++ ){
		inFile.seekg(1024+i*128+4*9, ios::beg );
		inFile.read((char *)(&(means(i))), sizeof (float) );
		if ( this->bigEndian == true ){
			vtkByteSwap::Swap4BE(&(means(i)));
		}
	}

	// Get alpha tilt angles from file
	angles.resize( this->Nz );
	for ( int i = 0; i < this->Nz; i++ ){
		inFile.seekg(1024+i*128, ios::beg );
		inFile.read((char *)(&(angles(i))), sizeof (float) );
		if ( this->bigEndian == true ){
			vtkByteSwap::Swap4BE(&(angles(i)));
		}
	}

	// get file size
	long fileSize;
	inFile.seekg(0, ios::beg);
	fileSize = inFile.tellg();
	inFile.seekg(0, ios::end);
	fileSize = (int)inFile.tellg() - fileSize;
	inFile.close();

	int headerSize;
	if ( this->type == 0 ){
		mrcReader->SetDataScalarTypeToUnsignedChar();
		headerSize = fileSize - this->Nx * this->Ny * this->Nz * 1;
	}
	if ( this->type == 1 ){
		mrcReader->SetDataScalarTypeToShort();
		headerSize = fileSize - this->Nx * this->Ny * this->Nz * 2;
	}
	if ( this->type == 2 ){
		mrcReader->SetDataScalarTypeToFloat();
		headerSize = fileSize - this->Nx * this->Ny * this->Nz * 4;
	}

	// compute header size
	//mrcReader->SetHeaderSize( 1024 );
	mrcReader->SetHeaderSize( headerSize );
	mrcReader->SetDataExtent(0,this->Nx-1,0,this->Ny-1,0,this->Nz-1);
	mrcReader->Update();
	
	image->ShallowCopy( mrcReader->GetOutput() );
	mrcReader->Delete();
}

void nbfMrcReader :: read( Array< float, 1 > & angles, Array< float, 1 > & means )
{
	// read header manually
	ifstream inFile( this->fileName, ios::in | ios::binary );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		return;		
	}

	// get dimensions
	inFile.read((char *)&this->Nx, 4);
	inFile.read((char *)&this->Ny, 4);
	inFile.read((char *)&this->Nz, 4);
	inFile.read((char *)&this->type, 4);

	if ( this->bigEndian == true ){
		vtkByteSwap::Swap4BE(&this->Nx);
		vtkByteSwap::Swap4BE(&this->Ny);
		vtkByteSwap::Swap4BE(&this->Nz);
		vtkByteSwap::Swap4BE(&this->type);
	}

	// After the first 1024 bytes comes following data (blocks of 128 bytes, one for each tilt):
	//
	// SIZE  	DATA  	NAME  	DESCRIPTION
	// 4 	float 	a_tilt 	Alpha tilt (deg)
	// 4 	float 	b_tilt 	Beta tilt (deg)
	// 4 	float 	x_stage 	Stage x position (Unit=m. But if value>1, unit=µm)
	// 4 	float 	y_stage 	Stage y position (Unit=m. But if value>1, unit=µm)
	// 4 	float 	z_stage 	Stage z position (Unit=m. But if value>1, unit=µm)
	// 4 	float 	x_shift 	Image shift x (Unit=m. But if value>1, unit=µm)
	// 4 	float 	y_shift 	Image shift y (Unit=m. But if value>1, unit=µm)
	// 4 	float 	z_shift 	Image shift z (Unit=m. But if value>1, unit=µm)
	// 4 	float 	defocus 	Defocus Unit=m. But if value>1, unit=µm)
	// 4 	float 	exp_time 	Exposure time (s)
	// 4 	float 	mean_int 	Mean value of image
	// 4 	float 	tilt_axis 	Tilt axis (deg)
	// 4 	float 	pixel_size 	Pixel size of image (m)
	// 4 	float 	magnification 	Magnification used
	// 4 	float 	remainder 	Not used (filling up to 128 bytes)

	// Get mean values from file
	means.resize( this->Nz );
	for ( int i = 0; i < this->Nz; i++ ){
		inFile.seekg(1024+i*128+4*9, ios::beg );
		inFile.read((char *)(&(means(i))), sizeof (float) );
		if ( this->bigEndian == true ){
			vtkByteSwap::Swap4BE(&(means(i)));
		}
	}

	// Get alpha tilt angles from file
	angles.resize( this->Nz );
	for ( int i = 0; i < this->Nz; i++ ){
		inFile.seekg(1024+i*128, ios::beg );
		inFile.read((char *)(&(angles(i))), sizeof (float) );
		if ( this->bigEndian == true ){
			vtkByteSwap::Swap4BE(&(angles(i)));
		}
	}

	inFile.close();
}

void nbfMrcReader :: read( vtkImageData * image )
{
	vtkImageReader * mrcReader = vtkImageReader::New();
	mrcReader->SetFileName( this->fileName );
	if ( this->bigEndian == true ){
		mrcReader->SetDataByteOrderToBigEndian();
	}
	else{
		mrcReader->SetDataByteOrderToLittleEndian();
	}
	mrcReader->SetFileDimensionality(3);

	// read header manually
	ifstream inFile( this->fileName, ios::in | ios::binary );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		cout << "ERROR: Unable to open file " << this->fileName << endl;
		return;		
	}

	// get dimensions
	inFile.read((char *)&this->Nx, 4);
	inFile.read((char *)&this->Ny, 4);
	inFile.read((char *)&this->Nz, 4);
	inFile.read((char *)&this->type, 4);
    
    // cout << this->Nx << ", " << this->Ny << ", " << this->Nx << ", " << this->type << endl;


	if ( this->bigEndian == true ){
		vtkByteSwap::Swap4BE(&this->Nx);
		vtkByteSwap::Swap4BE(&this->Ny);
		vtkByteSwap::Swap4BE(&this->Nz);
		vtkByteSwap::Swap4BE(&this->type);
	}

	// After the first 1024 bytes comes following data (blocks of 128 bytes, one for each tilt):
	//
	// SIZE  	DATA  	NAME  	DESCRIPTION
	// 4 	float 	a_tilt 	Alpha tilt (deg)
	// 4 	float 	b_tilt 	Beta tilt (deg)
	// 4 	float 	x_stage 	Stage x position (Unit=m. But if value>1, unit=µm)
	// 4 	float 	y_stage 	Stage y position (Unit=m. But if value>1, unit=µm)
	// 4 	float 	z_stage 	Stage z position (Unit=m. But if value>1, unit=µm)
	// 4 	float 	x_shift 	Image shift x (Unit=m. But if value>1, unit=µm)
	// 4 	float 	y_shift 	Image shift y (Unit=m. But if value>1, unit=µm)
	// 4 	float 	z_shift 	Image shift z (Unit=m. But if value>1, unit=µm)
	// 4 	float 	defocus 	Defocus Unit=m. But if value>1, unit=µm)
	// 4 	float 	exp_time 	Exposure time (s)
	// 4 	float 	mean_int 	Mean value of image
	// 4 	float 	tilt_axis 	Tilt axis (deg)
	// 4 	float 	pixel_size 	Pixel size of image (m)
	// 4 	float 	magnification 	Magnification used
	// 4 	float 	remainder 	Not used (filling up to 128 bytes)
	// get file size
	long fileSize;
	inFile.seekg(0, ios::beg);
	fileSize = inFile.tellg();
	inFile.seekg(0, ios::end);
	fileSize = (int)inFile.tellg() - fileSize;
	inFile.close();

	// compute header size
	int headerSize;
	if ( this->type == 0 ){
		mrcReader->SetDataScalarTypeToUnsignedChar();
		headerSize = fileSize - this->Nx * this->Ny * this->Nz * 1;
	}
	if ( this->type == 1 ){
		mrcReader->SetDataScalarTypeToShort();
		headerSize = fileSize - this->Nx * this->Ny * this->Nz * 2;
	}
	if ( this->type == 2 ){
		mrcReader->SetDataScalarTypeToFloat();
		headerSize = fileSize - this->Nx * this->Ny * this->Nz * 4;
	}
    else {
        cout << "ERROR - don't understand data type " << this->type << endl;
        return;
    }

	//headerSize = fileSize - this->Nx * this->Ny * this->Nz * sizeof( float ) - 64*4;
	headerSize = 1024;

	// compute header size
	mrcReader->SetHeaderSize( headerSize );
	mrcReader->SetDataExtent(0,this->Nx-1,0,this->Ny-1,0,this->Nz-1);
	//mrcReader->SetDataExtent(0,this->Nx*this->Ny*this->Nz-1,0,0,0,0);
	mrcReader->Update();

	image->ShallowCopy( mrcReader->GetOutput() );

	mrcReader->Delete();
}
