#pragma once

#include <vtkImageData.h>
#include <vtkTransform.h>
#include <vtkImageResample.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageMedian3D.h>

#include <io/nbfMrcReader.h>
#include <io/nbfMrcWriter.h>
#include <em/nbfWedgedImage3D.h>

/** Interface for VTK-like input-output pipeline filters.
	Update state is kept internally so execution is only done when needed.
	User is responsible for changing the state when doing changes that affect the filter's output.
*/
template< class Pixel >
class nbfWedgedSubImage3D : public nbfWedgedImage3D< Pixel >
{
public:

	nbfWedgedSubImage3D();

	// copy constructor
	nbfWedgedSubImage3D(const nbfWedgedSubImage3D&);

	virtual ~nbfWedgedSubImage3D();

	nbfWedgedSubImage3D< Pixel > & operator= ( const nbfWedgedSubImage3D< Pixel > & );

	/// Redefine from father
	void getWedgeImage( Array< Pixel, 3 > &, vtkTransform * = NULL );
	void getSphericalWedgeImage( Array< Pixel, 2 > &, vtkTransform * = NULL, nbfProjectionRotationMetric3D< Pixel > * = NULL );

	void getWedgeImageHalf( Array< Pixel, 3 > &, vtkTransform * = NULL );

	/// Redefine from father. Get extracted sub volume.
	void getImage( vtkImageData *, vtkTransform * = NULL, bool = false );
	void getImage( Array< Pixel, 3 > &, vtkTransform * = NULL, bool = false );
	void setImage( vtkImageData * );

	/// Provide input tomogram data to avoid reading multiple copies
	vtkImageData * getVolume(){ return this->tomogram; }

	TinyVector< Pixel, 3 > getDimensions(){ return ( this->geometry * this->magnification ); }
	void setDimensions( TinyVector< int, 3 > & d ){ this->geometry = d; }

	/// Set sampling magnification (same for all three dimensions)
	void setMagnification( TinyVector< Pixel, 3 > & p ){ this->magnification = p; }

	/// Read and write to file
	void serialize( stringstream & );
	void unserialize( stringstream & );

	static void writeHeader( stringstream & );
	static void readHeader( stringstream & );

	static void write( string, vector< nbfWedgedSubImage3D< Pixel > > & );
	static void read( string, vector< nbfWedgedSubImage3D< Pixel > > & );

	/// Access to attributes

	string getFileName(){ return this->fileName; }
	
	void setFileName( char * a ){ 
		this->fileName = a; 
	}

	void setFileName( string a ){ 
		this->fileName = a; 
	}

	TinyVector< Pixel, 3 > getPosition(){ return this->position; }
	void setPosition( TinyVector< Pixel, 3 > & a ){ this->position = a; }

	TinyVector< Pixel, 3 > getNormal(){ return this->normal; }
	void setNormal( TinyVector< Pixel, 3 > & a ){ 
		Pixel n[3]; n[0]=a[0], n[1]=a[1], n[2]=a[2]; 
		// vtkMath::Normalize(n); 
		this->normal = TinyVector<Pixel,3>(n[0],n[1],n[2]); 
	}

	//TinyVector< Pixel, 3 > getTranslation(){ return this->translation; }
	//void setTranslation( TinyVector< Pixel, 3 > & a ){ this->translation = a; }

	//TinyVector< Pixel, 3 > getRotation(){ return this->rotation; }
	//void setRotation( TinyVector< Pixel, 3 > & a ){ this->rotation = a; }

	/// VOIs are rooted at the specified point locations along the normal direction.
	/// A positive offset will displace the cut in the negative direction of the normal (below the cut position).
	/// A negative offset will displace the cut in the direction of the normal (above the cut position).
	void setCutOffset( Pixel p ){ this->cutOffset = p; }
	Pixel getCutOffset(){ return this->cutOffset; }

	void setTransform( vtkTransform * t ){
		if ( t != NULL ){
			vtkMatrix4x4::DeepCopy( this->matrix, t->GetMatrix() );
		} else {
			vtkMatrix4x4::Identity( this->matrix );
		}
	}
	
	void setTransform( vtkMatrix4x4 * m ){
		if ( m != NULL ){
			vtkMatrix4x4::DeepCopy( this->matrix, m );
		} else {
			vtkMatrix4x4::Identity( this->matrix );
		}
	}

	void getTransform( vtkTransform * m ){
		if ( m != NULL ){
			vtkTransform * t = vtkTransform::New();
			m->DeepCopy( t );
			t->Delete();
			m->Concatenate( this->matrix );
		}
	}

	int getTypeId(){ return NBF_WEDGED_SUB_IMAGE_3D; }

	void setFixedImage( Array< Pixel, 3 > & A ){
		if ( this->fixedImage == NULL ){
			this->fixedImage = vtkImageData::New();
		}
		nbfVTKInterface::blitzToVtk( A, this->fixedImage );
		TinyVector< int, 3 > shape( A.shape() );
		this->setDimensions( shape );
		TinyVector< Pixel, 3 > pos = shape / 2.0;
		this->setPosition( pos );
		this->fileName = "NA";
	}

	void setFixedImage( vtkImageData * data ){
		if ( this->fixedImage == NULL ){
			this->fixedImage = vtkImageData::New();
		}
		this->fixedImage->DeepCopy( data );
		int dims[3];
		this->fixedImage->GetDimensions(dims);
		TinyVector< int, 3 > shape( dims[0], dims[1], dims[2] );
		this->setDimensions( shape );
		TinyVector< Pixel, 3 > pos = shape / 2.0;
		this->setPosition( pos );
		this->fileName = "NA";
	}

protected:

	/// Atributes.

	// tomogram file name
	string fileName;

	// position of feature in tomogram
	TinyVector< Pixel, 3 > position;

	// volume dimensions
	TinyVector< Pixel, 3 > geometry;

	// normal orientation of feature in tomogram
	TinyVector< Pixel, 3 > normal;

	double matrix[16];

	// magnification
	TinyVector< Pixel, 3 > magnification;

	// store offset value to cutout volumes
	Pixel cutOffset;

	// pipeline for cutting volumes
	static vtkImageData * tomogram;
	static vtkImageData * box;
	static vtkImageChangeInformation * change;
	static vtkImageCast * cast;
	static vtkImageChangeInformation * changeGeometry;
	static vtkImageReslice * reslice;
	static vtkImageResample * resample;

	static string staticFileName;
	static int instances;

	vtkImageData * fixedImage;
};

// initialize static attributes

#define WEDGEDSUBIMAGE3D( type ) \
template<> int nbfWedgedSubImage3D< type > :: instances = 0; \
template<> string nbfWedgedSubImage3D< type > :: staticFileName = ""; \
template<> vtkImageData * nbfWedgedSubImage3D< type > :: tomogram = NULL; \
template<> vtkImageData * nbfWedgedSubImage3D< type > :: box = NULL; \
template<> vtkImageChangeInformation * nbfWedgedSubImage3D< type > :: change = NULL; \
template<> vtkImageCast * nbfWedgedSubImage3D< type > :: cast = NULL; \
template<> vtkImageChangeInformation * nbfWedgedSubImage3D< type > :: changeGeometry = NULL; \
template<> vtkImageReslice * nbfWedgedSubImage3D< type > :: reslice = NULL; \
template<> vtkImageResample * nbfWedgedSubImage3D< type > :: resample = NULL;

WEDGEDSUBIMAGE3D(float)
WEDGEDSUBIMAGE3D(double)

template< class Pixel >
nbfWedgedSubImage3D< Pixel > :: nbfWedgedSubImage3D()
{
	if ( this->box == NULL ){
		this->box = vtkImageData::New();
	}
	if ( this->change == NULL ){
		this->change = vtkImageChangeInformation::New();
	}
	if ( this->cast == NULL ){
		this->cast = vtkImageCast::New();
	}
	if ( this->changeGeometry == NULL ){
		this->changeGeometry = vtkImageChangeInformation::New();
	}
	if ( this->reslice == NULL ){
		this->reslice = vtkImageReslice::New();
	}
	if ( this->resample == NULL ){
		this->resample = vtkImageResample::New();
	}
	if ( this->tomogram == NULL ){
		this->tomogram = vtkImageData::New();
	}

	this->fileName = "NA";

	// initialization
	this->reslice->SetInput( this->change->GetOutput() );
	this->reslice->SetInformationInput( this->changeGeometry->GetOutput() );
	this->reslice->SetInterpolationModeToCubic();
	
	this->cast->SetOutputScalarTypeToDouble();
	this->cast->SetInput( this->reslice->GetOutput() );	

	this->resample->SetInput( this->cast->GetOutput() );
	this->resample->SetDimensionality( 3 );
	this->changeGeometry->SetInput( this->box );

	this->instances++;

	this->position = 0;
	this->geometry = 0;
	this->normal = 0; // this->normal[2] = 1;

	vtkMatrix4x4::Identity( this->matrix );
	
	this->magnification = 1;
	
	this->cutOffset = 0;

	this->fixedImage = NULL;
}

template< class Pixel >
nbfWedgedSubImage3D< Pixel > :: nbfWedgedSubImage3D( const nbfWedgedSubImage3D & c ){
	this->instances++;
	this->fixedImage = NULL;
	(*this) = c;
}
	
template< class Pixel >
nbfWedgedSubImage3D< Pixel > :: ~nbfWedgedSubImage3D()
{
	this->instances--;

	if ( this->instances == 0 ){
		if ( this->tomogram != NULL ){
			this->tomogram->Delete();
			this->tomogram = NULL;
		}
		if ( this->box != NULL ){
			this->box->Delete();
			this->box = NULL;
		}
		if ( this->change != NULL ){
			this->change->Delete();
			this->change = NULL;
		}
		if ( this->cast != NULL ){
			this->cast->Delete();
			this->cast = NULL;
		}
		if ( this->changeGeometry != NULL ){
			this->changeGeometry->Delete();
			this->changeGeometry = NULL;
		}
		if ( this->reslice != NULL ){
			this->reslice->Delete();
			this->reslice = NULL;
		}
		if ( this->resample != NULL ){
			this->resample->Delete();
			this->resample = NULL;
		}
	}

	if ( this->fixedImage != NULL ){
		this->fixedImage->Delete();
		this->fixedImage = NULL;
	}
}

template< class Pixel >
nbfWedgedSubImage3D< Pixel > & nbfWedgedSubImage3D< Pixel > :: operator= ( const nbfWedgedSubImage3D< Pixel > & param )
{
	this->fileName = param.fileName;
	this->position = param.position;
	this->geometry = param.geometry;
	this->normal = param.normal;
	this->magnification = param.magnification;
	this->cutOffset = param.cutOffset;

	for ( int i = 0; i < 16; i++ ){
		this->matrix[i] = param.matrix[i];
	}

	this->wedge = param.wedge;

	if ( param.fixedImage != NULL ){
		if ( this->fixedImage == NULL ){
			this->fixedImage = vtkImageData :: New();
		}
		this->fixedImage->DeepCopy( param.fixedImage );
	} else {
		if ( this->fixedImage != NULL ){
			this->fixedImage->Delete();
			this->fixedImage = NULL;
		}
	}
		
	//if ( this->fixedImage != NULL ){
	//	this->fixedImage->Delete();
	//	this->fixedImage = NULL;
	//}

	return *this;
}

template< class Pixel >
void nbfWedgedSubImage3D< Pixel > :: writeHeader( stringstream & output){
	output << "number\t"
		   << "lwedge\t"
		   << "uwedge\t"
		   << "posX\t"
		   << "posY\t"
		   << "posZ\t"
		   << "geomX\t"
		   << "geomY\t"
		   << "geomZ\t"
		   << "normalX\t"
		   << "normalY\t"
		   << "normalZ\t"
		   << "matrix[0]\t"
		   << "matrix[1]\t"
		   << "matrix[2]\t"
		   << "matrix[3]\t"
		   << "matrix[4]\t"
		   << "matrix[5]\t"
		   << "matrix[6]\t"
		   << "matrix[7]\t"
		   << "matrix[8]\t"
		   << "matrix[9]\t"
		   << "matrix[10]\t"
		   << "matrix[11]\t"
		   << "matrix[12]\t"
		   << "matrix[13]\t"
		   << "matrix[14]\t"
		   << "matrix[15]\t"
		   << "magnification[0]\t"
		   << "magnification[1]\t"
		   << "magnification[2]\t"
		   << "cutOffset\t"
		   << "filename" << endl;
}

template< class Pixel >
void nbfWedgedSubImage3D< Pixel > :: readHeader( stringstream & input ){
	string tmp;
	for ( int i = 0; i < 33; i++ ){ 
		input >> tmp;
	}
}

template< class Pixel >
void nbfWedgedSubImage3D< Pixel > :: serialize( stringstream & output ){
	output << this->wedge.lwedge << "\t"
		   << this->wedge.uwedge << "\t"
		   << this->position[0] << "\t"
		   << this->position[1] << "\t"
		   << this->position[2] << "\t"
		   << this->geometry[0] << "\t"
		   << this->geometry[1] << "\t"
		   << this->geometry[2] << "\t"
		   << this->normal[0] << "\t"
		   << this->normal[1] << "\t"
		   << this->normal[2] << "\t"
		   << this->matrix[0] << "\t"
		   << this->matrix[1] << "\t"
		   << this->matrix[2] << "\t"
		   << this->matrix[3] << "\t"
		   << this->matrix[4] << "\t"
		   << this->matrix[5] << "\t"
		   << this->matrix[6] << "\t"
		   << this->matrix[7] << "\t"
		   << this->matrix[8] << "\t"
		   << this->matrix[9] << "\t"
		   << this->matrix[10] << "\t"
		   << this->matrix[11] << "\t"
		   << this->matrix[12] << "\t"
		   << this->matrix[13] << "\t"
		   << this->matrix[14] << "\t"
		   << this->matrix[15] << "\t"
		   << this->magnification[0] << "\t"
		   << this->magnification[1] << "\t"
		   << this->magnification[2] << "\t"
		   << this->cutOffset << "\t"
		   << this->fileName << endl;
}

template< class Pixel >
void nbfWedgedSubImage3D< Pixel > :: unserialize( stringstream & input )
{
	Pixel lwedge, uwedge;
	input >> lwedge;
	input >> uwedge;
	this->wedge.set( lwedge, uwedge );
	input >> this->position[0];
	input >> this->position[1];
	input >> this->position[2];
	input >> this->geometry[0];
	input >> this->geometry[1];
	input >> this->geometry[2];
	input >> this->normal[0];
	input >> this->normal[1];
	input >> this->normal[2];
	input >> this->matrix[0];
	input >> this->matrix[1];
	input >> this->matrix[2];
	input >> this->matrix[3];
	input >> this->matrix[4];
	input >> this->matrix[5];
	input >> this->matrix[6];
	input >> this->matrix[7];
	input >> this->matrix[8];
	input >> this->matrix[9];
	input >> this->matrix[10];
	input >> this->matrix[11];
	input >> this->matrix[12];
	input >> this->matrix[13];
	input >> this->matrix[14];
	input >> this->matrix[15];
	input >> this->magnification[0];
	input >> this->magnification[1];
	input >> this->magnification[2];
	input >> this->cutOffset;
	input >> this->fileName;
}

template< class Pixel >
void nbfWedgedSubImage3D< Pixel > :: write( string file, vector< nbfWedgedSubImage3D< Pixel > > & vector ){
	
	// open file
	ofstream output( file.c_str(), ios::out );
	if ( output.is_open() != 1 ){
		output.close();
		return;		
	}
	
	// write header
	stringstream outputStream;
	nbfWedgedSubImage3D< Pixel > :: writeHeader( outputStream );
	
	// write vector elements
	for ( unsigned int i = 0; i < vector.size(); i++ ){
		outputStream << i+1 << "\t";
		vector[i].serialize( outputStream );

		// save wedge images if neccesary
		if ( vector[i].wedge.wedge.size() > 0 ){
			stringstream wedgeFile;
			wedgeFile <<  vector[i].fileName << ".wedge.mrc";
			nbfMrcWriter w;
			w.setFileName( wedgeFile.str().c_str() );
			vtkImageData * data = vtkImageData :: New();
			nbfVTKInterface :: blitzToVtk( vector[i].wedge.wedge, data );
			w.write( data );
			data->Delete();
		}
		if ( vector[i].wedge.sphericalWedge.size() > 0 ){
			stringstream sphericalWedgeFile;
			sphericalWedgeFile <<  vector[i].fileName << ".spherical.wedge.mrc";
			nbfMrcWriter w;
			w.setFileName( sphericalWedgeFile.str().c_str() );
			Array< Pixel, 3 > sW( vector[i].wedge.sphericalWedge.rows(), vector[i].wedge.sphericalWedge.cols(), 1 );
			sW( Range :: all(), Range :: all(), 0 ) = vector[i].wedge.sphericalWedge;
			vtkImageData * data = vtkImageData :: New();
			nbfVTKInterface :: blitzToVtk( sW, data );
			w.write( data );
			data->Delete();
		}
	}

	output << outputStream.str();

	// close file
	output.close();
}

template< class Pixel >
void nbfWedgedSubImage3D< Pixel > :: read( string file, vector< nbfWedgedSubImage3D< Pixel > > & vector )
{
	// open file
	ifstream input( file.c_str(), ios::in );
	if ( input.is_open() != 1 ){
		input.close();
		return;		
	}

	// get file size
	input.seekg( 0, ios::end );
	int lenght = input.tellg();
	input.seekg( 0, ios::beg );

	char * buffer;
	buffer = new char [lenght];
	input.read (buffer,lenght);

	stringstream inputStream;
	inputStream << buffer;

	// read and discard header
	nbfWedgedSubImage3D< Pixel > :: readHeader( inputStream );

	int count = 0;
	
	vector.clear();

#if 0
	// BACK COMPATIBILITY
	while ( inputStream.tellg() < lenght - 2 - count ){ // all 'endl's and 'eof' characters
		nbfWedgedSubImage3D< Pixel > volume;
		inputStream >> count;
		volume.unserialize( inputStream );
		vector.push_back( volume );
	}
#else

	vtkImageData * data = vtkImageData :: New();

	while ( true ){
		nbfWedgedSubImage3D< Pixel > volume;
		inputStream >> count;
		volume.unserialize( inputStream );

		if ( inputStream.fail() ){
			break;
		}

		// attempt to read wedge image
		if ( abs( volume.wedge.lwedge ) + abs( volume.wedge.uwedge ) == 0 ){
			stringstream wedgeFile;
			wedgeFile <<  volume.fileName << ".wedge.mrc";
			nbfMrcReader r;
			r.setFileName( wedgeFile.str().c_str() );
			data->Initialize();
			r.read( data );
			if ( data->GetNumberOfPoints() > 0 ){
				Array< Pixel, 3 > W;
				nbfVTKInterface :: vtkToBlitzReference( data, W );
				volume.wedge.setImage( W );
			}

			// attempt to read spherical wedge image
			stringstream sphericalWedgeFile;
			sphericalWedgeFile <<  volume.fileName << ".spherical.wedge.mrc";
			r.setFileName( sphericalWedgeFile.str().c_str() );
			data->Initialize();
			r.read( data );
			if ( data->GetNumberOfPoints() > 0 ){
				Array< Pixel, 3 > W;
				nbfVTKInterface :: vtkToBlitzReference( data, W );
				Array< Pixel, 2 > myW( W( Range :: all(), Range :: all(), 0 ) );
				volume.wedge.setSphericalImage( myW );
			}
		}

		vector.push_back( volume );
	}
#endif

	data->Delete();

	delete [] buffer;

	// close file
	input.close();
}

template< class Pixel >
void nbfWedgedSubImage3D< Pixel > :: getImage( vtkImageData * volume, vtkTransform * t, bool normalize )
{
	// if first time or if file changed: read volume from file
	if ( ( this->tomogram->GetNumberOfPoints() == 0 ) || ( strcmp( this->staticFileName.c_str(), this->fileName.c_str() ) != 0 ) ){
		if ( strcmp( this->fileName.c_str(), "NA" ) != 0 ){
			nbfMrcReader reader;
			reader.setFileName( this->fileName.c_str() );
			reader.read( this->tomogram );
			//Array< short, 3 > T;
			//nbfVTKInterface::vtkToBlitzReference( this->tomogram, T );
			//cout << "File " << this->fileName.c_str() << " read succesful, size = " << T.shape() << endl;
			//nbfMrcWriter w;
			//w.setFileName("mrc.mrc");
			//w.write( this->tomogram );

			// normalize to zero mean and unit variance
			if ( normalize == true ){
				Array< float, 3 > T;
				if ( this->tomogram->GetScalarType() == VTK_FLOAT ){
					nbfVTKInterface :: vtkToBlitzReference( this->tomogram, T );
					Array< float, 3 > T_bin4( T( Range(fromStart, toEnd, 4), Range(fromStart, toEnd, 4), Range(fromStart, toEnd, 4) ) );
					Pixel meanInside = mean( T_bin4 );
					Pixel varianceInside = sqrt( sum( pow2( T_bin4 - meanInside ) ) / T_bin4.size() );
					T = ( T - meanInside ) / varianceInside;
					cout << "NOTICE - Volume " << this->fileName.c_str() << " is being normalized to zero mean and unit variance before cutting." << endl;
				} 
				//else {
				//	if ( this->tomogram->GetScalarType() == VTK_SHORT ){
				//		Array< short, 3 > Ts;
				//		nbfVTKInterface :: vtkToBlitzReference( this->tomogram, Ts );
				//		Array< short, 3 > T_bin4( Ts( Range(fromStart, toEnd, 4), Range(fromStart, toEnd, 4), Range(fromStart, toEnd, 4) ) );
				//		Pixel meanInside = mean( T_bin4 );
				//		Pixel varianceInside = sqrt( sum( pow2( T_bin4 - meanInside ) ) / T_bin4.size() );
				//		Ts = ( Ts - meanInside ) / varianceInside;
				//		// nbfVTKInterface :: blitzToVtk( Ts, this->tomogram );
				//	}
				//}
				//nbfMatlabWriter w;
				//w.setFileName("p.matlab");
				//w.write( T );
			}
		}
		this->staticFileName = this->fileName;
	}

	int dims[3];
	if ( this->fixedImage != NULL ){
		this->fixedImage->GetDimensions(dims);
	} else {
		this->tomogram->GetDimensions(dims);
	}

	int currentDataType;

	// if volume size equal to file size -> assume we are dealing with full volumes
	//if ( ( dims[0] == this->geometry[0] ) && ( dims[1] == this->geometry[1] ) && ( dims[2] == this->geometry[2] ) ){
	if ( this->fixedImage != NULL ){
		//if ( this->fixedImage != NULL ){
		//	this->tomogram->DeepCopy( this->fixedImage );
		//}
		////if ( this->fixedImage == NULL ){
		////	this->fixedImage = vtkImageData::New();
		////}
		//////if ( this->fixedImage->GetNumberOfPoints() == 0 ){
		////	this->fixedImage->DeepCopy( this->tomogram );
		//////}
		//int dims[3];
		////this->fixedImage->GetDimensions(dims);
		//this->tomogram->GetDimensions(dims);
		//TinyVector< int, 3 > D( dims[0], dims[1], dims[2] );
		//this->setDimensions( D );
		////TinyVector< Pixel, 3 > Dn = ( D - 1.0 ) / 2.0;
		//TinyVector< Pixel, 3 > Dn = ( D - 0.0 ) / 2.0;
		//this->setPosition( Dn );
		//TinyVector< Pixel, 3 > T( 0, 0, 1 );
		////this->setNormal( T );
		////this->setCutOffset( ( dims[2] - 1.0 ) / 2.0 );
		//this->setCutOffset( ( dims[2] - 0.0 ) / 2.0 );
		////this->change->SetInput( this->fixedImage );
		////currentDataType = this->fixedImage->GetScalarType();
		this->change->SetInput( this->fixedImage );
		currentDataType = this->fixedImage->GetScalarType();
	}
	else {
		this->change->SetInput( this->tomogram );
		currentDataType = this->tomogram->GetScalarType();
	}
	
	Pixel lowerBoundForType;
	if ( true ){
		switch ( currentDataType ){
		case VTK_INT:
			lowerBoundForType = - numeric_limits< int > :: max();
			break;
		case VTK_SHORT:
			lowerBoundForType = - numeric_limits< short > :: max();
			break;
		case VTK_FLOAT:
			lowerBoundForType = - numeric_limits< float > :: max();
			break;
		case VTK_DOUBLE:
			lowerBoundForType = - numeric_limits< double > :: max();
			break;
		}
	} else {
		lowerBoundForType = 0.0;
	}
	this->reslice->SetBackgroundLevel( lowerBoundForType );

	this->change->SetOriginTranslation( - this->position[0], - this->position[1], - this->position[2] );
	this->change->Update();

	// set relative position of template wrt the image volume
#if 1
	this->box->SetDimensions( this->geometry[0], this->geometry[1], this->geometry[2] );
#else		
	if ( this->box->GetNumberOfPoints() != 0 ){
		this->box->SetDimensions( this->geometry[0], this->geometry[1], this->geometry[2] );
		this->box->AllocateScalars();
	}
#endif

	this->changeGeometry->SetOriginTranslation( - this->geometry[0] / 2.0, - this->geometry[1] / 2.0, - this->geometry[2] / 2.0 - this->cutOffset );
	this->changeGeometry->Update();

	vtkTransform * euler = vtkTransform::New();

	// compute rotation angles from normal orientation
	//Pixel angles[3];
	//angles[0] = 0;
	//angles[1] = atan2( sqrt( pow2(this->normal[0]) + pow2(this->normal[1]) ), this->normal[2] ) * vtkMath::RadiansToDegrees();
	//angles[2] = - atan2( -this->normal[1], this->normal[0] ) * vtkMath::RadiansToDegrees();

	//// THE ORDER OF THE ROTATIONS IS IMPORTANT!! //
	//euler->RotateZ( angles[2] );
	//euler->RotateY( angles[1] );
	//euler->RotateX( angles[0] );

	euler->RotateZ( - this->normal[2] );
	euler->RotateX( - this->normal[0] );
	euler->RotateZ( - this->normal[1] );

	// transform1->Concatenate( this->matrix );
	vtkMatrix4x4 * local = vtkMatrix4x4 :: New();
	vtkMatrix4x4 * current = vtkMatrix4x4 :: New();
	current->DeepCopy( this->matrix );
	vtkMatrix4x4::Multiply4x4( euler->GetMatrix(), current, local );
	euler->Delete();
	current->Delete();

	vtkMatrix4x4 * final = vtkMatrix4x4 :: New();
	if ( t != NULL ){
		vtkMatrix4x4::Multiply4x4( local, t->GetMatrix(), final );
		//vtkTransform * test = vtkTransform::New();
		//test->SetInput( t );
		//transform1->Concatenate( test );
		//test->Delete();
	} else {
		final->DeepCopy( local );
	}
	local->Delete();

	// check if final transform is valid
	if ( final->Determinant() == 0 ){
		cerr << "In file " << __FILE__ << ", line " << __LINE__ << "\nWARNING - Transformed image is out of range. Check that transformations are properly set.\nIn trying to transform file " << this->fileName << "." << endl;
		cerr << *final << endl;
		final->Identity();
		//vtkTransform * tmp = vtkTransform::New();
		//transform1->DeepCopy( tmp );
		//tmp->Delete();
	}

	vtkTransform * tfinal = vtkTransform :: New();
	tfinal->SetMatrix( final );
	final->Delete();
	this->reslice->SetResliceTransform( tfinal );
	this->reslice->Update();
	this->cast->Update();

	//cout << *tfinal->GetMatrix() << endl;

	//// normalize to zero mean and unit variance
	//if ( normalize ){
	//	Array< double, 3 > T;
	//	nbfVTKInterface::vtkToBlitzReference( this->cast->GetOutput(), T );
	//	Pixel normalization = sum( where( T > lowerBoundForType, 1, 0 ) );
	//	Pixel meanInside = sum( where( T > lowerBoundForType, T, 0 ) ) / normalization;
	//	Pixel varianceInside = sqrt( sum( where( T > lowerBoundForType, pow2( T - meanInside ), 0 ) ) / normalization );
	//	T = where( T > lowerBoundForType, ( T - meanInside ) / varianceInside, 0 );
	//}

	Array< double, 3 > T;
	nbfVTKInterface :: vtkToBlitzReference( this->cast->GetOutput(), T );
	Pixel normalization = sum( where( T > lowerBoundForType, 1, 0 ) );
	Pixel meanInside = sum( where( T > lowerBoundForType, T, 0 ) ) / normalization;
	if ( normalize == false ){
		T = where( T > lowerBoundForType, T - meanInside, 0 );
	} else {
		T = where( T > lowerBoundForType, T, meanInside );
	}

	//if ( ( this->magnification[0] != 1 ) || ( this->magnification[1] != 1 ) || ( this->magnification[2] != 1 ) ){
	//	this->resample->SetAxisMagnificationFactor( 0, this->magnification[0] );
	//	this->resample->SetAxisMagnificationFactor( 1, this->magnification[1] );
	//	this->resample->SetAxisMagnificationFactor( 2, this->magnification[2] );
	//	this->resample->Update();
	//	volume->ShallowCopy( this->resample->GetOutput() );
	//} else {
		volume->ShallowCopy( this->cast->GetOutput() );
	//}
	
	tfinal->Delete();

	//nbfVTKInterface :: vtkToBlitzReference( volume, T );
#if 0
	//vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	//writer->SetFileName("vol.vtk");
	//writer->SetInput( volume );
	//writer->Write();
	//writer->Delete();
	nbfMatlabWriter w;
	w.setFileName("p.matlab");
	nbfVTKInterface :: vtkToBlitz( volume, T );
	w.write(T);
#endif
}

template< class Pixel >
void nbfWedgedSubImage3D< Pixel > :: getImage( Array< Pixel, 3 > & A, vtkTransform * t, bool normalize )
{
	vtkImageData * data = vtkImageData::New();
	this->getImage( data, t, normalize );
	nbfVTKInterface::vtkToBlitz( data, A );
	data->Delete();
}


template< class Pixel >
void nbfWedgedSubImage3D< Pixel > :: getWedgeImage(  Array< Pixel, 3 > & volume, vtkTransform * trans )
{
	// new with euler angles
	vtkTransform * transform = vtkTransform::New();
	transform->RotateZ( - this->normal[2] );
	transform->RotateX( - this->normal[0] );
	transform->RotateZ( - this->normal[1] );

	vtkMatrix4x4 * local = vtkMatrix4x4 :: New();
	vtkMatrix4x4 * current = vtkMatrix4x4 :: New();
	current->DeepCopy( this->matrix );
	vtkMatrix4x4::Multiply4x4( transform->GetMatrix(), current, local );
	transform->Delete();
	current->Delete();

	vtkMatrix4x4 * final = vtkMatrix4x4 :: New();
	if ( trans != NULL ){
		vtkMatrix4x4::Multiply4x4( local, trans->GetMatrix(), final );
	} else {
		final->DeepCopy( local );
	}
	local->Delete();

	final->Invert();

	vtkTransform * tfinal = vtkTransform :: New();
	tfinal->SetMatrix( final );
	final->Delete();

	this->wedge.getImage( volume, tfinal );
	tfinal->Delete();

#if 0
	// compute rotation angles from normal orientation
	Pixel angles[3];
	angles[0] = 0;
	angles[1] = atan2( sqrt( pow2(this->normal[0]) + pow2(this->normal[1]) ), this->normal[2] ) * vtkMath::RadiansToDegrees();
	angles[2] = - atan2( -this->normal[1], this->normal[0] ) * vtkMath::RadiansToDegrees();

	//Pixel offset[3];
	//if ( trans != NULL ){
	//	trans->GetOrientation(offset);
	//} else{
	//	offset[0] = offset[1] = offset[2] = 0;
	//}

	// THE ORDER OF THE ROTATIONS IS IMPORTANT!! //
	vtkTransform * transform = vtkTransform::New();
	transform->RotateX( - angles[0] );
	transform->RotateY( - angles[1] );
	transform->RotateZ( - angles[2] );
	//transform->RotateX( - angles[0] - this->rotation[0] - offset[0] );
	//transform->RotateY( - angles[1] - this->rotation[1] - offset[1] );
	//transform->RotateZ( - angles[2] - this->rotation[2] - offset[2] );
	//transform->Concatenate( this->matrix );

	//// WORKS
	//vtkTransform * transt = vtkTransform::New();
	////transt->SetMatrix( this->matrix );
	//transt->Concatenate( this->matrix );
	//transt->Inverse();
	//transform->Concatenate( transt );
	//transt->Delete();

	// EXPERIMENTAL
	vtkTransform * transt = vtkTransform::New();
	transt->SetMatrix( this->matrix );
	transt->Inverse();
	transform->PostMultiply();
	transform->Concatenate( transt );
	transt->Delete();
	//if ( trans != NULL ){
	//	vtkTransform * invTrans = vtkTransform::New();
	//	//invTrans->SetMatrix( trans->GetMatrix() );
	//	invTrans->SetInput( trans );
	//	//invTrans->Inverse();
	//	transform->Concatenate( invTrans );
	//	invTrans->Delete();
	//}

	if ( trans != NULL ){
		vtkTransform * invTrans = vtkTransform::New();
		invTrans->SetInput( trans );
		invTrans->Inverse();
		transform->PostMultiply();
		transform->Concatenate( invTrans );
		invTrans->Delete();
	}

	//if ( trans != NULL ){
	//	transform->Concatenate( trans );
	//}

	//volume.resize( this->geometry );
	this->wedge.getImage( volume, transform );
	transform->Delete();
#endif
}

template< class Pixel >
void nbfWedgedSubImage3D< Pixel > :: getWedgeImageHalf(  Array< Pixel, 3 > & volume, vtkTransform * trans )
{
	// new with euler angles
	vtkTransform * transform = vtkTransform::New();
	transform->RotateZ( - this->normal[2] );
	transform->RotateX( - this->normal[0] );
	transform->RotateZ( - this->normal[1] );

	vtkMatrix4x4 * local = vtkMatrix4x4 :: New();
	vtkMatrix4x4 * current = vtkMatrix4x4 :: New();
	current->DeepCopy( this->matrix );
	vtkMatrix4x4::Multiply4x4( transform->GetMatrix(), current, local );
	transform->Delete();
	current->Delete();

	vtkMatrix4x4 * final = vtkMatrix4x4 :: New();
	if ( trans != NULL ){
		vtkMatrix4x4::Multiply4x4( local, trans->GetMatrix(), final );
	} else {
		final->DeepCopy( local );
	}
	local->Delete();

	final->Invert();

	vtkTransform * tfinal = vtkTransform :: New();
	tfinal->SetMatrix( final );
	final->Delete();

	this->wedge.getImageHalf( volume, tfinal );
	tfinal->Delete();
#if 0
	// compute rotation angles from normal orientation
	Pixel angles[3];
	angles[0] = 0;
	angles[1] = atan2( sqrt( pow2(this->normal[0]) + pow2(this->normal[1]) ), this->normal[2] ) * vtkMath::RadiansToDegrees();
	angles[2] = - atan2( -this->normal[1], this->normal[0] ) * vtkMath::RadiansToDegrees();

	vtkTransform * transform = vtkTransform::New();
	transform->RotateX( - angles[0] );
	transform->RotateY( - angles[1] );
	transform->RotateZ( - angles[2] );

	// EXPERIMENTAL
	vtkTransform * transt = vtkTransform::New();
	transt->SetMatrix( this->matrix );
	transt->Inverse();
	transform->PostMultiply();
	transform->Concatenate( transt );
	transt->Delete();

	if ( trans != NULL ){
		vtkTransform * invTrans = vtkTransform::New();
		invTrans->SetInput( trans );
		invTrans->Inverse();
		transform->PostMultiply();
		transform->Concatenate( invTrans );
		invTrans->Delete();
	}
	this->wedge.getImageHalf( volume, transform );
	transform->Delete();
#endif
}


template< class Pixel >
void nbfWedgedSubImage3D< Pixel > :: getSphericalWedgeImage(  Array< Pixel, 2 > & volume, vtkTransform * trans, nbfProjectionRotationMetric3D< Pixel > * metric )
{
	if ( this->wedge.sphericalWedge.size() > 0 ){
		this->wedge.getSphericalImage( volume, trans, metric );
	} else {
		// new with euler angles
		vtkTransform * transform = vtkTransform::New();
		transform->RotateZ( - this->normal[2] );
		transform->RotateX( - this->normal[0] );
		transform->RotateZ( - this->normal[1] );

		vtkMatrix4x4 * local = vtkMatrix4x4 :: New();
		vtkMatrix4x4 * current = vtkMatrix4x4 :: New();
		current->DeepCopy( this->matrix );
		vtkMatrix4x4::Multiply4x4( transform->GetMatrix(), current, local );
		transform->Delete();
		current->Delete();

		vtkMatrix4x4 * final = vtkMatrix4x4 :: New();
		if ( trans != NULL ){
			vtkMatrix4x4::Multiply4x4( local, trans->GetMatrix(), final );
		} else {
			final->DeepCopy( local );
		}
		local->Delete();

		final->Invert();

		vtkTransform * tfinal = vtkTransform :: New();
		tfinal->SetMatrix( final );
		final->Delete();

		this->wedge.getSphericalImage( volume, tfinal, metric );
		tfinal->Delete();
	}
#if 0
	// compute rotation angles from normal orientation
	Pixel angles[3];
	angles[0] = 0;
	angles[1] = atan2( sqrt( pow2(this->normal[0]) + pow2(this->normal[1]) ), this->normal[2] ) * vtkMath::RadiansToDegrees();
	angles[2] = - atan2( -this->normal[1], this->normal[0] ) * vtkMath::RadiansToDegrees();

	//Pixel offset[3];
	//if ( trans != NULL ){
	//	trans->GetOrientation(offset);
	//} else{
	//	offset[0] = offset[1] = offset[2] = 0;
	//}

	// THE ORDER OF THE ROTATIONS IS IMPORTANT!! //
	vtkTransform * transform = vtkTransform::New();

	//transform->RotateZ( this->rotation[2] );
	//transform->RotateY( this->rotation[1] );
	//transform->RotateX( this->rotation[0] );

	//transform->RotateX( - angles[0] - this->rotation[0] );
	//transform->RotateY( - angles[1] - this->rotation[1] );
	//transform->RotateZ( - angles[2] - this->rotation[2] );

	transform->RotateX( - angles[0] );
	transform->RotateY( - angles[1] );
	transform->RotateZ( - angles[2] );

	vtkTransform * transt = vtkTransform::New();
	//transt->SetMatrix( this->matrix );
	transt->Concatenate( this->matrix );
	transt->Inverse();
	transform->Concatenate( transt );
	transt->Delete();

	//if ( trans != NULL ){
	//	vtkTransform * invTrans = vtkTransform::New();
	//	//invTrans->SetMatrix( trans->GetMatrix() );
	//	invTrans->SetInput( trans );
	//	invTrans->Inverse();
	//	transform->Concatenate( invTrans->GetMatrix() );
	//	invTrans->Delete();
	//}

	if ( trans != NULL ){
		vtkTransform * invTrans = vtkTransform::New();
		invTrans->SetInput( trans );
		invTrans->Inverse();
		transform->PostMultiply();
		transform->Concatenate( invTrans->GetMatrix() );
		invTrans->Delete();
	}

	//transform->RotateZ( angles[2] );
	//transform->RotateY( angles[1] );
	//transform->RotateX( angles[0] );

	//transform->RotateZ( offset[2] );
	//transform->RotateY( offset[1] );
	//transform->RotateX( offset[0] );

	//transform->Concatenate( trans->GetMatrix() );


	//volume.resize( volume.rows(), volume.depth() );
	this->wedge.getSphericalImage( volume, transform );
	transform->Delete();
#endif
}