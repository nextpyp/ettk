#define BZ_GENERATE_GLOBAL_INSTANCES

#ifndef WIN32
	#include "mpi.h"
#endif

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <string.h>

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkImageReader.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageReslice.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkImageNonMaximumSuppression.h>
#include <vtkImageMathematics.h>
#include <vtkPolyData.h>
#include <vtkImageCast.h>
#include <vtkImageExtractComponents.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMODReader.h>
#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMrcWriter.h>

#include <em/nbfImageMetric.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfWedgedAverageImage3D.h>

#include <em/nbfExtractPointsAndNormals3D.h>

#define PIXEL double

#define USING_PHANTOM_LIST 0
#define LOOP_CREATE_VOLUME_LIST_UNDEF -1
#define LOOP_CREATE_VOLUME_LIST_POS 0
#define LOOP_CREATE_VOLUME_LIST_EUL 1
#define LOOP_CREATE_VOLUME_LIST_MRC 2
#define LOOP_CREATE_VOLUME_LIST_MOD 3

int main( int argc, char ** argv )
{
	cout << argv[0] << "\n";
	cout << " Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
#if USING_PHANTOM_LIST
	
	if ( argc == 1 ){
		cout << "USAGE: wedgeLimit ground_truth.matlab list_of_files volumeList.txt" << endl;
	}

	double wedge = atof(argv[1]);

	Array< PIXEL, 2 > ground_truth;
	nbfMatlabReader reader;
	reader.setFileName( argv[2] );
	reader.read( ground_truth );

	vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
	nbfWedgedSubImage3D< PIXEL > :: read( argv[argc-1], volumeList );

	for ( int tomograms = 0; tomograms < argc - 4; tomograms++ ){	
		stringstream tomogramFile;
		tomogramFile << argv[ tomograms + 3 ] << ".rec";
		nbfWedgedSubImage3D< PIXEL > currentVolume;
		currentVolume.getWedge()->set( -wedge, wedge );
		currentVolume.setFileName( tomogramFile.str().c_str() );
		//currentVolume.imageToFullVolumeOn();
		currentVolume.setDimensions( TinyVector< int, 3 >(64,64,64) );
		vtkImageData * data = vtkImageData::New();
		currentVolume.getImage(data);
		data->Delete();

		vtkTransform * t = vtkTransform::New();
		t->RotateZ( ground_truth( tomograms, 0 ) );
		t->RotateY( ground_truth( tomograms, 1 ) );
		t->RotateZ( ground_truth( tomograms, 2 ) );
		t->Inverse();
		currentVolume.setTransform( t );
		t->Delete();
		volumeList.push_back( currentVolume );
	}

	nbfWedgedSubImage3D< PIXEL > :: write( argv[argc-1], volumeList );
#else

	if ( argc < 10 ){
		cout << "USAGE:\nLoopCreatePhantomVolumeList magnification tilts tomogram positions surface levelset list\n" << endl;
		cout << "magnification - binning factor applied to get spike coordinates." << endl;
		cout << "lower tilt limit. e.g. -70." << endl;
		cout << "upper tilt limit. e.g. 70." << endl;
		cout << "tomogram - un-binned mrc file containing single virus." << endl;
		cout << "positions - file containing x,y,x spike locations." << endl;
		cout << "offset - along normal direction (in pixels)." << endl;
		cout << "surface - mrc file obtained from membrane segmentation." << endl;
		cout << "levelset - surface level corresponding to membrane." << endl;
		cout << "list - output txt file containing extracted spikes." << endl;
		return 0;
	}

	// PARSE PARAMETERS

	double mag = atof( argv[1] );
	
	PIXEL wedgel = atof( argv[2] );
	PIXEL wedgeu = atof( argv[3] );

	//PIXEL wedge = atof( argv[2] );

	//stringstream tiltSeriesFile;
	//tiltSeriesFile << argv[ 2 ];

	stringstream tomogramFile;
	tomogramFile << argv[ 4 ];

	stringstream positionsFile;
	positionsFile << argv[ 5 ];

	PIXEL offset = atof( argv[ 6 ] );

	stringstream surfaceFile;
	surfaceFile << argv[ 7 ];

	PIXEL surfaceLevelSet = atof( argv[8] );

	stringstream volumeListFile;
	volumeListFile << argv[ 9 ];

	// UNDOCUMENTED FEATURE
	PIXEL angle_th = 90;
	if ( argc > 10 ){
		angle_th = atof( argv[9] );
		cout << "**WARNING** UNDOCUMENTED USE OF COMMAND" << endl;
		cout << "            SKIPPING POSTIONS THAT EXCEED ANGULAR THRESHOLD." << endl;
	}

	//////////////////////////////
	nbfWedgedSubImage3D< PIXEL > cutVolume;
	//cutVolume.setDimensions( TinyVector<int,3>( volumeSize, volumeSize, volumeSize ) );
	//cutVolume.setCutOffset( volumeSize / 2.0 );
	//cutVolume.setCutOffset( 8.0 );
	//cutVolume.setCutOffset( 50.0 );
	//cutVolume.setMagnification( TinyVector< PIXEL, 3 >(2.0,2.0,2.0) );

	// MISSING WEDGE
	Array< float, 1 > angles;
	nbfMrcReader reader;
#if 0
	// read aligned tilt series to extract wedge geometry
	reader.setFileName( tiltSeriesFile.str().c_str() );
	Array< float, 1 > means;

	reader.read( angles, means );
	
	if ( angles.size() == 0 ){
		angles.resize(2);
		angles = numeric_limits< float > :: max();
		//cout << "Unable to get tilt range from tilt series " << argv[2] << ".\nPlease specify negative range: ";
		//cin >> angles(0);		
		//cout << "Please specify positive tilt range: ";
		//cin >> angles(1);
		angles(0) = - wedge;
		angles(1) = wedge;
	}
#else
	angles.resize(2);
	angles(0) = wedgel;
	angles(1) = wedgeu;
#endif

	// check angle consistency
	if ( angles(0) < -90 ){
		cout << "Invalid negative range, defaulting to -90." << endl;
		angles(0) = -90;
	}
	if ( angles(1) > 90 ){
		cout << "Invalid positive range, defaulting to 90." << endl;
		angles(1) = 90;
	}
	if ( angles(0) > angles(1) ){
		cout << "Invalid tilt range, negative limits cannot be bigger than positive limit." << endl;
		angles(0) = -90; angles(1) = 90;
	}
	cout << "Setting wedge geometry to [" << angles(0) << "," << angles(1) << "]." << endl;

	cutVolume.getWedge()->set( angles(0), angles(angles.ubound(firstDim)) );

	// PARTICLE SELECTION

	cout << "Reading positions from " << positionsFile.str() << ": ";

	int type = LOOP_CREATE_VOLUME_LIST_UNDEF;
	int lenght = positionsFile.str().size();
	if ( ( positionsFile.str()[ lenght - 4 ] == '.' ) &&
       ( positionsFile.str()[ lenght - 3 ] == 'p' ) &&
       ( positionsFile.str()[ lenght - 2 ] == 'o' ) &&
       ( positionsFile.str()[ lenght - 1 ] == 's' ) )
    {
		type = LOOP_CREATE_VOLUME_LIST_POS;
	}
	if ( ( positionsFile.str()[ lenght - 4 ] == '.' ) &&
       ( positionsFile.str()[ lenght - 3 ] == 'e' ) &&
       ( positionsFile.str()[ lenght - 2 ] == 'u' ) &&
       ( positionsFile.str()[ lenght - 1 ] == 'l' ) )
    {
		type = LOOP_CREATE_VOLUME_LIST_EUL;
	}
	if ( ( positionsFile.str()[ lenght - 4 ] == '.' ) &&
       ( positionsFile.str()[ lenght - 3 ] == 'm' ) &&
       ( positionsFile.str()[ lenght - 2 ] == 'o' ) &&
       ( positionsFile.str()[ lenght - 1 ] == 'd' ) )
    {
		type = LOOP_CREATE_VOLUME_LIST_MOD;
	}

	if ( type == LOOP_CREATE_VOLUME_LIST_UNDEF ){
		cout << "ERROR - a valid position file must be specified (.pos,.eul,.mod)." << endl;
	}

	nbfMODReader mr;
	mr.setFileName( positionsFile.str().c_str() );

	Array< float, 2 > points;

	switch ( type ){
			case LOOP_CREATE_VOLUME_LIST_POS:
				mr.readRaw( points );
				break;
			case LOOP_CREATE_VOLUME_LIST_EUL:
				mr.readRawEul( points );
				if ( points.cols() != 6 ){
					cout << "ERROR: eul file format not recognized." << endl;
					return 0;
				}
				break;
			case LOOP_CREATE_VOLUME_LIST_MOD:
				mr.read( points );
				// FIX READING ERRORS to match imod display
				points = where( points - floor(points) > .5, ceil(points), floor(points) );
				break;
	}

	if ( points.size() > 0 ){
		if ( ( type == LOOP_CREATE_VOLUME_LIST_POS ) || ( type == LOOP_CREATE_VOLUME_LIST_EUL ) ){
			cout << points.rows() << " locations to extract.\n";
		} else {
			cout << points.rows() << " locations to extract.\n";
		}
	}
	else{
		cout << "ERROR reading " << positionsFile.str() << ". No locations will be extracted." << endl;
		return 0;
	}

	// get dimensions of input tomogram so we can reverse second coordinate
	reader.setFileName( tomogramFile.str().c_str() );
	TinyVector< int, 3 > dims = reader.getDims();
	//int dimY = reader.getDimY();
	cout << "Extracting geometry from " << tomogramFile.str() << ":\n" << dims << endl;	

	// MAGNIFICATION: magnifiy original point coordinates to match full resolution image data
	// points( Range::all(), Range(0,2) ) = mag * ( points( Range::all(), Range(0,2) ) + 1 ) - mag / 2 - 1;
	points( Range::all(), Range(0,2) ) = ( mag * points( Range::all(), Range(0,2) ) );

	// EULER ANGLES FROM SURFACE
	
	vtkPolyData * vpoints = vtkPolyData :: New();
	if ( type == LOOP_CREATE_VOLUME_LIST_POS ){
		cout << "Extracting euler angles from membrane orientation" << endl;
		// read implicit surface
		reader.setFileName( surfaceFile.str().c_str() );
		vtkImageData * mySurface = vtkImageData :: New();
		reader.read( mySurface );
		nbfExtractPointsAndNormals3D< PIXEL > extract;
		extract.setSurface( mySurface, surfaceLevelSet );
		extract.setMagnification( mag );
		extract.execute( vpoints );
		mySurface->Delete();
	}


	// read volume list and append new data
	vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
	nbfWedgedSubImage3D< PIXEL > :: read( volumeListFile.str().c_str(), volumeList );

	// get coordinates and euler angles
	for ( int i = 0; i < points.rows(); i++ ){

		nbfWedgedSubImage3D< PIXEL > currentVolume;
		currentVolume = cutVolume;
		currentVolume.setFileName( tomogramFile.str().c_str() );

		TinyVector< PIXEL, 3 > pos( points(i,0), dims[1] - 1 - points(i,1), points(i,2) );
		//TinyVector< PIXEL, 3 > pos( points(i,0), points(i,1), points(i,2) );
		TinyVector< PIXEL, 3 > normal(0,0,0);

		// skip if not on side surface of virus
		TinyVector< PIXEL, 3 > spike( pos - dims / 2 );
		spike = spike / sqrt( sum( spike * spike ) );
		//cout << spike << endl;
		TinyVector< PIXEL, 3 > zdirection( 0, 0, 1 );
		// PIXEL angle = acos( dot( zdirection, spike ) ) * vtkMath::RadiansToDegrees() - 90;
		PIXEL angle = vtkMath::DegreesFromRadians( acos( dot( zdirection, spike ) ) ) - 90;
		if ( fabs( angle ) > angle_th ){
			cout << "WARNING - position " << i << " skipped. Angle " << fabs(angle) << " is greater than " << angle_th << endl;
			continue;
		}

		// euler angles
		switch ( type ){
				case LOOP_CREATE_VOLUME_LIST_EUL:
					normal = TinyVector< PIXEL, 3 >( points(i,4), points(i,5), points(i,3) );
					break;
				case LOOP_CREATE_VOLUME_LIST_POS:
				case LOOP_CREATE_VOLUME_LIST_MOD:
					// normal orientation specified by pair of points
					// TinyVector< PIXEL, 3 > normal( points(i,0) - points(0,0), - points(i,1) + points(0,1), points(i,2) - points(0,2) ) );
					// normal orientation from center of virus
					// TinyVector< PIXEL, 3 > normal( points(i,0) - dims[0] / 2.0, - points(i,1) + dims[1] / 2.0, points(i,2) - dims[2] / 2.0 );
					// read implicit surface
					// take specified level set

					if ( vpoints->GetNumberOfPoints() > 0 ){
						int id = vpoints->FindPoint( pos(0), pos(1), pos(2) );
						//int id = vpoints->FindPoint( points(i,0), points(i,1), points(i,2) );
						PIXEL p[3];
						vpoints->GetPoint(id,p);
						//cout << "searching position = " << pos(0) << ", " << pos(1) << ", " << pos(2) << endl;
						//cout << "closest position = " << p[0] << ", " << p[1] << ", " << p[2] << endl;
						PIXEL tmp[3];
						vpoints->GetPointData()->GetNormals()->GetTuple( id, tmp );

		TinyVector< PIXEL, 3 > n( tmp[0], tmp[1], tmp[2] );
		pos = floor( pos - offset * n );


						//tmp[1]*=-1;
						// convert normal orientation to Euler angles
						// cout << "normal = " << tmp[0] << ", " << tmp[1] << ", " << tmp[2] << endl;
						// PIXEL rotX = atan( fabs( tmp[0] / tmp[1] ) ) * vtkMath :: RadiansToDegrees();
						PIXEL rotX = vtkMath::DegreesFromRadians( atan( fabs( tmp[0] / tmp[1] ) ) );
						if ( ( tmp[0] < 0 ) && ( tmp[1] > 0 ) ){
							rotX = - rotX;
						}
						if ( ( tmp[0] < 0 ) && ( tmp[1] < 0 ) ){
							rotX = - ( 180 - rotX );
						}
						if ( ( tmp[0] > 0 ) && ( tmp[1] < 0 ) ){
							rotX = 180 - rotX;
						}
                        // PIXEL rotZ = atan( sqrt( tmp[0] * tmp[0] + tmp[1] * tmp[1] ) / fabs( tmp[2] ) ) * vtkMath :: RadiansToDegrees();
						PIXEL rotZ = vtkMath::DegreesFromRadians( atan( sqrt( tmp[0] * tmp[0] + tmp[1] * tmp[1] ) / fabs( tmp[2] ) ) );
						if ( ( tmp[1] < 0 ) && ( tmp[2] > 0 ) ){
							rotZ = - rotZ;
						}
						if ( ( tmp[1] < 0 ) && ( tmp[2] < 0 ) ){
							rotZ = - ( 180 - rotZ );
						}
						if ( ( tmp[1] > 0 ) && ( tmp[2] < 0 ) ){
							rotZ = 180 - rotZ;
						}
						rotZ = fabs(rotZ);
						normal[0] = rotZ;
						normal[1] = 0;
						normal[2] = rotX;
						break;
					}
		}

		//pos(1) = dims[1] - 1 - pos(1);
		currentVolume.setPosition( pos );
		currentVolume.setNormal( normal );			
		volumeList.push_back( currentVolume );
	}

	vpoints->Delete();

	cout << "Adding " << points.rows() << " spikes to " << volumeListFile.str().c_str() << endl;
	nbfWedgedSubImage3D< PIXEL > :: write( volumeListFile.str().c_str(), volumeList );
	
	//// cut-out enlarged individual volumes
	//vector< nbfWedgedSubImage3D< PIXEL > > cutVolumeList;
	//
	//for ( int i = 0; i < volumeList.size(); i++ ){

	//	// cut-out enlarged sub volume
	//	TinyVector< int, 3 > tsize( volumeSize, volumeSize, volumeSize );
	//	volumeList[i].setDimensions( tsize );
	//	volumeList[i].setCutOffset( volumeSize / 2.0 );
	//	
	//	vtkImageData * data = vtkImageData::New();
	//	volumeList[i].getImage(data);
	//	nbfMrcWriter writer;
	//	stringstream subTomogramFile;
	//	subTomogramFile << argv[ 3 ] << "_";
	//	if ( i+1 < 10 ){
	//		subTomogramFile << "000";
	//	} else if ( i+1 < 100 ){
	//		subTomogramFile << "00";
	//	} else if ( i+1 < 1000 ){
	//		subTomogramFile << "0";
	//	}
	//	subTomogramFile << i+1 << ".rec";
	//	writer.setFileName( subTomogramFile.str().c_str() );

	//	vtkImageCast * cast = vtkImageCast::New();
	//	cast->SetOutputScalarTypeToShort();
	//	cast->SetInput(data);
	//	cast->Update();

	//	writer.write( cast->GetOutput() );
	//	data->Delete();
	//	cast->Delete();

	//	volumeList[i].setPosition( TinyVector< PIXEL, 3 >( volumeSize / 2, volumeSize / 2, volumeSize / 2 ) );
	//	volumeList[i].setNormal( TinyVector< PIXEL, 3 >( 0,0,1 ) );
	//	volumeList[i].setFileName( subTomogramFile.str().c_str() );

	//	cutVolumeList.push_back( volumeList[i] );
	//}

	//stringstream cutVolumeListName;
	//cutVolumeListName << argv[argc-1] << ".cut.txt";
	//nbfWedgedSubImage3D< PIXEL > :: write( cutVolumeListName.str().c_str(), cutVolumeList );

#endif
	return 0;
}
