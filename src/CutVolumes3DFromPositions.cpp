#define BZ_GENERATE_GLOBAL_INSTANCES

#ifndef WIN32
	#include "mpi.h"
#endif

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

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
#include <vtkImageContinuousDilate3D.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfBlitzWriter.h>
#include <io/nbfBlitzReader.h>
#include <io/nbfMrcWriter.h>
#include <io/nbfMrcReader.h>
//#include <nbfVeselnessFilter.h>
#include <bs/nbfBordStrategyMirror.h>

#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>

#define PIXEL float

int main( int argc, char ** argv )
{
	cout << argv[0] << "\n";
	cout << " Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";

	//nbfMrcReader reader;
	//reader.setFileName( argv[1] );
	//vtkImageData * data = vtkImageData::New();
	//Array< float, 3 > A;
	//reader.read(data);
	//nbfVTKInterface::vtkToBlitzReference(data,A);

	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");
	//w.write(A);

	// argc[1] - input volume list
	// argc[2] - sub volume size
	// argc[3] - output file prefix
	// argc[4] - use unit normal
	// argc[5] - start position

	if ( argc < 5 ){
		cout << "USAGE:\n";
		cout << "CutVolumes3DFromPositions volumes cutsize prefix noeulers start\n";
		cout << "volumes - txt file obtained from LoopCreateVolumeList." << endl;
		cout << "cutsize - crop size of output volumes." << endl;
		cout << "prefix - string to preceed all volumes filename." << endl;
		cout << "noeulers - USE 1 (cropping done without using euler angles). 0 if cropping done using euler angles." << endl;
		cout << "start - index within the volumes list to start processing. Set to -1 to produce the output txt file alone (without generating cropped volumes)." << endl;
		return 0;
	}

	vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
	nbfWedgedSubImage3D< PIXEL > :: read( argv[1], volumeList );

	cout << volumeList.size() << " volumes to extract." << endl;

	// specify cropped region size
	int volumeSize = atoi(argv[2]);

	cout << "Volume size set to " << volumeSize << "." << endl;

	int unit_normal = atoi(argv[4]);
	int start_at = 0;
	if ( argc == 6 ){
		start_at = atoi(argv[5]);
	}
	bool alreadyCutVolumes = false;
	if ( start_at < 0 ){
		start_at = 0;
		alreadyCutVolumes = true;
	}
	cout << "Starting at position " << start_at << "." << endl;

	vtkImageData * vol = vtkImageData::New();

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToFloat();

	nbfMrcWriter writerm;

	vector< nbfWedgedSubImage3D< PIXEL > > cutVolumeList;
	stringstream cutVolumeListName;
	cutVolumeListName << argv[3] << ".txt";

	nbfWedgedSubImage3D< PIXEL > :: read( cutVolumeListName.str().c_str(), cutVolumeList );

	for ( int i = start_at; i < volumeList.size(); i++ ){

		cout << "Cutting volume " << i+1 << " at " << volumeList[i].getPosition() << endl;

		// set size and cutting offset
		TinyVector< int, 3 > tsize( volumeSize, volumeSize, volumeSize );
		volumeList[i].setDimensions( tsize );

		// save current offset
		PIXEL offset = volumeList[i].getCutOffset();
		volumeList[i].setCutOffset( 0.0 );

		// store normal
		TinyVector< PIXEL, 3 > normal = volumeList[i].getNormal();
		TinyVector< PIXEL, 3 > unitNormal( 0,0,0 );
		if ( unit_normal == 1 ){
			volumeList[i].setNormal( unitNormal );
		} 

		if ( alreadyCutVolumes == false ){
			nbfImageFilter< PIXEL, 3 > imf;
			volumeList[i].getImage(vol,(vtkTransform*)NULL,true); // virus-wise normalization

			//// fix geometry for MRC writer
			//Array< float, 3 > A, B;
			//nbfVTKInterface::vtkToBlitz( vol, A );
			//A.transposeSelf(thirdDim,secondDim,firstDim);
			//B.resize( A.shape() );
			//B = A.reverse(secondDim);
			//nbfVTKInterface::blitzToVtk(B,vol);
		}

		cast->SetInput( vol );
		cast->Update();

		stringstream fileName;
		if ( i < 9 ){
			fileName << argv[3] << "_000" << i+1 << ".mrc";
		}
		else if ( i < 99 ){
			fileName << argv[3] << "_00" << i+1 << ".mrc";
		}
		else if ( i < 999 ){
			fileName << argv[3] << "_0" << i+1 << ".mrc";
		}
		else{
			fileName << argv[3] << "_" << i+1 << ".mrc";
		}

		if ( alreadyCutVolumes == false ){
			writerm.setFileName( fileName.str().c_str() );
			writerm.write( cast->GetOutput(), true );
			//writerm.write( vol );
			cout << "File: " << fileName.str() << " written." << endl;
		}

		TinyVector< PIXEL, 3 > pos( volumeSize / 2, volumeSize / 2, volumeSize / 2 );
		volumeList[i].setPosition( pos );
	
		// always keep same normal
		// if ( unit_normal == 0 ){
			volumeList[i].setNormal( normal );
		//} else {
		//	TinyVector< PIXEL, 3 > unit( 0,0,0 );
		//	volumeList[i].setNormal( unit );
		//}
		volumeList[i].setFileName( fileName.str().c_str() );

		cutVolumeList.push_back( volumeList[i] );

		// reset offset to 0
		cutVolumeList.back().setCutOffset( offset );
	}	

	nbfWedgedSubImage3D< PIXEL > :: write( cutVolumeListName.str().c_str(), cutVolumeList );

	vol->Delete();
	cast->Delete();
}