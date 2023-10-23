#define NBF_VERBOSE 1

#include "mpi.h"

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
#include <vtkImageContinuousDilate3D.h>
#include <vtkImageCast.h>
#include <vtkImageExtractComponents.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMrcWriter.h>
#include <bs/nbfBordStrategyMirror.h>

#include <em/nbfImageMetric.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfWedgedAverageImage3D.h>

#include <em/nbfLoopClustering.h>

#define PIXEL float

int main( int argc, char ** argv )
{
	//  Initialize MPI.
	MPI_Init ( &argc, &argv );

	int my_id, num_procs;

	//  Get this processes's rank.
	MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

	//  Find out how many processes are available.
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	cout << "Process " << my_id << " is active.\n";

	nbfImageFilter< PIXEL, 3 > imfilter;
	imfilter.paddingOn(1);

	nbfFourierFilter< PIXEL, 3 > fffilter;
	//fffilter.bandPassOn(.01,.25,.001,.05); // groel
	fffilter.bandPassOn(.01,.1,.001,.05); // groel

	nbfProjectionRotationMetric3D< PIXEL > metric( &imfilter, &fffilter );
	metric.setNumberOfCandidatePeaksToSearch(10);

	if ( my_id == 0 ){ // Master process

		char * inputFileName = argv[1];
		char * outputFileName = argv[2];

		vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
		nbfWedgedSubImage3D< PIXEL > :: read( inputFileName, volumeList );

		if ( volumeList.size() < 2 ){
			cout << "ERROR - Volume list must contain at least two volumes." << endl;
		}

		// set volume size and cutting offset to center of subvolume
		int volumeSizeX = 256;
		int volumeSizeY = 256;
		int volumeSizeZ = 256;

		for ( int i = 0; i < volumeList.size(); i++ ){
			TinyVector< int, 3 > tsize( volumeSizeX, volumeSizeY, volumeSizeZ );
			volumeList[i].setDimensions( tsize );
			volumeList[i].setCutOffset( volumeSizeZ / 2.0 );
		}

		// reset original alignments
		for ( int i = 0; i < volumeList.size(); i++ ){
			volumeList[i].setTransform( (vtkTransform*)NULL );
		}

		Array< PIXEL, 3 > alignments;
		char * alignmentsFileName = argv[3];
		nbfMatlabReader mreader;
		mreader.setFileName( alignmentsFileName );
		mreader.read( alignments );

		if ( alignments.size() > 0 ){
			cout << "Read alignments " << alignments.shape() << endl;
			// make correct distance matrix in case file was corrupted
			metric.makeDistanceMatrix(alignments);
		} else {
			alignments.resize( volumeList.size(), volumeList.size(), 19 );
			alignments = -1;
			metric.getDistances( volumeList, alignments );
			nbfMatlabWriter mwriter;
			mwriter.setFileName( alignmentsFileName );
			mwriter.write( alignments );
		}

		nbfWedgedAverageImage3D< PIXEL > av;
		double matrix[16];
		for ( int k = 0; k < 16; k++ ){
			matrix[k] = alignments( 0, 1, k + 3 );
		}
		vtkMatrix4x4 * mat = vtkMatrix4x4::New();
		mat->DeepCopy( matrix );
		volumeList[ 1 ].setTransform( mat );
		mat->Delete();

		av.getVolumes().push_back( volumeList[0] );
		av.getWeights().push_back( 1.0 );

		av.getVolumes().push_back( volumeList[1] );
		av.getWeights().push_back( 1.0 );

		vtkImageData * output = vtkImageData::New();
		av.getImage(output);

		// change geometry to make it compatible with mrc file
		Array< double, 3 > A, B;
		nbfVTKInterface::vtkToBlitzReference( output, A );

		A.transposeSelf(thirdDim,secondDim,firstDim);
		B.resize( A.shape() );
		B = A.reverse(secondDim);
		nbfVTKInterface::blitzToVtk(B,output);

		vtkImageCast * cast = vtkImageCast::New();
		cast->SetOutputScalarTypeToFloat();
		cast->SetInput( output );
		cast->Update();

		nbfMrcWriter mrcw;
		mrcw.setFileName( outputFileName );
		mrcw.write( cast->GetOutput() );
		cast->Delete();
		output->Delete();

		metric.finalizeMPI();

	} else {
		metric.slaveMPI();
	}

	MPI_Finalize();

	return 0;
}