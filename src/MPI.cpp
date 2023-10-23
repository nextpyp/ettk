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
#include <bs/nbfBordStrategyMirror.h>

#include <em/nbfImageMetric.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfWedgedAverageImage3D.h>

#define PIXEL double

#define BUILD_VOLUME_LIST 0

void main( int argc, char ** argv )
{
	//  Initialize MPI.
	MPI::Init ( argc, argv );

	int my_id, num_procs;
	int source;
	int tag;
	int tag_target = 1;
	int tag_size = 2;
	int tag_data = 3;
	int tag_found = 4;
	int tag_done = 5;
	MPI::Status status;

	//  Get this processes's rank.
	MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

	//  Find out how many processes are available.
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	cout << "Process " << my_id << " is active.\n";

	// the metric is the only think all processes share

	nbfImageFilter< PIXEL, 3 > imfilter;

	nbfFourierFilter< PIXEL, 3 > fffilter;
	fffilter.bandPassOn(.05,.1,.01);
	
	nbfProjectionRotationMetric3D< PIXEL > metric( &imfilter, &fffilter );
	metric.setNumberOfCandidatePeaksToSearch(1);
	
	if ( my_id == 0 ){

		cout << "MPI - Master process:\n";
		cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
		cout << "  The number of processes is " << num_procs << "\n";

		vector< nbfWedgedSubImage3D< PIXEL > > volumeList;

		stringstream fileName;
		fileName << argv[1];

		nbfWedgedSubImage3D< PIXEL > :: read( fileName.str().c_str(), volumeList );

		// compute distance matrix
		Array< PIXEL, 2 > D, W;
		stringstream fileName1;
		fileName1 << fileName.str() << ".distances.matlab";

		nbfMatlabReader mreader;
		mreader.setFileName( fileName1.str().c_str() );
		mreader.read( D );

		stringstream fileName2;
		fileName2 << fileName.str() << ".overlaps.matlab";

		mreader.setFileName( fileName2.str().c_str() );
		mreader.read( W );

		if ( sum( D.shape() ) == 0 ){
			D.resize( volumeList.size(), volumeList.size() );
			D = -1;
			W.resize( D.shape() );
			W = 1;
		}

		for ( int i = 0; i < volumeList.size(); i++ ){
			volumeList[i].setTransform( (vtkTransform*) NULL );
		}

		// other code
		
		nbfWedgedAverageImage3D< PIXEL > average;
		average.getVolumes().push_back( volumeList[0] );
		average.getVolumes().push_back( volumeList[1] );
		average.getVolumes().push_back( volumeList[2] );

		vector< nbfWedgedAverageImage3D< PIXEL > > lista;
		lista.push_back( average );

		Array< PIXEL,3 > myD( lista.size(), volumeList.size(), 19 );
		myD = -1;
		metric.getDistances( lista, volumeList, myD );
		//metric.getDistances( volumeList, D );
		
		nbfMatlabWriter w;
		w.setFileName("joder.matlab");
		w.write(myD);

		// other code

		metric.finalizeMPI();

	} else {

		metric.slaveMPI();
	
	}

	MPI::Finalize();

	return;
}