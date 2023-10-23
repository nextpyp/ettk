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

	PIXEL results[4];

	//  Get this processes's rank.
	MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

	//  Find out how many processes are available.
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	cout << "Process " << my_id << " is active.\n";

	// RETRIEVE OR BUILD VOLUME LIST STRUCTURE

	vector< nbfWedgedSubImage3D< PIXEL > > volumeList;

	stringstream fileName;
	fileName << argv[1];

	nbfWedgedSubImage3D< PIXEL > :: read( fileName.str().c_str(), volumeList );

	nbfImageFilter< PIXEL, 3 > imfilter;

	nbfFourierFilter< PIXEL, 3 > fffilter;
	fffilter.bandPassOn(.05,.1,.01);
	
	//nbfCorrelationImageMetric< PIXEL, 3 > metric( &imfilter, &fffilter );
	nbfProjectionRotationMetric3D< PIXEL > metric( &imfilter, &fffilter );
	metric.setNumberOfCandidatePeaksToSearch(1);
	
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

	int jobsSent = 0;
	int jobsDone = 0;

	vector< TinyVector< int, 2 > > indexesVector;

	if ( my_id == 0 ){
		cout << "\n";
		cout << "SEARCH - Master process:\n";
		cout << "  A program using MPI, to search an array.\n";
		cout << "\n";
		cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
		cout << "\n";
		cout << "  The number of processes is " << num_procs << "\n";

		// store positions needing distance computation
		for ( int i = 0; i < volumeList.size(); i++ ){
			for ( int j = i + 1; j < volumeList.size(); j++ ){
				if ( D(i,j) == - 1 ){
					indexesVector.push_back( TinyVector< int, 2 >( i, j ) );
					cout << i << ", " << j << endl;
				}
			}
		}

		// start all available processes
		for ( int i = 1; i <  num_procs; i++ ){
			if ( jobsSent < volumeList.size() ){
				cout << "Process " << i << " processing unit " << jobsSent << endl;
				results[0] = indexesVector[jobsSent][0];
				results[1] = indexesVector[jobsSent][1];
				MPI::COMM_WORLD.Send ( &results, 4, MPI_DOUBLE, i, tag_found );
				jobsSent++;
			}
		}
	}

	if ( my_id == 0 ){

		// measure runtime
		double t1 = MPI::Wtime();

		while ( jobsDone < indexesVector.size() ){

			MPI::COMM_WORLD.Recv ( &results, 4, MPI_DOUBLE, MPI::ANY_SOURCE, MPI::ANY_TAG, status );

			jobsDone++;

			source = status.Get_source();
			tag = status.Get_tag();

			int x = floor( results[0] );
			int y = floor( results[1] );
			D( x, y ) = D( y, x ) = results[2];
			cout << "D(" << x << "," << y << ") = " << D(x,y) << ",\t";

			W( x, y )= W( y, x ) = results[3];
			cout << "W(" << x << "," << y << ") = " << W(x,y) << endl;

			nbfMatlabWriter w;
			w.setFileName( fileName1.str().c_str() );
			w.write(D);
			w.setFileName( fileName2.str().c_str() );
			w.write(W);

			if ( jobsSent < indexesVector.size() ){
				int          resultlen;
				char         hostname[MPI_MAX_PROCESSOR_NAME];
				MPI::Get_processor_name(hostname, resultlen);

				// send new job to the available node
				cout << "Process " << source << " processing unit " << jobsSent << " in " << hostname << endl;

				results[0] = indexesVector[jobsSent][0];
				results[1] = indexesVector[jobsSent][1];

				MPI::COMM_WORLD.Send ( &results, 4, MPI_DOUBLE, source, tag_found );
				jobsSent++;
			} else {
				// send termination message to all nodes
				for ( int i = 1; i < num_procs; i++ ){
					MPI::COMM_WORLD.Send ( &results, 4, MPI_DOUBLE, i, tag_done );
				}
			}
		}

		double t2 = MPI::Wtime();
		cout << "Runnign time : " << t2-t1 << " seconds" << endl;

		// fill diagonal
		for ( int i = 0; i < volumeList.size(); i++ ){
			D(i,i) = 0;
		}

		cout << "D = " << D << endl;

		cout << "W = " << W << endl;

		nbfMatlabWriter w;
		w.setFileName( fileName1.str().c_str() );
		w.write(D);
		w.setFileName( fileName2.str().c_str() );
		w.write(W);

	} else {
		source = 0;
		tag = tag_target;

		while ( tag != tag_done ){

			MPI::COMM_WORLD.Recv ( &results, 4, MPI::DOUBLE, source, MPI::ANY_TAG, status );

			tag = status.Get_tag();

			if ( tag != tag_done ){
				// do processing
				metric.setInput1( &( volumeList[ floor( results[0] ) ] ) );
				metric.setInput2( &( volumeList[ floor( results[1] ) ] ) );
				
				results[2] = metric.getDistance();
				results[3] = metric.getWedgeOverlap();

				// notify master node we have finished processing
				MPI::COMM_WORLD.Send ( &results, 4, MPI_DOUBLE, 0, tag_found );
			}
		}
	}

	MPI::Finalize();

	return;
}