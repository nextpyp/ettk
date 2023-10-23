#define NBF_PARALLEL_IMPLEMENTATION_MPI 1

#define NBF_AVERAGE_IN_RECIPROCAL_SPACE

#ifdef WIN32
	 #define NBF_VERBOSE
#endif

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

#include <em/nbfUnsupervisedLoopClustering.h>


#define PIXEL float

int main( int argc, char ** argv )
{
#if 1
	nbfMatlabReader mr;
	mr.setFileName( argv[1] );
	Array< float, 3 > A, B;
	mr.read(A);

	vtkImageData * data = vtkImageData::New();

	A.transposeSelf(thirdDim,secondDim,firstDim);
	B.resize( A.shape() );
	B = A.reverse(secondDim);
	nbfVTKInterface::blitzToVtk(B,data);

	nbfMrcWriter mrcw;
	mrcw.setFileName( argv[2] );
	mrcw.write( data );

	//// enlarge Jun's volumes
	//nbfMrcReader r;
	//r.setFileName(argv[1]);
	//vtkImageData * data = vtkImageData::New();
	//r.read(data);
	//
	//Array< double, 3 > As;
	//nbfVTKInterface::vtkToBlitz( data, As );

	//Array< double, 3 > enlarged( 64, 64, 64 );
	//enlarged = mean(As);
	//int xoffset = floor( As.rows() / 2.0 );
	//int yoffset = floor( As.cols() / 2.0 );
	//int zoffset = floor( As.depth() / 2.0 );
	//Range I( 32 - xoffset, 32 - xoffset + As.rows() - 1 );
	//Range J( 32 - yoffset, 32 - yoffset + As.cols() - 1 );
	//Range K( 32 - yoffset, 32 - zoffset + As.depth() - 1 );
	//enlarged( I, J, K ) = As;

	//// change geometry to make it compatible with mrc file
	//Array< double, 3 > B;

	//enlarged.transposeSelf(thirdDim,secondDim,firstDim);
	//B.resize( enlarged.shape() );
	//B = enlarged.reverse(secondDim);
	//nbfVTKInterface::blitzToVtk(B,data);

	//vtkImageCast * cast = vtkImageCast::New();
	//cast->SetOutputScalarTypeToFloat();
	//cast->SetInput( data );
	//cast->Update();

	//nbfMrcWriter mrcw;
	//mrcw.setFileName( argv[2] );
	//mrcw.write( cast->GetOutput() );
	//cast->Delete();

	//return 0;
#endif

	//  Initialize MPI.
	MPI_Init ( &argc, &argv );

	if ( argc != 7 ){
		cout << "Usage: subvolumesFile loopPattern hierarchicalCutOff iterations minNumElementsPerClass precomputedAlignments (if available)" << endl;
		exit(0);
	}

	int my_id, num_procs;

	//  Get this processes's rank.
	MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

	//  Find out how many processes are available.
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	cout << "Process " << my_id << " is active.\n";

	// the metric is the only thing all processes share

	nbfImageFilter< PIXEL, 3 > imfilter;
	imfilter.paddingOn(1);
	// imfilter.medianFilterOn();

	nbfFourierFilter< PIXEL, 3 > fffilter;
	//fffilter.bandPassOn(.05,.1,.01,.01); // phantoms
	//fffilter.bandPassOn(.01,.3,.05,.05); // mac239
	fffilter.bandPassOn(.05,.2,.01,1e-3); // mac239
	//fffilter.bandPassOn(0.01,.5,.001);
	//fffilter.bandPassOn(.06,.35,.001);
	//fffilter.bandPassOn(.01,.25,.001,.05); // groel
	//fffilter.bandPassOn(.05,.25,.001,.05); // groel - does not work as well
	
	//nbfCorrelationImageMetric< PIXEL, 3 > metric( &imfilter, &fffilter );
	nbfProjectionRotationMetric3D< PIXEL > metric( &imfilter, &fffilter );
	metric.setNumberOfCandidatePeaksToSearch(5);
	metric.setRotationSearchRestriction(25.0);
	metric.setTranslationSearchRestriction(15.0);
	metric.setNumberOfCandidates(3);
	
	if ( my_id == 0 ){ // Master process

		cout << "MPI - Master process:\n";
		cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
		cout << "  The number of processes is " << num_procs << "\n";

		char * spikesFileName = argv[1];
		char * patern = argv[2];

		// nuevos
		PIXEL hierarchicalThreshold = atof(argv[3]);

		int maxIter = atoi( argv[4] );

		// umbral usado como criterio para juntar clases - 
		int minElementsPerClass = atoi( argv[5] );

		vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
		nbfWedgedSubImage3D< PIXEL > :: read( spikesFileName, volumeList );

		cout << volumeList.size() << " volumes to process." << endl;

		// set volume size and cutting offset to center of subvolume
		int volumeSize = 64;
		for ( int i = 0; i < volumeList.size(); i++ ){
			TinyVector< int, 3 > tsize( volumeSize, volumeSize, volumeSize );
			volumeList[i].setDimensions( tsize );
			volumeList[i].setCutOffset( volumeSize / 2.0 );
		}

		//Array< PIXEL, 3 > lowerDimensionalRepresentation;
		//nbfFourierImageMetric< PIXEL, 3 > fmetric( &imfilter, &fffilter );
		//for ( int i = 0; i < volumeList.size(); i++ ){
		//	fmetric.setInput1( &volumeList[i] );
		//	Array< complex< PIXEL >, 1 > C;
		//	Array< PIXEL, 1 > W;
		//	fmetric.getLowDimensionRepresentation( C, W );
		//	if ( i == 0 ){
		//		lowerDimensionalRepresentation.resize( volumeList.size(), C.rows(), 4 );
		//		cout << "Initializing representation size to " << lowerDimensionalRepresentation.shape() << endl;
		//	}
		//	lowerDimensionalRepresentation( i, Range::all(), 0 ) = real(C);
		//	lowerDimensionalRepresentation( i, Range::all(), 1 ) = imag(C);
		//	lowerDimensionalRepresentation( i, Range::all(), 2 ) = sqrt( pow2( real(C) ) + pow2( imag(C) ) );
		//	lowerDimensionalRepresentation( i, Range::all(), 3 ) = W;
		//}

		//cout << "Done computing lower dimensional representation." << endl;

		//nbfMatlabWriter w;
		//w.setFileName("lowerDimensionalRepresentation.matlab");
		//w.write( lowerDimensionalRepresentation );
	
		//cout << "Now computing distance matrix." << endl;

		//// store distance matrix + overlap
		//Array< PIXEL, 3 > D( lowerDimensionalRepresentation.rows(), lowerDimensionalRepresentation.rows(), 2 );

		//for ( int i = 0; i < lowerDimensionalRepresentation.rows(); i++ ){
		//	for ( int j = i + 1; j < lowerDimensionalRepresentation.rows(); j++ ){
		//		Array< PIXEL, 1 > x1( lowerDimensionalRepresentation( i, Range::all(), 2 ) );
		//		Array< PIXEL, 1 > w1( lowerDimensionalRepresentation( i, Range::all(), 3 ) );
		//		Array< PIXEL, 1 > x2( lowerDimensionalRepresentation( j, Range::all(), 2 ) );
		//		Array< PIXEL, 1 > w2( lowerDimensionalRepresentation( j, Range::all(), 3 ) );
		//		D(i,j,0) = sqrt( sum( pow2( x1 - x2 ) * w1 * w2 ) );
		//		D(i,j,1) = sum( w1 * w2 );

		//		// symmetrize
		//		D(j,i,0) = D(i,j,0);
		//		D(j,i,1) = D(i,j,1);
		//	}
		//}

		//w.setFileName("lowerDimensionalRepresentationDistance.matlab");
		//w.write( D );
		//return 0;

		// reset original alignments
		for ( int i = 0; i < volumeList.size(); i++ ){
			volumeList[i].setTransform( (vtkTransform*)NULL );
		}

		Array< PIXEL, 3 > alignments;
		char * alignmentsFileName = argv[6];
		nbfMatlabReader mreader;
		mreader.setFileName( alignmentsFileName );
		mreader.read( alignments );

		Array< PIXEL, 4 > newAlignments( alignments.rows(), alignments.cols(), alignments.depth(), 1 );
		newAlignments( Range::all(), Range::all(), Range::all(), 0 ) = alignments;

		if ( newAlignments.size() > 0 ){
			cout << "Read alignments " << newAlignments.shape() << endl;
			// make correct distance matrix in case file was corrupted
			metric.makeDistanceMatrix(newAlignments);
		} else {
			newAlignments.resize( volumeList.size(), volumeList.size(), 19, metric.getNumberOfCandidates() );
			newAlignments = -1;
			// distance matrix with no refinement (option 3)
			metric.getDistances( volumeList, newAlignments, 3 ); // translation only

#if 0 // CODE FOR COMPUTING DISTANCE MATRIX IN 2 STAGES
			// two step distance computation:
			// 1. try only few peaks and no refinement
			// 2. refine only closest volumes

			nbfMatlabReader mreader;

			// retrieve previous state to be able to tell whether the first step was already computed
			stringstream stateFileName;
			stateFileName << "previous.state.matlab";
			mreader.setFileName( stateFileName.str().c_str() );
			Array< PIXEL, 2 > previousState;
			mreader.read( previousState );
			if ( previousState.size() == 0 ){
				previousState.resize(1,1);
				previousState = 0;
			}

			newAlignments.resize( volumeList.size(), volumeList.size(), 19, metric.getNumberOfCandidates() );
			newAlignments = -1;

			// if first pass not already done
			if ( previousState(0,0) == 0 ){

				// first try only few spherical points and no refinement (option 3)
				metric.getDistances( volumeList, newAlignments, 3 );

				mwriter.setFileName( "distances.estimate.matlab" );
				mwriter.write( newAlignments );

				// detect threshold and refine only relevant distances
				vector< PIXEL > distanceVector;
				for ( int i = 0; i < newAlignments.rows(); i++ ){
					for ( int j = i + 1; j < newAlignments.cols(); j++ ){
						distanceVector.push_back( alignments(i,j,0) );
					}
				}
				sort( distanceVector.begin(), distanceVector.end() );
				PIXEL threshold = distanceVector[floor( distanceVector.size() / 2.0 ) ];
				
				cout << "Refining only distances below " << threshold << endl;
				alignments( Range::all(), Range::all(), 0 ) = where( alignments( Range::all(), Range::all(), 0 ) <= threshold, -1, alignments( Range::all(), Range::all(), 0 ) );

				// save state as first pass already computed (=1)
				previousState = 1;
				mwriter.setFileName( stateFileName.str().c_str() );
				mwriter.write( previousState );
			}

			metric.getDistances( volumeList, alignments );

			// reset state as first pass not computed (=0)
			previousState = 0;
			mwriter.setFileName( stateFileName.str().c_str() );
			mwriter.write( previousState );
#endif // END CODE FOR COMPUTING DISTANCE MATRIX IN 2 STAGES

			nbfMatlabWriter mwriter;
			mwriter.setFileName( alignmentsFileName );
			mwriter.write( newAlignments );
		}

		if ( maxIter > 0 ){
			// Do clustering
			//nbfUnsupervisedLoopClustering< PIXEL > cluster;
			nbfLoopClustering< PIXEL > cluster;

			// Computes distanceMatrix if it's not set
			cluster.setFileHeader( patern, 0 );
			cluster.setMaxIterations( maxIter );
			cluster.setInput( volumeList );
			cluster.setMetric( &metric, newAlignments );
			cluster.setMinVolNumber( minElementsPerClass );
			cluster.setHierarchicalCutoff( hierarchicalThreshold ); // .25
			cluster.setHierarchicalMinOverlap(0.0);
			Array< PIXEL, 3 > classes;
			cluster.execute(classes);
		}

		// end MPI
		metric.finalizeMPI();
	} else {
		metric.slaveMPI();
	}
	MPI_Finalize();
	return 0;
}