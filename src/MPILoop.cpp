#define NBF_PARALLEL_IMPLEMENTATION_MPI 1

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

#define BUILD_VOLUME_LIST 0

int main( int argc, char ** argv )
{
	nbfMatlabReader r;
	stringstream rs;
	rs << argv[1] << ".blitz";
	r.setFileName( rs.str().c_str() );
	Array< float, 3 > A;
	r.read(A);
	nbfMrcWriter w;
	stringstream s;
	s << argv[1] << "_crop.mrc";
	w.setFileName( s.str().c_str() );
	w.write(A);
	return 0;

	//  Initialize MPI.
	MPI_Init ( &argc, &argv );

	if ( argc > 10 || argc == 1 ){
		cout << "Usage: spikesFileName alignmentFileName (optional)" << endl;
		exit(0);
	}

	int my_id, num_procs;

	//  Get this processes's rank.
	MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

	//  Find out how many processes are available.
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	cout << "Process " << my_id << " is active.\n";

	// the metric is the only think all processes share

	nbfImageFilter< PIXEL, 3 > imfilter;

	nbfFourierFilter< PIXEL, 3 > fffilter;
	fffilter.bandPassOn(.05,.1,.01);
	
	//nbfCorrelationImageMetric< PIXEL, 3 > metric( &imfilter, &fffilter );
	nbfProjectionRotationMetric3D< PIXEL > metric( &imfilter, &fffilter );
	metric.setNumberOfCandidatePeaksToSearch(10);
	
	if ( my_id == 0 ){ // Master process

		cout << "MPI - Master process:\n";
		cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
		cout << "  The number of processes is " << num_procs << "\n";

		char * spikesFileName = argv[1];
		char * patern = argv[2];
		// nuevos
		int maxIter = atoi( argv[3] );

		// parametro para calcular la "pertenencia" de los volumenes a las clases - B = 0  usa un valor calculado a partir de la matriz de distancias
		// A = exp(-d^2/B)
		double B = atof( argv[4] );
		// factor se usa para cortar el histograma de distancias para recalcular un promedio 
		// umbral_corte = factor * media (media es el promedio de las distancias de los volumenes a la media)
		double factor = atof( argv[5] );
		// umbral usado como criterio para asignar volumenes a clases - 
		double similarityTh = atof( argv[6] );
		// umbral usado como criterio para juntar clases - 
		double meanSimilarityTh = atof( argv[7] );

		cout << "Reading volumes...";

		vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
		nbfWedgedSubImage3D< PIXEL > :: read( spikesFileName, volumeList );

		//for ( int i = 0; i < volumeList.size(); i++ ){
		//	volumeList[i].setCutOffset( 8.0 );
		//}

		//for ( int i = 0; i < volumeList.size(); i++ ){
		//	volumeList[i].imageToFullVolumeOn();
		//}

		// reset original alignments
		for ( int i = 0; i < volumeList.size(); i++ ){
			volumeList[i].setTransform( (vtkTransform*)NULL );
		}

		cout << " done" << endl;

		Array< PIXEL, 3 > alignments;
		char * alignmentsFileName;
		if ( argc == 9 ){
			alignmentsFileName = argv[8];
			nbfMatlabReader mreader;
			mreader.setFileName( alignmentsFileName );
			mreader.read( alignments );
			cout << "Read alignments " << alignments.shape() << endl;
		} else {
			alignments.resize( volumeList.size(), volumeList.size(), 19 );
			alignments = -1;
			metric.getDistances( volumeList, alignments );
			nbfMatlabWriter mwriter;
			mwriter.setFileName( alignmentsFileName );
			mwriter.write( alignments );
		}

		// Do clustering
		nbfLoopClustering< PIXEL > cluster;

		// Computes distanceMatrix if it's not set
		cluster.setFileHeader( patern );
		cluster.setParameterB( B );
		cluster.setFactor( factor );
		cluster.setSimilarityThreshold( similarityTh );
		cluster.setMeanSimilarityThreshold( meanSimilarityTh );
		cluster.setMaxIterations( maxIter );
		cluster.setInput( volumeList );
		cluster.setMetric( &metric, alignments );
		cluster.setMinVolNumber(3);
		
		cluster.execute();

		// end MPI

		metric.finalizeMPI();

	} else {

		metric.slaveMPI();
	
	}

	MPI_Finalize();

	return 0;
}