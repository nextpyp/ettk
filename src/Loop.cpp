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

#include <em/nbfLoopClustering.h>

#include "vtkDebugLeaks.h"

#define PIXEL double

#define BUILD_VOLUME_LIST 0

void main( int argc, char ** argv )
{
	
	if ( argc > 8 || argc ==1 ){
		cout << "Usage: spikesFileName distanceMatrixFileName( optional )" << endl;
		exit(0);
	}

	char * spikesFileName = argv[1];
	char * patern = argv[2];
	// nuevos
	int maxIter = atoi( argv[3] );
	
	// parametro para calcular la "pertenencia" de los volumenes a las clases - B = 0  usa un valor calculado a partir de la matriz de distancias
	// A = exp(-d^2/B)
	double B = atof( argv[4] );
	// factor se usa para cortar el histograma de distancias para recalcular un promedio 
	// umbral_corte = factor * media - media es el promedio de las distancias de los volumenes a la media 
	double factor = atof( argv[5] );
	// umbral usado como criterio para juntar clases - 
	double similarityTh = atof( argv[6] );

	Array< PIXEL, 2 > distanceMatrix;
	char * distanceMatrixFileName;
	if ( argc == 8 ){
		distanceMatrixFileName = argv[7];
		nbfMatlabReader mreader;
		mreader.setFileName( distanceMatrixFileName );
		mreader.read( distanceMatrix );
	}
	
	cout << "Reading volumes...";
	
	vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
	nbfWedgedSubImage3D< PIXEL > :: read( spikesFileName , volumeList );

	//for ( int i = 0; i < volumeList.size(); i++ ){
	//	volumeList[i].setCutOffset( 8.0 );
	//}

	//for ( int i = 0; i < volumeList.size(); i++ ){
	//	volumeList[i].imageToFullVolumeOn();
	//}

	cout << " done" << endl;

	// Set metric
	nbfImageFilter< PIXEL, 3 > imfilter;
	//imfilter.windowOff();
	//imfilter.medianFilterOn();

	nbfFourierFilter< PIXEL, 3 > fffilter;
	fffilter.bandPassOn(.05,.1,.01);
	//fffilter.lowPassOn( TinyVector<PIXEL,3>(.2,.2,.2) );


	nbfCorrelationImageMetric< PIXEL, 3 > metric( &imfilter, &fffilter );
	//nbfProjectionRotationMetric3D< PIXEL > metric( &imfilter, &fffilter );
	//metric.setAngularSearchRestriction(10);
	
	cout << "D = [" << min(distanceMatrix) << "," << max(distanceMatrix) << "]" << endl;

	// Do clustering
	nbfLoopClustering< PIXEL > cluster;

	// Computes distanceMatrix if it's not set
	cluster.setFileHeader( patern );
	cluster.setParameterB( B );
	cluster.setFactor( factor );
	cluster.setSimilarityThreshold( similarityTh );
	cluster.setMetric( & metric , distanceMatrix);
	cluster.setMaxIterations( maxIter );

	vector< nbfWedgedSubImage3D< PIXEL > > volumeList2;
	
	cluster.setInput( volumeList );
	cluster.execute();



}