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
	int seed = atoi(argv[2]);
	//int volumeSize = atoi(argv[2]);
	//double wedge = atof(argv[3]);
	double cccth = atof(argv[3]);
	
	//nbfWedge3D< PIXEL > wedge3D;
	//wedge3D.set(-wedge,wedge);
	//nbfWedgedSubImage3D< PIXEL > reference;

	//nbfWedgedSubImage3D< PIXEL > cutVolume;
	//cutVolume.getWedge()->set( -wedge, wedge );
	//cutVolume.setDimensions( TinyVector<int,3>( volumeSize, volumeSize, volumeSize ) );
	//cutVolume.setCutOffset( 8.0 );
	//cutVolume.setCutOffset( 50.0 );
	//cutVolume.setMagnification( TinyVector< PIXEL, 3 >(2.0,2.0,2.0) );

	// RETRIEVE OR BUILD VOLUME LIST STRUCTURE

	vector< nbfWedgedSubImage3D< PIXEL > > volumeList;

#if BUILD_VOLUME_LIST
	for ( int tomograms = 0; tomograms < argc - 5; tomograms++ ){
	
		stringstream tomogramFile;
		tomogramFile << argv[ tomograms + 5 ] << ".rec";

		stringstream spikeLocationsFile;
		spikeLocationsFile << argv[ tomograms + 5 ] << "pairs.matlab";

		nbfMatlabReader mread;
		mread.setFileName( spikeLocationsFile.str().c_str() );

		Array< float, 2 > points;
		cout << "Reading " << spikeLocationsFile.str() << "...";
		mread.read( points );
		// FIX READING ERRORS
		points = where( points - floor(points) > .5, ceil(points), floor(points) );
		if ( points.size() > 0 ){
			cout << "done.\n";
		}
		else{
			cout << "error reading " << spikeLocationsFile.str() << endl;
		}

		nbfMrcReader reader;
		reader.setFileName( tomogramFile.str().c_str() );
		cout << "Reading " << tomogramFile.str() << "...";

		int dimY = reader.getDimY();
		cout << "done.\n";

		for ( int i = 1; i < points.rows(); i+=2 ){
		//for ( int i = 1; i < 2; i++ ){

			nbfWedgedSubImage3D< PIXEL > currentVolume;
			currentVolume = cutVolume;

			currentVolume.setFileName( tomogramFile.str().c_str() );
			currentVolume.setPosition( TinyVector< PIXEL, 3 >( points(i-1,0), dimY - 1 - points(i-1,1), points(i-1,2) ) );
			//currentVolume.setNormal( TinyVector< PIXEL, 3 >( 0,0,1 ) );
			//currentVolume.setNormal( TinyVector< PIXEL, 3 >( points(i,0) - points(i-1,0), - points(i,1) + points(i-1,1), 1 ) );
			currentVolume.setNormal( TinyVector< PIXEL, 3 >( points(i,0) - points(i-1,0), - points(i,1) + points(i-1,1), points(i,2) - points(i-1,2) ) );
			volumeList.push_back( currentVolume );
		}
	}

	nbfWedgedSubImage3D< PIXEL > :: write( argv[1], volumeList );
#else
	nbfWedgedSubImage3D< PIXEL > :: read( argv[1], volumeList );
#endif

	//stringstream tomogramFile;
	//tomogramFile << argv[5];

	nbfImageFilter< PIXEL, 3 > imfilter;
	//imfilter.windowOff();
	//imfilter.medianFilterOn();

	nbfFourierFilter< PIXEL, 3 > fffilter;
	//fffilter.lowPassOn( TinyVector<PIXEL,3>(.15,.15,.15) );
	fffilter.bandPassOn(.05,.1,.01);
	
	//nbfCorrelationImageMetric< PIXEL, 3 > metric( &imfilter, &fffilter );
	nbfProjectionRotationMetric3D< PIXEL > metric( &imfilter, &fffilter );
	metric.setNumberOfCandidatePeaksToSearch(5);
	
	// compute distance matrix
	Array< PIXEL, 2 > D, W;

	for ( int i = 0; i < volumeList.size(); i++ ){
		volumeList[i].setTransform( (vtkTransform*) NULL );
	}
	// metric.getMatrix( argv[1], volumeList, D, &W );
	
	//return;

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileTypeToBinary();

	PriorityQueue< TinyVector< PIXEL, 4 >, greaterPD< PIXEL, 3 > , vector< TinyVector< PIXEL, 4 > > > pqNB;

	//Array< short, 3 > A;
	//nbfVTKInterface::vtkToBlitzReference(referenceImage,A);

	int numberOfPasses = 1;

	nbfWedgedAverageImage3D< PIXEL > wedgedAverage;

	if ( ( seed < 0 ) || ( seed > volumeList.size() - 1 ) ){
		cerr << "Invalid seed index: " << seed << " , Range(0," << volumeList.size() - 1 << ")" << endl;
		exit(0);
	} else{
		wedgedAverage.getVolumes().push_back( volumeList[seed] );
	}

	nbfTimer timer;
	timer.start();

	for ( int j = 1; j <= numberOfPasses; j++ ){

		pqNB.clear();

		metric.setInput1(( &wedgedAverage ) );

		for ( int i = 0; i < volumeList.size(); i++ ){
		//for ( int i = 0; i < 2; i++ ){

			volumeList[i].setTransform( (vtkTransform*)NULL );
			metric.setInput2( &(volumeList[i]) );

			cout << i << ":";

			timer.start();
			PIXEL d = metric.getDistance();
			//PIXEL d = 0;
			timer.stop();

			cout << ", d=" << d <<  ", t = " << timer.elapsedSeconds() << endl;

			volumeList[i].setTransform( metric.getTransform() );

			pqNB.push( TinyVector< PIXEL, 4 >( i, i, i, d ) );
		}

		int count = 0;
		int maxNumberForAveraging = 500;

		// add to average if close enough to reference

		TinyVector< PIXEL, 4 > nextVolume = pqNB.top();

		// reset average image
		wedgedAverage.getVolumes().clear();

		// return;

		do {
			count++;

			cout << "Adding volume " << nextVolume[0] << " to average, d = " << nextVolume[3] << endl;

			wedgedAverage.getVolumes().push_back( volumeList[ nextVolume[0] ] );

			// save aligned image
			stringstream fileName1;
			if ( nextVolume[0] < 10 ){
				fileName1 << "volume.cut.00" << nextVolume[0] << ".vtk";
			}
			else {
				if ( nextVolume[0] < 100 ){
				fileName1 << "volume.cut.0" << nextVolume[0] << ".vtk";
				}
				else{
					fileName1 << "volume.cut." << nextVolume[0] << ".vtk";
				}
			}

			vtkImageData * aligned = vtkImageData::New();
			volumeList[ nextVolume[0] ].getImage( aligned );
			writer->SetFileName( fileName1.str().c_str() );
			writer->SetInput( aligned );
			// writer->Write();

			aligned->Delete();

			pqNB.pop();
			nextVolume = pqNB.top();

		} while ( ( pqNB.size() > 0 ) && ( nextVolume[3] < cccth ) && ( count < maxNumberForAveraging ) );

		//return;

		stringstream fileName3;
		fileName3 << argv[5] << ".average.real." << j << ".vtk";
		vtkImageData * averageVtk = vtkImageData::New();

		//wedgedAverage.setCutSize( 70 );
		//wedgedAverage.setCutOffset( 20 );
		wedgedAverage.getImage( averageVtk );
		//wedgedAverage.setCutSize( cutVolume.getDimensions()[0] );
		//wedgedAverage.setCutOffset( cutVolume.getCutOffset() );

		//nbfVTKInterface::blitzToVtk( average, averageVtk );
		writer->SetInput( averageVtk );
		writer->SetFileName(fileName3.str().c_str());
		writer->Write();
		averageVtk->Delete();

		return;
	}

	timer.stop();
	cout << timer.elapsedSeconds() << " seconds elapsed." << endl;

	nbfWedgedSubImage3D< PIXEL > :: write( "allspikesaligned.txt", volumeList );
	//referenceImage->Delete();
	//writer->Delete();
}