#pragma once

//using namespace blitz;

#include <em/nbfClustering.h>
#include <em/nbfHierarchicalClustering.h>
#include <em/nbfFourierImageMetricCore.h>

extern "C" {
#include <svdlib.h>
}

#include "KMeans.h"
#include "KMterm.h"
#include "KMdata.h"
#include "KMfilterCenters.h"
#include "KMlocal.h"
#include "KMrand.h"


/** Iterative Clustering Method "Loop".
*/
template< class Pixel >
class nbfLoopClustering : public nbfClustering< Pixel >
{
public:

	nbfLoopClustering();

	~nbfLoopClustering(){};

	void setFileHeader( char * h, int i ){ this->fileHeader = h; this->starting_iteration = i; };

	void setMaxIterations( int max ){ this->maxIterations = max; };

	void doClass();
	int running_mode;
};


template< class Pixel >
void nbfLoopClustering< Pixel > :: doClass()
{
	this->running_mode = NBF_LOOP_CLUSTERING_CLASS;

	stringstream currentInputFile;

	if ( this->starting_iteration - 1 < 10 ){
		currentInputFile << this->fileHeader << "_iteration_00" << this->starting_iteration - 1;
	} else if ( this->starting_iteration - 1 < 100 ){
		currentInputFile << this->fileHeader << "_iteration_0" << this->starting_iteration - 1;
	} else {
		currentInputFile << this->fileHeader << "_iteration_" << this->starting_iteration - 1;
	}
	this->loadState( currentInputFile );

	// this->computeNewReferencesMSA( this->starting_iteration );

	this->metric->get2DRepresentations( volumesFile, R, this->useRealRepresentation, this->binFactorForClassification );
	
	this->doSVDDecomposition( R, 12 );

	this->computeInitialReferencesKmeans(R,iteration);


	// Construct class averages in parallel to minimize disk access
	cout << "Generating new references..." << endl;
	for ( int i = 0; i < this->references.size(); i++ ){
		cout << "\tGenerating reference " << i << endl;
		this->metric->getImage( this->references[i] );
		if ( this->metric->getMissingWedgeCompensation() == true ){
			cout << "\tForcing class average " << i << " to have no missing-wedge" << endl;
			// SET TO HAVE NO MISSING WEDGE
			Array< Pixel, 3 > dummy3D( this->references[i].average.shape() );
			dummy3D = 1;
			this->references[i].setAccumulatedWedgeImage(dummy3D);
			//cout << "\tUpdating missing wedge for reference " << i << endl;
			//this->references[i].updateAccumulatedWedgeImage();
			//// this->metric->updateAccumulatedWedgeImage( this->references[i] );
			TinyVector< int, 2 > size;
			size = 2 * reinterpret_cast<nbfProjectionRotationMetric3D< Pixel >*>(this->metric)->B;
			//cout << "\tUpdating spherical missing wedge for reference " << i << " (size=" << size << ")" << endl;
			//this->metric->updateSphericalWedgeImage( this->references[i], size );
			Array< Pixel, 2 > dummy2D( size );
			dummy2D = 1;
			this->references[i].setSphericalWedgeImage( dummy2D );
		}
	}

	// save state
	stringstream currentOutputFile;

	if ( this->starting_iteration < 10 ){
		currentOutputFile << this->fileHeader << "_iteration_00" << this->starting_iteration;
	} else if ( this->starting_iteration < 100 ){
		currentOutputFile << this->fileHeader << "_iteration_0" << this->starting_iteration;
	} else {
		currentOutputFile << this->fileHeader << "_iteration_" << this->starting_iteration;
	}
	cout << "State being saved with prefix: " << currentOutputFile.str() << endl;
	this->saveState( currentOutputFile );
}
