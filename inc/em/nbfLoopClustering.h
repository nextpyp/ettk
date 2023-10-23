#pragma once

//using namespace blitz;

#include <em/nbfClustering.h>
#include <em/nbfHierarchicalClustering.h>
#include <em/nbfFourierImageMetricCore.h>

#include <nbfPolarDomain.h>

extern "C" {
#include <svdlib.h>
}

#include "KMeans.h"
#include "KMterm.h"
#include "KMdata.h"
#include "KMfilterCenters.h"
#include "KMlocal.h"
#include "KMrand.h"

#define NBF_LOOP_CLUSTERING_AUTO	-1
#define NBF_LOOP_CLUSTERING_CENTER	0
#define NBF_LOOP_CLUSTERING_CLASS	1
#define NBF_LOOP_CLUSTERING_REFINE	2
#define NBF_LOOP_CLUSTERING_MRA		3
#define NBF_LOOP_CLUSTERING_CTF		4

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
	
	void setMinDistanceBetweenReferences( double Th ){ this->minDistanceBetweenReferences = Th;};

	void saveState( stringstream &, bool = false );
	void loadState( stringstream & );

	void setMinVolNumber( int a ){ this->minVolumeNumber = a; }

	void setHierarchicalClasses( int c, int d ){ this->hierarchicalClasses = c; this->hierarchicalClassesHigh = d; };
	void setHierarchicalToBuildEntireTree( Pixel t ){ this->hierarchicalToBuildEntireTree = t; }

	void setHierarchicalMinOverlap( double o ){ this->hierarchicalMinOverlap = o;};

	void setInitialClasses( char * );

	void setSymmetryFactor( int i = 0 ){ this->symmetryFactor = i; }

	// set top percentage of volumes to include in average (\in [0,1])
	void setDistanceTopCutoff( Pixel g, Pixel p, Pixel q ){
		this->distanceTopCutoffPre = blitz :: extrema :: max( 0.0, blitz :: extrema :: min( g, 1.0 ) );
		this->distanceTopCutoff = blitz :: extrema :: max( 0.0, blitz :: extrema :: min( p, 1.0 ) ); 
		this->distanceTopCutoffSelection = blitz :: extrema :: max( 0.0, blitz :: extrema :: min( q, 1.0 ) ); 
	} 

	void setRefinementIterations( int b ){ 
		this->refinementIterations = b;
		if ( b < 0 ){
			this->refinementOnly = true; 
		}
	}

	void useRealRepresentationOn( int d = 20 ){ this->useRealRepresentation = d; }
	void useRealRepresentationOff(){ this->useRealRepresentation = 0; }

	void symmetrize( nbfWedgedAverageImage3D< Pixel > &, Pixel = 0, bool = true );

	void setRunningMode( int i ){ this->running_mode = i; };

	void doClass();
	void doRefine( bool = false, bool = false, int = 2, int volumeCutOffset = 0 );
	void doMra(bool,Pixel=0,Pixel=0,Pixel=0,int=0);
	void doClassAndSelect();

	void doCTFClass();

	void do2DClass();
	void do2DRefine(){};
	void do2DMra(){};


	Pixel binFactorForClassification;

protected:

	void doClustering( Array< Pixel, 3 > & );

	void denoise( vector< nbfWedgedSubImage3D< Pixel > > &, Array< Pixel, 3 > &, vector< nbfWedgedAverageImage3D< Pixel > > &, Pixel, Pixel);

	void computeInitialReferences( Array< Pixel, 2 > &, int );
	void computeInitialReferencesKmeans( Array< Pixel, 3 > &, int );
	void computeInitialArbitraryReferences();

	void asignBestAlignmentAndClean();

	void computeNewReferences(int);
	void computeNewReferencesHierarchical(int);
	void computeNewReferencesMSA(int);
	void computeNewReferencesSpectral(int);
	void compareAndMergeReferences();
	
	void alignReferencesToStrongest( bool = false );
	void alignReferencesToStrongestBundle( bool = false );
	void alignReferencesToStrongestBundleExternal( bool = false );
	void alignReferencesToCommonFrame();
	void alignReferencesInBundles( bool = false, int = 2 );

	void buildInitialReference( nbfWedgedAverageImage3D< Pixel > &, int );
	void buildInitialReferences();
	
	void bundleAlignment( vector< nbfWedgedAverageImage3D< Pixel > > &, int = 2 );
	void bundleAlignment( vector< nbfWedgedSubImage3D< Pixel > > & );

	void alignVolumesToReferences( bool = false, int = 0 );
	void alignVolumesToExternalReferences( int = 0 );

	void doPCADecomposition( Array< Pixel, 3 > &, int = 10 );
	void doSVDDecomposition( Array< Pixel, 3 > &, int = 10 );

	void doSpectralSVDDecomposition( Array< Pixel, 3 > &, int = 10 );

	void imputeInReciprocalSpace( Array< Pixel, 3 > &, Array< Pixel, 2 > & );

	Array< Pixel, 4 > alignmentToReferences;
	Array< Pixel, 4 > alignmentBetweenReferences;

	Array< Pixel, 1 > referenceRadii;

	vector< nbfWedgedAverageImage3D< Pixel > > references;
	vector< nbfWedgedSubImage3D< Pixel > > precomputedReferences;
	vector< int > referencesIndexes;

	Array< Pixel , 2 > weights;

	//Array< Pixel, 3 > oldClasses;

	vector< vector< int > > classesSelected;

	int minVolumeNumber;

	int maxIterations;

	char * fileHeader;

	double minDistanceBetweenReferences;

	bool convergence;

	int hierarchicalClasses;

	int hierarchicalClassesHigh;

	Pixel hierarchicalToBuildEntireTree;

	double hierarchicalMinOverlap;

	bool usingExternalReferences;

	vector< nbfWedgedSubImage3D< Pixel > > externalReferences;

	int starting_iteration;

	int symmetryFactor;

	Pixel distanceTopCutoffPre, distanceTopCutoff, distanceTopCutoffSelection;

	int refinementIterations;

	bool refinementOnly;

	int useRealRepresentation;

	int running_mode;
};

template< class Pixel >
nbfLoopClustering< Pixel > :: nbfLoopClustering()
: nbfClustering< Pixel >()
{
	this->convergence = false;
	this->maxIterations = 5;

	// By default: don't merge means
	this->minDistanceBetweenReferences = 0;
	this->fileHeader = "loop";
	this->minVolumeNumber = 2;
	this->hierarchicalMinOverlap = 0;

	this->usingExternalReferences = false;
	this->starting_iteration = 0;

	this->hierarchicalClasses = 0;
	this->hierarchicalClassesHigh = 0;
	this->hierarchicalToBuildEntireTree = 0;
	this->symmetryFactor = 0;

	this->distanceTopCutoffPre = 1;
	this->distanceTopCutoff = 1;
	this->distanceTopCutoffSelection = 1;

	this->refinementIterations = 0;

	this->binFactorForClassification = 1;

	this->refinementOnly = false;

	this->useRealRepresentation = 0;

	this->running_mode = NBF_LOOP_CLUSTERING_AUTO;
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: setInitialClasses( char * initial_classes )
{
	if ( initial_classes != NULL ){
		this->usingExternalReferences = true;

		// vector< nbfWedgedSubImage3D< Pixel > > externalReferences;
		nbfWedgedSubImage3D< Pixel > :: read( initial_classes, this->externalReferences );

		this->references.clear();
		for ( int i = 0; i < externalReferences.size(); i++ ){
			nbfWedgedAverageImage3D< Pixel > av( externalReferences );
			Array< Pixel, 3 > A( externalReferences.size(), 17, 1 );
			A = 0;
			A( i, 0, Range::all() ) = 1;
			// set identity as the transformation
			double matrix[16];
			vtkMatrix4x4::Identity(matrix);
			for ( int i = 0; i < A.rows(); i++ ){
				for ( int j = 0; j < 16; j++ ){
					A( i, 1 + j, 0 ) = matrix[j];
				}
				//A( i, Range(1,toEnd), 0 ) = this->alignments( 0, 0, Range(3,toEnd), 0 );
			}
			av.setAlignments( A );
			this->references.push_back( av );
		}

		this->classes.resize( this->references.size(), this->input.size(), 1 );
	}
}

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

	//// set input alignments to closest reference from alignmentsToReferences
	//if ( this->metric->getMissingWedgeCompensation() == true ){
	//	this->computeNewReferencesSpectral( this->starting_iteration );
	//	//this->computeNewReferencesMSA( this->starting_iteration );
	////	this->computeNewReferencesHierarchical( this->starting_iteration );
	//} else {
	//	this->computeNewReferencesMSA( this->starting_iteration );
	//}
	this->computeNewReferencesMSA( this->starting_iteration );

	// Construct references in parallel to minimize disk access
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
			// cout << "\tUpdating missing wedge for reference " << i << endl;
			// this->references[i].updateAccumulatedWedgeImage();
			// // this->metric->updateAccumulatedWedgeImage( this->references[i] );
			TinyVector< int, 2 > size;
			size = 2 * reinterpret_cast<nbfProjectionRotationMetric3D< Pixel >*>(this->metric)->B;
			// cout << "\tUpdating spherical missing wedge for reference " << i << " (size=" << size << ")" << endl;
			// this->metric->updateSphericalWedgeImage( this->references[i], size );
			Array< Pixel, 2 > dummy2D( size );
			dummy2D = 1;
			this->references[i].setSphericalWedgeImage( dummy2D );
		} else {
			cout << "\tForcing class average " << i << " to have no missing-wedge" << endl;
			// SET TO HAVE NO MISSING WEDGE
			Array< Pixel, 3 > dummy3D( this->references[i].average.shape() );
			dummy3D = 1;
			this->references[i].setAccumulatedWedgeImage(dummy3D);
			TinyVector< int, 2 > size;
			size = 2 * reinterpret_cast<nbfProjectionRotationMetric3D< Pixel >*>(this->metric)->B;
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

template< class Pixel >
void nbfLoopClustering< Pixel > :: doClassAndSelect()
{
	// do first classification
	this->doClass();

	////////////////////////
	// CLASSES OF CLASSES //
	////////////////////////

	//// store class sizes so we can then figure out which classes to use as references
	//Array< Pixel, 1 > classSizes( this->classes.rows() );
	//Array< Pixel, 2 > vClasses( this->classes( Range::all(), Range::all(), 0 ) );
	//secondIndex j;
	//classSizes = sum( vClasses, j );

	//// figure out the second order class memberships and save it in *_averages.txt

	//if ( this->hierarchicalClasses != this->hierarchicalClassesHigh ){

	//	// re-assign classification parameters
	//	this->hierarchicalClasses = this->hierarchicalClassesHigh; 
	//	this->distanceTopCutoff = 1; // ( this->hierarchicalClasses - 1 ) / this->hierarchicalClasses; // do not exclude outliers
	//	this->distanceTopCutoffPre = 1; // ( this->hierarchicalClasses - 1 ) / this->hierarchicalClasses; // do not exclude outliers

	//	// Compute reduced representations
	//	Array< Pixel, 3 > R;
	//	stringstream volumesFile;
	//	int iteration = this->starting_iteration;
	//	volumesFile << this->fileHeader << "_iteration_00" << iteration << "_averages.txt";
	//	this->metric->getRepresentations( volumesFile, R, this->useRealRepresentation, this->binFactorForClassification, this->symmetryFactor );
	//	// Reduce dimensionality by SVD decomposition
	//	this->doSVDDecomposition(R,16);
	//	this->computeInitialReferencesKmeans(R,iteration);

	//	nbfWedgedSubImage3D< Pixel > :: read( volumesFile.str().c_str(), this->precomputedReferences );

	//	// determine reference with most volumes for each cluster of references
	//	Array< int, 1 > biggest( this->classes.rows() );
	//	Array< int, 1 > biggestIndex( this->classes.rows() );
	//	biggest = 0;
	//	for ( int i = 0; i < this->classes.rows(); i++ ){
	//		for ( int j = 0; j < this->classes.cols(); j++ ){
	//			if ( this->classes( i, j, 0 ) > 0 ){
	//				if ( classSizes(j) > biggest(i) ){
	//					biggest(i) = classSizes(j);
	//					biggestIndex(i) = j;
	//				}
	//			}
	//		}
	//	}
	//	// determine class to use as global reference (-1)
	//	TinyVector< int, 1 > globalMax = maxIndex( classSizes );

	//	// encode memberships into cut offset field
	//	int cummulativeOffset = 2;
	//	for ( int i = 0; i <  this->classes.rows(); i++ ){
	//		Pixel offset = cummulativeOffset;
	//		if ( biggestIndex(i) == globalMax[0] ){
	//			offset = 1;
	//		} else {
	//			cummulativeOffset++;
	//		}
	//		for ( int j = 0; j <  this->classes.cols(); j++ ){
	//			if ( this->classes( i, j, 0 ) > 0 ){
	//				if ( biggestIndex(i) == j ){
	//					this->precomputedReferences[j].setCutOffset( - offset );
	//				} else {
	//					this->precomputedReferences[j].setCutOffset( offset );
	//				}
	//			}
	//		}
	//	}

	//	nbfWedgedSubImage3D< Pixel > :: write( volumesFile.str().c_str(), this->precomputedReferences );
	//}

	///////////////////////////////
	// AUTOMATIC CLASS SELECTION //
	///////////////////////////////

	// Sort classes by increasing distance to reference of previous iteration
	
	// Attempt to read selected global average from previous iteration
	stringstream latestReferenceFile;
	if ( this->starting_iteration - 1 < 10 ){
		latestReferenceFile << this->fileHeader << "_iteration_00" << this->starting_iteration - 1 << "_refined_selected_global_average.mrc";
	} else if ( this->starting_iteration - 1 < 100 ){
		latestReferenceFile << this->fileHeader << "_iteration_0" << this->starting_iteration - 1 << "_refined_selected_global_average.mrc";
	} else {
		latestReferenceFile << this->fileHeader << "_iteration_" << this->starting_iteration - 1 << "_refined_selected_global_average.mrc";
	}
	nbfMrcReader reader;
	vtkImageData * data = vtkImageData :: New();
	reader.setFileName( latestReferenceFile.str().c_str() );
	reader.read( data );

	// If not succesful, try reading reference 0 from previous iteration
	stringstream latestReferenceFileAlternate;
	bool usingAlternate = false;
	if ( data->GetDimensions()[0] == 0 ){
		if ( this->starting_iteration - 1 < 10 ){
			latestReferenceFileAlternate << this->fileHeader << "_iteration_00" << this->starting_iteration - 1 << "_refined_selected_average_0.mrc";
		} else if ( this->starting_iteration - 1 < 100 ){
			latestReferenceFileAlternate << this->fileHeader << "_iteration_0" << this->starting_iteration - 1 << "_refined_selected_average_0.mrc";
		} else {
			latestReferenceFileAlternate << this->fileHeader << "_iteration_" << this->starting_iteration - 1 << "_refined_selected_average_0.mrc";
		}
		reader.setFileName( latestReferenceFileAlternate.str().c_str() );
		reader.read( data );
		usingAlternate = true;
	}

	// UPDATE THE METRIC BEFORE DOING THE ALIGNMENTS
	// We are currently using the 'classification' metric (mode 1) but we now need to
	// switch to the 'refine' metric (mode 2) to compute the distances properly.
	// Since we do not have all the 'refine' parameters, all we do in increase the
	// mask apodization which is usually set to 0 in the classification mode.
	if ( this->metric->imageFilter->getSigma() == 0 ){
		Pixel x, y, z, sigma; bool cylinder;
		this->metric->imageFilter->getMaskSize( x, y, z, sigma, cylinder );
		sigma = 4;
		this->metric->imageFilter->setMaskSize( x, y, z, sigma, cylinder );
	}
	// Set translation tolerance to 0 (no search for shift)
	this->metric->setTranslationSearchRestriction( 0 );

	// retrieve current classes from file
	vector< nbfWedgedSubImage3D< Pixel > > latestClasses;
	stringstream volumesFile;
	if ( this->starting_iteration < 10 ){
		volumesFile << this->fileHeader << "_iteration_00" << this->starting_iteration << "_averages.txt";
	} else if ( this->starting_iteration < 100 ){
		volumesFile << this->fileHeader << "_iteration_0" << this->starting_iteration << "_averages.txt";
	} else {
		volumesFile << this->fileHeader << "_iteration_" << this->starting_iteration << "_averages.txt";
	}
	cout << "Reading references from file " << volumesFile.str().c_str() << endl;
	nbfWedgedSubImage3D< Pixel > :: read( volumesFile.str().c_str(), latestClasses );	

	// If reference from previous iteration available
	if ( data->GetDimensions()[0] > 0 ){

		vector< nbfWedgedImage3D< Pixel > * > list2;
		for ( int i = 0; i < latestClasses.size(); i++ ){
			list2.push_back( &( latestClasses[i] ) );
		}

		nbfWedgedSubImage3D< Pixel > latestReference;
		latestReference = latestClasses[0];
		if ( usingAlternate == false ){
			latestReference.setFileName( latestReferenceFile.str().c_str() );
		} else {
			latestReference.setFileName( latestReferenceFileAlternate.str().c_str() );
		}

		vector< nbfWedgedImage3D< Pixel > * > list1;
		list1.push_back( & latestReference );

		vector< TinyVector< int, 2 > > pos;
		for ( int i = 0; i < latestClasses.size(); i++ ){
			TinyVector< int, 2 > p( 0, i );
			pos.push_back(p);
		}

		Array< Pixel, 3 > refAlignments( list2.size(), 17, 1 );
		refAlignments = -1;
		this->metric->getDistances( list1, list2, pos, refAlignments, 2 );

		cout << "Distances to " << latestReference.getFileName() << ": " << refAlignments( Range :: all(), 0, 0 ) << endl;

		// Sort distances and apply selection threshold
		vector< Pixel > sortedDistances;
		for ( int i = 0; i < latestClasses.size(); i++ ){
			sortedDistances.push_back( refAlignments( i, 0, 0 ) );
		}
		sort( sortedDistances.begin(), sortedDistances.end() );
		Pixel th = sortedDistances[ this->distanceTopCutoffSelection * ( sortedDistances.size() - 1 ) ];

		bool firstDone = false;
		for ( int i = 0; i < latestClasses.size(); i++ ){
			if ( refAlignments( i, 0, 0 ) == sortedDistances[0] ){
				latestClasses[i].setCutOffset( -1 );
				cout << "Class " << i << " selected as next reference." << endl;
			} else {
				if ( refAlignments( i, 0, 0 ) <= th ){
					latestClasses[i].setCutOffset( 1 );
				} else {
					latestClasses[i].setCutOffset( 0 );
				}
			}
		}

		// // keep all classes
		// for ( int i = 0; i < latestClasses.size(); i++ ){
			// latestClasses[i].setCutOffset( -(i+1) );
		// }
		
		// Save class selection based on correlation to latest reference
		nbfWedgedSubImage3D< Pixel > :: write( volumesFile.str().c_str(), latestClasses );

	} else {

		// When the latest reference cannot be retrieved from a file we assume that this is the initial classification.
		// The criteria for class selection in this case is based in how trimeric each class is.
		// 3-folded-ness is estimated as the inverse distance to a 60-degree rotated copy of each class.

		vector< Pixel > symmetric( latestClasses.size() );
		for ( int i = 0; i < latestClasses.size(); i++ ){
			symmetric[i] = 0;
			for ( int f = 1; f < 2; f++ ){
		
				// initialize metric to measure distances
				nbfFourierImageMetric< Pixel, 3 > fMetric( this->metric->imageFilter, this->metric->fourierFilter );
				fMetric.setTranslationSearchRestriction( 0 );
				fMetric.setToComputeOverlapNormalizedDistances( this->metric->getToComputeOverlapNormalizedDistances() );
				fMetric.setToUseMutualCorrelation( this->metric->getToUseMutualCorrelation() );
				fMetric.setInput1( &latestClasses[i] );
				fMetric.setInput2( &latestClasses[i] );

				// evaluate distance at 60 degree-rotation
				vtkTransform * t = vtkTransform :: New();
				t->RotateZ( f * 360.0 / 6.0 );
				fMetric.executeFourierNewHalf( t );

				// assign negative distance as measure of 3-foldedness
				symmetric[i] = - fMetric.getDistance();
				t->Delete();
			}
			cout << "3-foldedness(" << i << ") = " << - symmetric[i] << endl;
		}

		vector< Pixel > sortedSymmetry;
		for ( int i = 0; i < latestClasses.size(); i++ ){
			sortedSymmetry.push_back( symmetric[i] );
		}
		sort( sortedSymmetry.begin(), sortedSymmetry.end() );
		Pixel th = sortedSymmetry[ this->distanceTopCutoffSelection * ( sortedSymmetry.size() - 1 ) ];

		bool firstDone = false;
		for ( int i = 0; i < latestClasses.size(); i++ ){
			if ( symmetric[i] == sortedSymmetry[0] ){
				latestClasses[i].setCutOffset( -1 );
				cout << "Class " << i << " identified as most trimeric class and selected as reference." << endl;
			} else {
				if ( symmetric[i] <= th ){
					latestClasses[i].setCutOffset( 1 );
				} else {
					latestClasses[i].setCutOffset( 0 );
				}
			}
		}

		// Save class selection based on correlation to latest reference
		nbfWedgedSubImage3D< Pixel > :: write( volumesFile.str().c_str(), latestClasses );
	}
	data->Delete();
}


template< class Pixel >
void nbfLoopClustering< Pixel > :: doCTFClass()
{
	stringstream volumesFile;
	volumesFile << "ctf_msa_stacks.txt";

	nbfFourierImageMetric< Pixel, 3 > fmetric( this->metric->imageFilter, this->metric->fourierFilter );
	vector< nbfWedgedSubImage3D< Pixel > > volumeList;
	nbfWedgedSubImage3D< Pixel > :: read( volumesFile.str(), volumeList );

	// Array< complex< Pixel >, 1 > C;
	// Array< Pixel, 2 > M;
	// Array< Pixel, 3 > Ri;
	// for ( int i = 0; i < volumeList.size(); i++ ){
		// fmetric.setInput1( &volumeList[i] );
		// fmetric.get2DimensionRepresentationHalf( C, M, true );
		
		// nbfPolarDomain< Pixel, 2 > polar;
		// TinyVector< Pixel , 2 > center( (M.rows()+1)/2, M.cols()-2 );
		// polar.setCenter( center );
		// polar.setMinRho( M.cols() * this->metric->fourierFilter->getBandLowCut() );
		// polar.setMaxRho( M.cols() * this->metric->fourierFilter->getBandHighCut()  );
		// // polar.setMaxRho( M.cols() );
		// polar.setResRho( 2 * M.cols() );
		// polar.setResTheta( 360 );
		// Array< Pixel, 2 > P;
		// Array< bool, 2 > B;
		// polar.cartesian2polar( M, P, B );
		// Array< Pixel, 1 > Ps( P.rows() );
		// secondIndex j;
		// Ps = sum( P( Range :: all(), Range( Ps.cols()/2, toEnd ) ), j );
		
		// if ( Ri.size() == 0 ){
			// // Ri.resize( P.rows(), P.cols(), volumeList.size() );
			// Ri.resize( Ps.rows(), volumeList.size(), 1 );
		// }
		// // Ri( Range::all(), Range::all(), i ) = P;
		// Ri( Range::all(), i, 0 ) = Ps;

		// // if ( Ri.size() == 0 ){
			// // Ri.resize( M.rows(), M.cols(), volumeList.size() );
		// // }
		// // Ri( Range::all(), Range::all(), i ) = M;
	// }
	nbfMrcWriter w;
	// w.setFileName("ctf_msa_fft.mrc");
	// w.write(Ri);
	// return;
	
	Array< Pixel, 3 > R;
	cout << "get2DRepresentations " << endl;
	this->metric->get2DRepresentations( volumesFile, R, this->useRealRepresentation, this->binFactorForClassification, this->symmetryFactor );
	
	w.setFileName("ctf_msa_data.mrc");
	w.write(R);

	cout << "doSVDDecomposition of " << R.shape() << endl;
	this->doSVDDecomposition( R, 16 );

	cout << "computeInitialReferencesKmeans " << endl;
	// classify and store in attribute: this->classes
	this->computeInitialReferencesKmeans( R, this->starting_iteration );

	// Construct class averages in parallel to minimize disk access
	// vector< nbfWedgedSubImage3D< Pixel > > volumeList;
	// nbfWedgedSubImage3D< Pixel > :: read( volumesFile.str(), volumeList );
	Array< Pixel, 3 > averages;
	
	this->metric->get2DAverages( volumesFile, this->classes, averages );
	
	// average normalization by number of images in each class
	for ( int i = 0; i < averages.depth(); i++ ){
		averages( Range :: all(), Range :: all(), i ) /= sum( this->classes( i, Range :: all() ) );
	}

	// match tomoctffind's convention
	averages.transposeSelf(secondDim,firstDim,thirdDim);
	
	w.setFileName("ctf_msa_classes.mrc");
	w.write( averages );

	vector< int > sizes;
	for ( int i = 0; i < volumeList.size(); i++ ){
		nbfMrcReader reader;
		reader.setFileName( volumeList[i].getFileName().c_str() );
		sizes.push_back( reader.getDims()[2] );
	}

	// // save class memberships to file
	// ofstream output( "classes.txt", ios::out );
	// if ( output.is_open() != 1 ){
		// output.close();
		// return;		
	// }	
	FILE * pFile;
	pFile = fopen ("ctf_msa_classes.txt","w");
	
	// output << "Class\tStack\tParticle" << endl;
	for ( int i = 0; i < this->classes.rows(); i++ ){
		for ( int j = 0; j < this->classes.cols(); j++ ){
			if ( this->classes( i, j ) > 0 ){
				int stackIndex = 0;
				int lastStackCounter = 0;
				int stackCounter = sizes[stackIndex];
				while ( j >= stackCounter ){
					stackIndex++;
					lastStackCounter = stackCounter;
					stackCounter += sizes[stackIndex];
				}
				fprintf(pFile,"%d\tStack%04d\tIndex%04d\n",i,stackIndex,j-lastStackCounter);
				// output << i << "\tStack" << stackIndex << "\tIndex" << j-lastStackCounter << endl;
			}
		}
	}

	// close file
	fclose (pFile);
	// output.close();

	// // Construct class averages in parallel to minimize disk access
	// cout << "Generating new references..." << endl;
	// for ( int i = 0; i < this->references.size(); i++ ){
		// cout << "\tGenerating reference " << i << endl;
		// this->metric->getImage( this->references[i] );
		// if ( this->metric->getMissingWedgeCompensation() == true ){
			// cout << "\tForcing class average " << i << " to have no missing-wedge" << endl;
			// // SET TO HAVE NO MISSING WEDGE
			// Array< Pixel, 3 > dummy3D( this->references[i].average.shape() );
			// dummy3D = 1;
			// this->references[i].setAccumulatedWedgeImage(dummy3D);
			// //cout << "\tUpdating missing wedge for reference " << i << endl;
			// //this->references[i].updateAccumulatedWedgeImage();
			// //// this->metric->updateAccumulatedWedgeImage( this->references[i] );
			// TinyVector< int, 2 > size;
			// size = 2 * reinterpret_cast<nbfProjectionRotationMetric3D< Pixel >*>(this->metric)->B;
			// //cout << "\tUpdating spherical missing wedge for reference " << i << " (size=" << size << ")" << endl;
			// //this->metric->updateSphericalWedgeImage( this->references[i], size );
			// Array< Pixel, 2 > dummy2D( size );
			// dummy2D = 1;
			// this->references[i].setSphericalWedgeImage( dummy2D );
		// }
	// }

	// // save state
	// stringstream currentOutputFile;

	// if ( this->starting_iteration < 10 ){
		// currentOutputFile << this->fileHeader << "_iteration_00" << this->starting_iteration;
	// } else if ( this->starting_iteration < 100 ){
		// currentOutputFile << this->fileHeader << "_iteration_0" << this->starting_iteration;
	// } else {
		// currentOutputFile << this->fileHeader << "_iteration_" << this->starting_iteration;
	// }
	// cout << "State being saved with prefix: " << currentOutputFile.str() << endl;
	// this->saveState( currentOutputFile );
}


template< class Pixel >
void nbfLoopClustering< Pixel > :: do2DClass()
{
	stringstream volumesFile;
	volumesFile << "ctf_msa_stacks.txt";

	nbfFourierImageMetric< Pixel, 3 > fmetric( this->metric->imageFilter, this->metric->fourierFilter );
	vector< nbfWedgedSubImage3D< Pixel > > volumeList;
	nbfWedgedSubImage3D< Pixel > :: read( volumesFile.str(), volumeList );

	// Array< complex< Pixel >, 1 > C;
	// Array< Pixel, 2 > M;
	// Array< Pixel, 3 > Ri;
	// for ( int i = 0; i < volumeList.size(); i++ ){
		// fmetric.setInput1( &volumeList[i] );
		// fmetric.get2DimensionRepresentationHalf( C, M, true );
		
		// nbfPolarDomain< Pixel, 2 > polar;
		// TinyVector< Pixel , 2 > center( (M.rows()+1)/2, M.cols()-2 );
		// polar.setCenter( center );
		// polar.setMinRho( M.cols() * this->metric->fourierFilter->getBandLowCut() );
		// polar.setMaxRho( M.cols() * this->metric->fourierFilter->getBandHighCut()  );
		// // polar.setMaxRho( M.cols() );
		// polar.setResRho( 2 * M.cols() );
		// polar.setResTheta( 360 );
		// Array< Pixel, 2 > P;
		// Array< bool, 2 > B;
		// polar.cartesian2polar( M, P, B );
		// Array< Pixel, 1 > Ps( P.rows() );
		// secondIndex j;
		// Ps = sum( P( Range :: all(), Range( Ps.cols()/2, toEnd ) ), j );
		
		// if ( Ri.size() == 0 ){
			// // Ri.resize( P.rows(), P.cols(), volumeList.size() );
			// Ri.resize( Ps.rows(), volumeList.size(), 1 );
		// }
		// // Ri( Range::all(), Range::all(), i ) = P;
		// Ri( Range::all(), i, 0 ) = Ps;

		// // if ( Ri.size() == 0 ){
			// // Ri.resize( M.rows(), M.cols(), volumeList.size() );
		// // }
		// // Ri( Range::all(), Range::all(), i ) = M;
	// }
	nbfMrcWriter w;
	// w.setFileName("ctf_msa_fft.mrc");
	// w.write(Ri);
	// return;
	
	Array< Pixel, 3 > R;
	cout << "get2DRepresentations " << endl;
	this->metric->get2DRepresentations( volumesFile, R, this->useRealRepresentation, this->binFactorForClassification, this->symmetryFactor );
	
	w.setFileName("ctf_msa_data.mrc");
	w.write(R);

	cout << "doSVDDecomposition of " << R.shape() << endl;
	this->doSVDDecomposition( R, 16 );

	cout << "computeInitialReferencesKmeans " << endl;
	// classify and store in attribute: this->classes
	this->computeInitialReferencesKmeans( R, this->starting_iteration );

	// Construct class averages in parallel to minimize disk access
	// vector< nbfWedgedSubImage3D< Pixel > > volumeList;
	// nbfWedgedSubImage3D< Pixel > :: read( volumesFile.str(), volumeList );
	Array< Pixel, 3 > averages;
	
	this->metric->get2DAverages( volumesFile, this->classes, averages );
	
	// average normalization by number of images in each class
	for ( int i = 0; i < averages.depth(); i++ ){
		averages( Range :: all(), Range :: all(), i ) /= sum( this->classes( i, Range :: all() ) );
	}

	// match tomoctffind's convention
	averages.transposeSelf(secondDim,firstDim,thirdDim);
	
	w.setFileName("ctf_msa_classes.mrc");
	w.write( averages );

	vector< int > sizes;
	for ( int i = 0; i < volumeList.size(); i++ ){
		nbfMrcReader reader;
		reader.setFileName( volumeList[i].getFileName().c_str() );
		sizes.push_back( reader.getDims()[2] );
	}

	// // save class memberships to file
	// ofstream output( "classes.txt", ios::out );
	// if ( output.is_open() != 1 ){
		// output.close();
		// return;		
	// }	
	FILE * pFile;
	pFile = fopen ("ctf_msa_classes.txt","w");
	
	// output << "Class\tStack\tParticle" << endl;
	for ( int i = 0; i < this->classes.rows(); i++ ){
		for ( int j = 0; j < this->classes.cols(); j++ ){
			if ( this->classes( i, j ) > 0 ){
				int stackIndex = 0;
				int lastStackCounter = 0;
				int stackCounter = sizes[stackIndex];
				while ( j >= stackCounter ){
					stackIndex++;
					lastStackCounter = stackCounter;
					stackCounter += sizes[stackIndex];
				}
				fprintf(pFile,"%d\tStack%04d\tIndex%04d\n",i,stackIndex,j-lastStackCounter);
				// output << i << "\tStack" << stackIndex << "\tIndex" << j-lastStackCounter << endl;
			}
		}
	}

	// close file
	fclose (pFile);
	// output.close();

	// // Construct class averages in parallel to minimize disk access
	// cout << "Generating new references..." << endl;
	// for ( int i = 0; i < this->references.size(); i++ ){
		// cout << "\tGenerating reference " << i << endl;
		// this->metric->getImage( this->references[i] );
		// if ( this->metric->getMissingWedgeCompensation() == true ){
			// cout << "\tForcing class average " << i << " to have no missing-wedge" << endl;
			// // SET TO HAVE NO MISSING WEDGE
			// Array< Pixel, 3 > dummy3D( this->references[i].average.shape() );
			// dummy3D = 1;
			// this->references[i].setAccumulatedWedgeImage(dummy3D);
			// //cout << "\tUpdating missing wedge for reference " << i << endl;
			// //this->references[i].updateAccumulatedWedgeImage();
			// //// this->metric->updateAccumulatedWedgeImage( this->references[i] );
			// TinyVector< int, 2 > size;
			// size = 2 * reinterpret_cast<nbfProjectionRotationMetric3D< Pixel >*>(this->metric)->B;
			// //cout << "\tUpdating spherical missing wedge for reference " << i << " (size=" << size << ")" << endl;
			// //this->metric->updateSphericalWedgeImage( this->references[i], size );
			// Array< Pixel, 2 > dummy2D( size );
			// dummy2D = 1;
			// this->references[i].setSphericalWedgeImage( dummy2D );
		// }
	// }

	// // save state
	// stringstream currentOutputFile;

	// if ( this->starting_iteration < 10 ){
		// currentOutputFile << this->fileHeader << "_iteration_00" << this->starting_iteration;
	// } else if ( this->starting_iteration < 100 ){
		// currentOutputFile << this->fileHeader << "_iteration_0" << this->starting_iteration;
	// } else {
		// currentOutputFile << this->fileHeader << "_iteration_" << this->starting_iteration;
	// }
	// cout << "State being saved with prefix: " << currentOutputFile.str() << endl;
	// this->saveState( currentOutputFile );
}


template< class Pixel >
void nbfLoopClustering< Pixel > :: doRefine( bool apply_rotational_symmetry, bool generate_volumes_for_fsc, int alignment_mode, int volumeCutOffset )
{
	this->running_mode = NBF_LOOP_CLUSTERING_REFINE;

	stringstream currentInputFile;
	if ( this->starting_iteration < 10 ){
		currentInputFile << this->fileHeader << "_iteration_00" << this->starting_iteration;
	} else if ( this->starting_iteration < 100 ){
		currentInputFile << this->fileHeader << "_iteration_0" << this->starting_iteration;
	} else {
		currentInputFile << this->fileHeader << "_iteration_" << this->starting_iteration;
	}
	this->loadState( currentInputFile );

	// bundle alignment within references
#if 0
	if ( this->refinementIterations > 0 ){
		// refine before massive multi reference alignment
		cout << "Refining within selected references:" << endl;
		this->bundleAlignment( this->references, alignment_mode );

		// generate new references
		cout << "Generating new references" << endl;
		for ( int i = 0; i < this->references.size(); i++ ){
			this->metric->getImage( this->references[i] );
		}

		// save bundle-refined references
		for ( int i = 0; i < this->references.size(); i++ ){

			stringstream fileName;

			if ( i < 10 ){
				fileName << this->fileHeader << "_iteration_00" << this->starting_iteration << "_bundle_level_" << this->hierarchicalClasses << "_average_00" << i << ".mrc";
			} else if ( i < 100 ){
				fileName << this->fileHeader << "_iteration_0" << this->starting_iteration << "_bundle_level_" << this->hierarchicalClasses << "_average_0" << i << ".mrc";
			} else {
				fileName << this->fileHeader << "_iteration_" << this->starting_iteration << "_bundle_level_" << this->hierarchicalClasses << "_average_" << i << ".mrc";
			}

			// save image data (apply symmetry if neccesary)
			vtkImageData * averageVtk = vtkImageData::New();
			this->references[i].getImage( averageVtk );
			//this->metric->imageFilter->symmetrize( averageVtk, this->symmetryFactor );
			Array< double, 3 > A;
			nbfVTKInterface :: vtkToBlitzReference( averageVtk, A );
			A *= -1;

			vtkImageCast * cast = vtkImageCast::New();
			cast->SetOutputScalarTypeToFloat();
			cast->SetInput( averageVtk );
			cast->Update();

			nbfMrcWriter mrcw;
			mrcw.setFileName( fileName.str().c_str() );
			if ( mrcw.write( cast->GetOutput() ) == false ){
				cerr << "ERROR: Failed to write " << fileName.str() << endl;
			}

			cast->Delete();
			averageVtk->Delete();

			this->precomputedReferences[i].setFileName( fileName.str().c_str() );
		}


		//// update precomputed references
		//for ( int i = 0; i < this->precomputedReferences.size(); i++ ){
		//	vtkImageData * data = vtkImageData :: New();
		//	this->references[i].getImage( data );
		//	Array< double, 3 > A;
		//	nbfVTKInterface :: vtkToBlitzReference( data, A );
		//	A *= -1;
		//	this->precomputedReferences[i].setFixedImage( data );
		//	data->Delete();

		//	TinyVector< int, 3 > currentVolumeDimensions = this->input[0].getDimensions();
		//	this->precomputedReferences[i].setDimensions( currentVolumeDimensions );
		//	this->precomputedReferences[i].setCutOffset( 0 );

		//	if ( this->references[i].wedgeImage.size() > 0 ){
		//		this->precomputedReferences[i].wedge.setImage( this->references[i].wedgeImage );
		//	}
		//	if ( this->references[i].sphericalWedgeImage.size() > 0 ){
		//		this->precomputedReferences[i].wedge.setSphericalImage( this->references[i].sphericalWedgeImage );
		//	}
		//}
	}
#endif
	// bundle alignment among references
	this->alignReferencesInBundles( apply_rotational_symmetry, alignment_mode );

	//// bundle alignment
	//if ( this->refinementIterations > 0 ){
	//	// refine before massive multi reference alignment
	//	cout << "Refining shifts within selected references:" << endl;
	//	this->bundleAlignment( this->references, alignment_mode );

	//	// generate new references
	//	cout << "Generating new references" << endl;
	//	for ( int i = 0; i < this->references.size(); i++ ){
	//		this->metric->getImage( this->references[i] );
	//	}
	//}

	//// align all references to one with most volumes
	//cout << "Aligning references to common frame..." << endl;
	//this->alignReferencesToStrongestBundleExternal( apply_rotational_symmetry );

	if ( false ){
		// apply penczek algorithm to the entire set
		cout << "Global refinement..." << endl;
		for ( int iterat = 0; iterat < this->refinementIterations; iterat++ ){

			// compute global average
			nbfWedgedAverageImage3D< Pixel > globalAverage( this->input );
			Array< Pixel, 3 > globalAlignments( this->input.size(), 17, 1 );
			globalAlignments = 1;

			// eliminate outlier volumes from global average
			for ( int c = 0; c < globalAlignments.rows(); c++ ){
				globalAlignments( c, 0, 0 ) = sum( this->classes( Range::all(), c, 0 ) );
			}

			double matrix[16];
			for ( int i = 0; i < this->input.size(); i++ ){
				vtkTransform * t = vtkTransform::New();
				this->input[i].getTransform(t);
				vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
				t->Delete();
				for ( int e = 0; e < 16; e++ ){
					globalAlignments( i, 1 + e, 0 ) = matrix[e];
				}
			}
			globalAverage.setAlignments( globalAlignments );
			this->metric->getImage( globalAverage );

			vector< nbfWedgedAverageImage3D< Pixel > > reference;
			reference.push_back( globalAverage );
			Array< Pixel, 4 > alignments( 1, this->input.size(), 17, 1 );
			alignments = -1;
			// align with refinement only
			// align with shift only (=2)
			// align with rotation but no refinement (=3)
			this->metric->getDistances( reference, this->input, alignments, 1 );

			cout << "Global score iter " << iterat << " = " << sum( alignments( Range::all(), Range::all(), 0, 0 ) ) << endl;

			// apply alignments
			for ( int i = 0; i < alignments.cols(); i++ ){
				// retrieve correction
				for ( int j = 0; j < 16; j++ ){
					matrix[j] = alignments( 0, i, 1 + j, 0 );
				}
				vtkMatrix4x4 * matriz = vtkMatrix4x4::New();
				matriz->DeepCopy(matrix);

				// retrieve original alignments
				vtkTransform * t = vtkTransform::New();
				this->input[i].getTransform(t);

				// apply correction
				vtkMatrix4x4 * mat3 = vtkMatrix4x4 :: New();
				vtkMatrix4x4 ::Multiply4x4( t->GetMatrix(), matriz, mat3 );

				this->input[i].setTransform( mat3 );
				matriz->Delete();
				mat3->Delete();
				t->Delete();
			}
		}
	}

	if ( false ){
		// re-compute averages with new transformations
		for ( int i = 0; i < this->references.size(); i++ ){
			Array< Pixel, 3 > newAlignments( this->references[i].getVolumes().size(), 17, this->references[i].weights.cols() );
			newAlignments = 0;
			newAlignments( Range::all(), 0, Range::all() ) = this->references[i].weights;

			// concatenate new transform to existing alignments
			for ( int j = 0; j < this->references[i].getVolumes().size(); j++ ){
				// retrieve alignment of current reference to common frame
				if ( this->references[i].weights(j,0) > 0 ){
					double matrix[16];
					vtkTransform * t = vtkTransform :: New();
					this->input[j].getTransform(t);
					vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
					t->Delete();
					for ( int k = 0; k < 16; k++ ){
						newAlignments( j, 1 + k, 0 ) = matrix[k];
					}
				}
			}
			// set alignments and re-compute average
			this->references[i].setAlignments( newAlignments );
			this->metric->getImage( this->references[i] );
		}
	}


	//// generate new references
	//if ( this->symmetryFactor > 1 ){
	//	cout << "Symmetrizing new references to " << this->symmetryFactor << "-fold" << endl;
	//	for ( int i = 0; i < this->references.size(); i++ ){
	//		this->symmetrize( this->references[i], this->symmetryFactor );
	//	}
	//}
		
	stringstream refinedFile;
	if ( this->starting_iteration < 10 ){
		refinedFile << this->fileHeader << "_iteration_00" << this->starting_iteration << "_refined";
	} else if ( this->starting_iteration < 100 ){
		refinedFile << this->fileHeader << "_iteration_0" << this->starting_iteration << "_refined";
	} else {
		refinedFile << this->fileHeader << "_iteration_" << this->starting_iteration << "_refined";
	}

	cout << "State being saved with prefix: " << refinedFile.str() << endl;
	this->saveState( refinedFile, generate_volumes_for_fsc );

	//#####################################################
	// SET CLASS MEMBERSHIPS SO NO INTERACTION IS NEEDED //
	//#####################################################

	vector< nbfWedgedSubImage3D< Pixel > > latestClasses;
	stringstream volumesFile;
	volumesFile << refinedFile.str() << "_averages.txt";
	nbfWedgedSubImage3D< Pixel > :: read( volumesFile.str().c_str(), latestClasses );
	// set first to be the reference
	latestClasses[0].setCutOffset( -1 );
	for ( int i = 1; i < latestClasses.size(); i++ ){
		//latestClasses[i].setCutOffset( -i-1 );
		latestClasses[i].setCutOffset(1);
	}
	// overwrite file with new memberships
	nbfWedgedSubImage3D< Pixel > :: write( volumesFile.str().c_str(), latestClasses );

	//#####################################################
	// SELECT CLASSES BY CORRELATION TO LATEST REFERENCE //
	//#####################################################

	if ( false ){
		// compute distance to previous reference and get rid of lowest correlating classes
		stringstream latestReferenceFile;
		latestReferenceFile << this->fileHeader << "_iteration_00" << this->starting_iteration - 1 << "_refined_selected_average_0.mrc";
		nbfMrcReader reader;
		vtkImageData * data = vtkImageData :: New();
		reader.setFileName( latestReferenceFile.str().c_str() );
		reader.read( data );
		// cout << "Attempting to read latest reference from " << latestReferenceFile.str().c_str() << " ...";
		if ( data->GetDimensions()[0] > 0 ){
			vector< nbfWedgedSubImage3D< Pixel > > latestClasses;
			stringstream volumesFile;
			volumesFile << refinedFile.str() << "_averages.txt";
			cout << "reading references from file " << volumesFile.str().c_str() << endl;
			nbfWedgedSubImage3D< Pixel > :: read( volumesFile.str().c_str(), latestClasses );

			vector< nbfWedgedImage3D< Pixel > * > list2;
			for ( int i = 0; i < latestClasses.size(); i++ ){
				latestClasses[i].setCutOffset(0);
				list2.push_back( &( latestClasses[i] ) );
			}

			nbfWedgedSubImage3D< Pixel > latestReference;
			latestReference = latestClasses[0];
			latestReference.setFileName( latestReferenceFile.str().c_str() );

			vector< nbfWedgedImage3D< Pixel > * > list1;
			list1.push_back( &latestReference );

			vector< TinyVector< int, 2 > > pos;
			for ( int i = 0; i < latestClasses.size(); i++ ){
				TinyVector< int, 2 > p( 0, i );
				pos.push_back(p);
			}

			Array< Pixel, 3 > refAlignments( list2.size(), 17, 1 );
			refAlignments = -1; 
			// this->metric->setTranslationSearchRestriction(0);
			this->metric->getDistances( list1, list2, pos, refAlignments, 2 );

			cout << "Distances to latest reference = " << refAlignments( Range :: all(), 0, 0 ) << endl;

			vector< Pixel > sortedDistances;
			for ( int i = 0; i < latestClasses.size(); i++ ){
				sortedDistances.push_back( refAlignments( i, 0, 0 ) );
			}
			sort( sortedDistances.begin(), sortedDistances.end() );
			Pixel th = sortedDistances[ this->distanceTopCutoffSelection * ( sortedDistances.size() - 1 ) ];

			bool firstDone = false;
			for ( int i = 0; i < latestClasses.size(); i++ ){
				if ( refAlignments( i, 0, 0 ) <= th ){
					if ( firstDone == false ){
						latestClasses[i].setCutOffset( -1 );
						firstDone = true;
					} else {
						latestClasses[i].setCutOffset( 1 );
					}
				} else {
					latestClasses[i].setCutOffset( 0 );
				}
			}
			nbfWedgedSubImage3D< Pixel > :: write( volumesFile.str().c_str(), latestClasses );
		}
		data->Delete();
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: doMra( bool reset, Pixel corr_x, Pixel corr_y, Pixel corr_z, int alignment_mode )
{
	this->running_mode = NBF_LOOP_CLUSTERING_MRA;

	vector< nbfWedgedSubImage3D< Pixel > > completeListOfVolumes;
	completeListOfVolumes = this->input;

	if ( this->usingExternalReferences == false ){
		stringstream currentFile;

		if ( this->starting_iteration < 10 ){
			currentFile << this->fileHeader << "_iteration_00" << this->starting_iteration << "_refined";
		} else if ( this->starting_iteration < 100 ){
			currentFile << this->fileHeader << "_iteration_0" << this->starting_iteration << "_refined";
		} else {
			currentFile << this->fileHeader << "_iteration_" << this->starting_iteration << "_refined";
		}
		
		this->loadState( currentFile );

		// external euler angle correction
		vtkTransform * externalTransform = vtkTransform :: New();
		externalTransform->RotateZ( - corr_z );
		externalTransform->RotateY( - corr_x );
		externalTransform->RotateZ( - corr_y );
		vtkMatrix4x4 * external = vtkMatrix4x4 :: New();
		external->DeepCopy( externalTransform->GetMatrix() );
		externalTransform->Delete();

		// assign most current alignments
		for ( int i = 0; i < this->references.size(); i++ ){
			this->references[i].rotate( external, reinterpret_cast< nbfProjectionRotationMetric3D<Pixel>* >(this->metric) );
		}
		external->Delete();

		//// symmetrize and/or generate new references
		//if ( this->symmetryFactor > 1 ){
		//	cout << "Symmetrizing new references to " << this->symmetryFactor << "-fold" << endl;
		//	for ( int i = 0; i < this->references.size(); i++ ){
		//		this->symmetrize( this->references[i], this->symmetryFactor, false );
		//	}
		//}

		// convert new references to precomputed references
		for ( int i = 0; i < this->precomputedReferences.size(); i++ ){
			// assign new rotated image and corresponding wedges
			vtkImageData * data = vtkImageData :: New(); 
			this->references[i].getImage( data );
			Array< Pixel, 3 > A;
			nbfVTKInterface :: vtkToBlitz( data, A );
			this->precomputedReferences[i].setFixedImage( A );
			// transfer missing wedge if not using symmetry
			if ( this->symmetryFactor < 2 ){
				if ( this->references[i].isWedgeUpToDate() == true ){
					this->precomputedReferences[i].wedge.setImage( this->references[i].wedgeImage );
				}
				if ( this->references[i].isSphericalWedgeUpToDate() == true ){
					this->precomputedReferences[i].wedge.setSphericalImage( this->references[i].sphericalWedgeImage );
				}
			}
		}

		// compute pre-computed references
		Array< Pixel, 4 > finalAlignmentsExternal( this->precomputedReferences.size(), 17, 1, this->classesSelected.size() );
		finalAlignmentsExternal = 0;

		// merge classes
		for ( int i = 0; i < this->classesSelected.size(); i++ ){
			for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
				int currentIndex = this->classesSelected[i][j];
				finalAlignmentsExternal( currentIndex, 0, 0, i ) = 1.0 / this->classesSelected[i].size();
				vtkMatrix4x4 * id = vtkMatrix4x4 :: New();
				id->Identity();
				double matrix[16];
				vtkMatrix4x4::DeepCopy( matrix, id );
				for ( int k = 0; k < 16; k++ ){
					finalAlignmentsExternal( currentIndex, 1 + k, 0, i ) = matrix[k];
				}
				id->Delete();
			}
		}

		nbfWedgedAverageImage3D< Pixel > newAverage( this->precomputedReferences );
		vector< nbfWedgedAverageImage3D< Pixel > > newAverages;
		for ( int i = 0; i < this->classesSelected.size(); i++ ){
			newAverages.push_back( newAverage );
		}

		// this->precomputedReferences.resize( this->classesSelected.size() );
		vtkImageData * data = vtkImageData :: New();
		for ( int i = 0; i < newAverages.size(); i++ ){
			
			Array< Pixel, 3 > A( finalAlignmentsExternal( Range :: all(), Range :: all(), Range :: all(), i ) );
			newAverages[i].setAlignments( A );

			// cannot submit to MPI because each nbfWedgedSubImage has an externally assiged volume that will not be transmited by serialize()
			newAverages[i].getImage( data );
			if ( ( this->metric->getMissingWedgeCompensation() == true ) && ( this->symmetryFactor < 2 ) ){
				newAverages[i].updateAccumulatedWedgeImage();
				if ( this->metric->getId() == NBF_IMAGE_METRIC_PROJECTION ){
					TinyVector< int, 2 > size = newAverages[i].getVolumesRO()[0].wedge.sphericalWedge.shape();
					newAverages[i].updateSphericalWedgeImage( size, reinterpret_cast< nbfProjectionRotationMetric3D< Pixel > * >(this->metric) );
				}
			}
		}
		data->Delete();

		this->references = newAverages;

		// DONE applying external euler angles to correct vertical orientation

		// APPLY SYMMETRY TO AVERAGE OF CLASS AVERAGES
		
		// symmetrize and/or generate new references
		if ( this->symmetryFactor > 1 ){
			cout << "Symmetrizing new references to " << this->symmetryFactor << "-fold" << endl;
			for ( int i = 0; i < this->references.size(); i++ ){
				this->symmetrize( this->references[i], this->symmetryFactor, false );
			}
		}

		// write new references to mrc files
		for ( int i = 0; i < this->references.size(); i++ ){
			stringstream fileName;
			fileName << currentFile.str().c_str() << "_selected_average_" << i << ".mrc";
			vtkImageData * averageVtk = vtkImageData::New();
			this->references[i].getImage( averageVtk );

			Array< double, 3 > A;
			nbfVTKInterface :: vtkToBlitzReference( averageVtk, A );
			A *= -1;

			vtkImageCast * cast = vtkImageCast::New();
			cast->SetOutputScalarTypeToFloat();
			cast->SetInput( averageVtk );
			cast->Update();

			nbfMrcWriter mrcw;
			mrcw.setFileName( fileName.str().c_str() );
			if ( mrcw.write( cast->GetOutput() ) == false ){
				cerr << "ERROR: Failed to write " << fileName.str() << endl;
			}

			// also write filtered and masked version
			this->metric->fourierFilter->execute( averageVtk );
			this->metric->imageFilter->execute( averageVtk );
			this->metric->imageFilter->symmetrize( averageVtk, this->symmetryFactor );
			vtkImageCast * castFlt = vtkImageCast::New();
			castFlt->SetOutputScalarTypeToFloat();
			castFlt->SetInput( averageVtk );
			castFlt->Update();

			stringstream fileNameFlt;
			fileNameFlt << currentFile.str().c_str() << "_selected_average_" << i << "_filtered.mrc";
			mrcw.setFileName( fileNameFlt.str().c_str() );
			if ( mrcw.write( castFlt->GetOutput() ) == false ){
				cerr << "ERROR: Failed to write " << fileName.str() << endl;
			}

			castFlt->Delete();
			cast->Delete();
			averageVtk->Delete();
		}

		// compute and save global average by merging all references
		if ( this->references.size() > 1 ){
			nbfWedgedAverageImage3D< Pixel > globalAverage( this->precomputedReferences );
			Array< Pixel, 3 > globalA( this->precomputedReferences.size(), 17, 1 );
			vtkMatrix4x4 * id = vtkMatrix4x4 :: New();
			id->Identity();
			double matrix[16];
			vtkMatrix4x4::DeepCopy( matrix, id );
			for ( int k = 0; k < 16; k++ ){
				globalA( Range :: all(), 1 + k, 0 ) = matrix[k];
			}
			id->Delete();
			globalA( Range::all(), 0, 0 ) = - 1.0 / this->precomputedReferences.size();
			
			globalAverage.setAlignments( globalA );
			vtkImageData * averageVtk = vtkImageData::New();
			globalAverage.getImage( averageVtk );
			this->metric->imageFilter->symmetrize( averageVtk, this->symmetryFactor );
		
			vtkImageCast * cast = vtkImageCast::New();
			cast->SetOutputScalarTypeToFloat();
			cast->SetInput( averageVtk );
			cast->Update();

			nbfMrcWriter mrcw;
			stringstream fileNameGlobal;
			fileNameGlobal << currentFile.str().c_str() << "_selected_global_average.mrc";
			mrcw.setFileName( fileNameGlobal.str().c_str() );
			if ( mrcw.write( cast->GetOutput() ) == false ){
				cerr << "ERROR: Failed to write " << fileNameGlobal.str() << endl;
			}

			this->metric->fourierFilter->execute( averageVtk );
			this->metric->imageFilter->execute( averageVtk );
			this->metric->imageFilter->symmetrize( averageVtk, this->symmetryFactor );

			stringstream fileNameGlobalFiltered;
			fileNameGlobalFiltered << currentFile.str().c_str() << "_selected_global_average_filtered.mrc";
			mrcw.setFileName( fileNameGlobalFiltered.str().c_str() );
			cast->Modified();
			cast->Update();
			if ( mrcw.write( cast->GetOutput() ) == false ){
				cerr << "ERROR: Failed to write " << fileNameGlobalFiltered.str() << endl;
			}

			averageVtk->Delete();
			cast->Delete();
		}


		// compute global average and save result
		if ( false ){
			Array< Pixel, 3 > aligns( this->input.size(), 17, 1 );
			aligns = 0;
			for ( int i = 0; i < this->references.size(); i++ ){
				for ( int k = 0; k < this->input.size(); k++ ){
					if ( abs( this->references[ i ].weights(k,0) ) > 0 ){
						aligns( k, 0, 0 ) = this->references[i].weights(k,0);
						aligns( k, Range(1,toEnd), 0 ) = this->references[i].multipleAlignments( k, Range :: all(), 0 );
					}
				}
			}
			nbfWedgedAverageImage3D< Pixel > total( this->input );
			total.setAlignments( aligns );
			this->metric->getImage( total );

			if ( this->symmetryFactor > 1 ){
				this->symmetrize( total, this->symmetryFactor, false );
			}

			stringstream fileName;
			fileName << currentFile.str().c_str() << "_selected_global_average.mrc";
			vtkImageData * averageVtk = vtkImageData::New();
			total.getImage( averageVtk );

			Array< double, 3 > A;
			nbfVTKInterface :: vtkToBlitzReference( averageVtk, A );
			A *= -1;

			vtkImageCast * cast = vtkImageCast::New();
			cast->SetOutputScalarTypeToFloat();
			cast->SetInput( averageVtk );
			cast->Update();

			nbfMrcWriter mrcw;
			mrcw.setFileName( fileName.str().c_str() );
			if ( mrcw.write( cast->GetOutput() ) == false ){
				cerr << "ERROR: Failed to write " << fileName.str() << endl;
			}
			cast->Delete();

			stringstream fileNameFiltered;
			fileNameFiltered << currentFile.str().c_str() << "_selected_global_average_filtered.mrc";
			total.getImage( averageVtk );

			this->metric->fourierFilter->execute( averageVtk );
			this->metric->imageFilter->execute( averageVtk );
			//this->metric->imageFilter->symmetrize( averageVtk, this->symmetryFactor );

			nbfVTKInterface :: vtkToBlitzReference( averageVtk, A );
			A *= -1;

			vtkImageCast * castF = vtkImageCast::New();
			castF->SetOutputScalarTypeToFloat();
			castF->SetInput( averageVtk );
			castF->Update();

			mrcw.setFileName( fileNameFiltered.str().c_str() );
			if ( mrcw.write( castF->GetOutput() ) == false ){
				cerr << "ERROR: Failed to write " << fileNameFiltered.str() << endl;
			}
			castF->Delete();

			averageVtk->Delete();
		}

		// switch back to the complete volume list to give all volumes a chance
		this->input.clear();
		this->input = completeListOfVolumes;
		completeListOfVolumes.clear();

		cout << "Aligning " << this->input.size() << " volumes to " << this->references.size() << " references:" << endl;
		this->alignVolumesToReferences( reset, alignment_mode );
	} else {

		// set input alignments to NULL and compute new transformations
		cout << "Aligning " << this->input.size() << " volumes to " << this->externalReferences.size() << " references:" << endl;
		this->alignVolumesToExternalReferences( alignment_mode );
	}

	//nbfMatlabWriter w;
	//w.setFileName("timer_start");
	//Array< Pixel, 1 > A;
	//w.write(A);
	//this->alignVolumesToReferences( reset );
	//w.setFileName("timer_stop");
	//w.write(A);


	// select best alignment based on overall distance

	if ( this->references.size() > 1 ){
	
		int rangeSize = 1000;
		int currentOffset = 0;
		while ( currentOffset < this->alignmentToReferences.cols() ){
	
			if ( currentOffset + rangeSize > this->alignmentToReferences.cols() ){
				rangeSize = this->alignmentToReferences.cols() - currentOffset;
			}
			
			// store each input volume with all possible transformations
			vector< nbfWedgedSubImage3D< Pixel > > duplicateVolumes;
			// for ( int i = 0; i < this->alignmentToReferences.cols(); i++ ){
			for ( int i = 0; i < rangeSize; i++ ){
				double matrix[16];
				int offset = 1;
				if ( this->alignmentToReferences.depth() > 1 + 16 ){
					offset = 2;
				} else if ( this->alignmentToReferences.depth() > 2 + 16 ){
					offset = 3;
				}
				for ( int j = 0; j < this->alignmentToReferences.rows(); j++ ){
					for ( int k = 0; k < 16; k++ ){
						matrix[k] = this->alignmentToReferences( j, currentOffset + i, k + offset, 0 );
					}
					vtkMatrix4x4 * mat = vtkMatrix4x4::New();
					mat->DeepCopy( matrix );
					this->input[ currentOffset + i ].setTransform( mat );
					duplicateVolumes.push_back( this->input[ currentOffset + i ] );
					mat->Delete();
				}
			}

			// evaluate distances
			Array< Pixel, 4 > alignments( this->references.size(), duplicateVolumes.size(), 17, 1 );
			alignments = -1;

			// evaluate correlation at position zero (we are not looking for the translation, just evaluating the distance)
			int storeTranslationTolerance = this->metric->getTranslationSearchRestriction();
			this->metric->setTranslationSearchRestriction( 0 );
			this->metric->getDistances( this->references, duplicateVolumes, alignments, 2 );
			this->metric->setTranslationSearchRestriction( storeTranslationTolerance );

			// cout << "## DEBUG ##: Done computing distances." << endl;
			
			Array< Pixel, 1 > partialSums( duplicateVolumes.size() );
			duplicateVolumes.clear();
			firstIndex i;
			secondIndex j;
			Array< Pixel, 2 > extractedDistances( alignments( Range :: all(), Range :: all(), 0, 0 ) );
			partialSums = sum( extractedDistances( j, i ), j );

			// cout << "## DEBUG ##: Reshaping partial sums in " << partialSums.shape() << endl;
			// cout << "## DEBUG ##: Reshaping partial sums to " << shape( this->alignmentToReferences.cols(), this->references.size() ) << endl;

			// Array< Pixel, 2 > reshapedPartialSums( partialSums.data(), shape( this->alignmentToReferences.cols(), this->references.size() ) );
			Array< Pixel, 2 > reshapedPartialSums( partialSums.data(), shape( rangeSize, this->references.size() ) );

			// apply alignment of best match to input volumes
			// for ( int i = 0; i < this->alignmentToReferences.cols(); i++ ){
			for ( int i = 0; i < rangeSize; i++ ){
				TinyVector< int, 1 > best = minIndex( reshapedPartialSums( i, Range::all() ) );
				this->alignmentToReferences( best[0], currentOffset + i, 0, 0 ) *= -1;
			}

			currentOffset += rangeSize;
			
			// cout << "## DEBUG ##: Done assigning optimal alignments." << endl;
		}		
	}

	// save alignments to references
	stringstream alignmentsFile;

	if ( this->starting_iteration < 10 ){
		alignmentsFile << this->fileHeader << "_iteration_00" << this->starting_iteration << "_alignmentToReferences.matlab";
	} else if ( this->starting_iteration < 100 ){
		alignmentsFile << this->fileHeader << "_iteration_0" << this->starting_iteration << "_alignmentToReferences.matlab";
	} else {
		alignmentsFile << this->fileHeader << "_iteration_" << this->starting_iteration << "_alignmentToReferences.matlab";
	}
	
	nbfMatlabWriter w;
	w.setFileName( alignmentsFile.str().c_str() );
	if ( w.write( this->alignmentToReferences ) == false ){
		cerr << "ERROR: Failed to write file " << alignmentsFile.str() << endl;
	}

	// save one txt file with the alignments to each reference
	for ( int i = 0; i < this->alignmentToReferences.rows(); i++ ){
		for ( int j = 0; j < this->alignmentToReferences.cols(); j++ ){
			
			// store distance in OFFSET field
			this->input[j].setCutOffset( fabs( this->alignmentToReferences( i, j, 0, 0 ) ) );
			
			int offset = 1;
			if ( this->alignmentToReferences.depth() > 1 + 16 ){
				offset = 2;
			} else if ( this->alignmentToReferences.depth() > 2 + 16 ){
				offset = 3;
			}

			// store alignment to current reference
			double matrix[16];
			for ( int k = 0; k < 16; k++ ){
				matrix[k] = this->alignmentToReferences( i, j, k + offset, 0 );
			}
			vtkMatrix4x4 * mat = vtkMatrix4x4::New();
			mat->DeepCopy( matrix );
			this->input[ j ].setTransform( mat );
			mat->Delete();
		}
		stringstream volumesFile;

		if ( this->starting_iteration < 10 ){
			volumesFile << this->fileHeader << "_iteration_00" << this->starting_iteration << "_alignments_to_reference_" << i << ".txt";
		} else if ( this->starting_iteration < 100 ){
			volumesFile << this->fileHeader << "_iteration_0" << this->starting_iteration << "_alignments_to_reference_" << i << ".txt";
		} else {
			volumesFile << this->fileHeader << "_iteration_" << this->starting_iteration << "_alignments_to_reference_" << i << ".txt";
		}
		nbfWedgedSubImage3D< Pixel > :: write( volumesFile.str().c_str(), this->input );
	}

}

template< class Pixel >
void nbfLoopClustering< Pixel > :: doClustering( Array< Pixel, 3 > & classification )
{
	// retrieve saved state
	stringstream currentFile;
	if ( this->running_mode == NBF_LOOP_CLUSTERING_MRA ){
		this->starting_iteration++;

		if ( this->starting_iteration < 10 ){
			currentFile << this->fileHeader << "_iteration_00" << this->starting_iteration << "_refined";
		} else if ( this->starting_iteration < 100 ){
			currentFile << this->fileHeader << "_iteration_0" << this->starting_iteration << "_refined";
		} else {
			currentFile << this->fileHeader << "_iteration_" << this->starting_iteration << "_refined";
		}
	} else {
		if ( this->running_mode == NBF_LOOP_CLUSTERING_REFINE ){
			this->starting_iteration++;
		}
		if ( this->starting_iteration < 10 ){
			currentFile << this->fileHeader << "_iteration_00" << this->starting_iteration;
		} else if ( this->starting_iteration < 100 ){
			currentFile << this->fileHeader << "_iteration_0" << this->starting_iteration;
		} else {
			currentFile << this->fileHeader << "_iteration_" << this->starting_iteration;
		}
	}

	// load previous state only if not using external references
	if ( this->usingExternalReferences == false ){
		this->loadState( currentFile );
	}

	int iterations = this->starting_iteration;

	do {

		iterations++;

		cout << "\nExecuting iteration " << this->fileHeader << "_iteration_" << iterations << endl;

		// align raw volumes to references
		if ( this->references.size() > 0 ){
			
			stringstream refinedFile;
			if ( iterations - 1 < 10 ){
				refinedFile << this->fileHeader << "_iteration_00" << iterations - 1 << "_refined";
			} else if ( iterations - 1 < 100 ){
				refinedFile << this->fileHeader << "_iteration_0" << iterations - 1 << "_refined";
			} else {
				refinedFile << this->fileHeader << "_iteration_" << iterations - 1 << "_refined";
			}
			
			if ( ( ( this->running_mode == NBF_LOOP_CLUSTERING_AUTO ) && ( this->alignmentToReferences.size() == 0 ) ) || 
				   ( this->running_mode == NBF_LOOP_CLUSTERING_REFINE ) || 
				   ( this->running_mode == NBF_LOOP_CLUSTERING_MRA ) ){

				// if not using external reference
				if ( ( this->usingExternalReferences == false ) && ( this->running_mode == NBF_LOOP_CLUSTERING_REFINE ) ){

					// BYPASS bundle alignment
					if ( false ){
						// refine before massive multi reference alignment
						cout << "Refining shifts within selected references:" << endl;
						this->bundleAlignment( this->references );

						// generate new references
						cout << "Generating new references" << endl;
						for ( int i = 0; i < this->references.size(); i++ ){
							this->metric->getImage( this->references[i] );
						}
					}

					// align all references to one with most volumes
					cout << "Aligning references to common frame..." << endl;
					this->alignReferencesToStrongest();

					if ( false ){
						// apply penczek algorithm to the entire set
						cout << "Global refinement..." << endl;
						for ( int iterat = 0; iterat < this->refinementIterations; iterat++ ){

							// compute global average
							nbfWedgedAverageImage3D< Pixel > globalAverage( this->input );
							Array< Pixel, 3 > globalAlignments( this->input.size(), 17, 1 );
							globalAlignments = 1;

							// eliminate outlier volumes from global average
							for ( int c = 0; c < globalAlignments.rows(); c++ ){
								globalAlignments( c, 0, 0 ) = sum( this->classes( Range::all(), c, 0 ) );
							}

							double matrix[16];
							for ( int i = 0; i < this->input.size(); i++ ){
								vtkTransform * t = vtkTransform::New();
								this->input[i].getTransform(t);
								vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
								t->Delete();
								for ( int e = 0; e < 16; e++ ){
									globalAlignments( i, 1 + e, 0 ) = matrix[e];
								}
							}
							globalAverage.setAlignments( globalAlignments );
							this->metric->getImage( globalAverage );

							vector< nbfWedgedAverageImage3D< Pixel > > reference;
							reference.push_back( globalAverage );
							Array< Pixel, 4 > alignments( 1, this->input.size(), 17, 1 );
							alignments = -1;
							// align with refinement only
							// align with shift only (=2)
							// align with rotation but no refinement (=3)
							this->metric->getDistances( reference, this->input, alignments, 1 );

							cout << "Global score iter " << iterat << " = " << sum( alignments( Range::all(), Range::all(), 0, 0 ) ) << endl;

							// apply alignments
							for ( int i = 0; i < alignments.cols(); i++ ){
								// retrieve correction
								for ( int j = 0; j < 16; j++ ){
									matrix[j] = alignments( 0, i, 1 + j, 0 );
								}
								vtkMatrix4x4 * matriz = vtkMatrix4x4::New();
								matriz->DeepCopy(matrix);

								// retrieve original alignments
								vtkTransform * t = vtkTransform::New();
								this->input[i].getTransform(t);

								// apply correction
								vtkMatrix4x4 * mat3 = vtkMatrix4x4 :: New();
								vtkMatrix4x4 ::Multiply4x4( t->GetMatrix(), matriz, mat3 );

								this->input[i].setTransform( mat3 );
								matriz->Delete();
								mat3->Delete();
								t->Delete();
							}
						}
					}

					if ( false ){
						// re-compute averages with new transformations
						for ( int i = 0; i < this->references.size(); i++ ){
							Array< Pixel, 3 > newAlignments( this->references[i].getVolumes().size(), 17, this->references[i].weights.cols() );
							newAlignments = 0;
							newAlignments( Range::all(), 0, Range::all() ) = this->references[i].weights;

							// concatenate new transform to existing alignments
							for ( int j = 0; j < this->references[i].getVolumes().size(); j++ ){
								// retrieve alignment of current reference to common frame
								if ( this->references[i].weights(j,0) > 0 ){
									double matrix[16];
									vtkTransform * t = vtkTransform :: New();
									this->input[j].getTransform(t);
									vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
									t->Delete();
									for ( int k = 0; k < 16; k++ ){
										newAlignments( j, 1 + k, 0 ) = matrix[k];
									}
								}
							}
							// set alignments and re-compute average
							this->references[i].setAlignments( newAlignments );
							this->metric->getImage( this->references[i] );
						}
					}

					// generate new references
					if ( this->symmetryFactor > 1 ){
						cout << "Symmetrizing new references to " << this->symmetryFactor << "-fold" << endl;
						for ( int i = 0; i < this->references.size(); i++ ){
							this->symmetrize( this->references[i], this->symmetryFactor );
						}
					}

					cout << "State being saved with prefix: " << refinedFile.str() << endl;
					this->saveState( refinedFile, true );
				} else {
					// use external reference for first iteration *only*
					this->usingExternalReferences = false;
				}

				if ( ( ( this->running_mode == NBF_LOOP_CLUSTERING_AUTO ) && ( this->refinementIterations > 0 ) ) ||
				       ( this->running_mode == NBF_LOOP_CLUSTERING_MRA ) )
				{
					// set input alignments to NULL and compute new transformations
					cout << "Aligning " << this->input.size() << " volumes to " << this->references.size() << " references:" << endl;
					this->alignVolumesToReferences();

					// save alignments to references
					stringstream alignmentsFile;
					if ( iterations - 1 < 10 ){
						alignmentsFile << this->fileHeader << "_iteration_00" << iterations - 1 << "_alignmentToReferences.matlab";
					} else if ( iterations - 1 < 100 ){
						alignmentsFile << this->fileHeader << "_iteration_0" << iterations - 1 << "_alignmentToReferences.matlab";
					} else {
						alignmentsFile << this->fileHeader << "_iteration_" << iterations - 1 << "_alignmentToReferences.matlab";
					}
					nbfMatlabWriter w;
					w.setFileName( alignmentsFile.str().c_str() );
					if ( w.write( this->alignmentToReferences ) == false ){
						cerr << "ERROR: Failed to write file " << alignmentsFile.str() << endl;
					}
				} else {
					// abort to let user re-select classes
					return;
				}
			} else {
				if ( this->usingExternalReferences == false ){
					cout << "Using pre-computed alignments to refined references from " << refinedFile.str().c_str() << endl;
					this->loadState( refinedFile );
				}
			}
		}

		// Classification
		if ( ( this->running_mode == NBF_LOOP_CLUSTERING_CLASS ) || ( this->running_mode == NBF_LOOP_CLUSTERING_AUTO ) ){

			// Recompute references
			cout << "Running classification step..." << endl;

			// load external mask (if available)
			//stringstream maskFile;
			//maskFile << this->fileHeader << "_iteration_" << iterations << "_mask.matlab";
			//this->metric->imageFilter->maskOn( maskFile );

			// set input alignments to closest reference from alignmentsToReferences
			this->computeNewReferencesHierarchical(iterations);

			// this->metric->imageFilter->maskFile.clear();

			if ( this->references.size() == 0 ){
				break;
			}

			// Construct references in parallel to minimize disk access
			cout << "Generating new references..." << endl;
			for ( int i = 0; i < this->references.size(); i++ ){
				this->metric->getImage( this->references[i] );
			}

			// save state
			stringstream currentFile;
			if ( iterations < 10 ){
				currentFile << this->fileHeader << "_iteration_00" << iterations;
			} else if ( iterations < 100 ){
				currentFile << this->fileHeader << "_iteration_0" << iterations;
			} else {
				currentFile << this->fileHeader << "_iteration_" << iterations;
			}
			cout << "State being saved with prefix: " << currentFile.str() << endl;
			this->saveState( currentFile );

			// reset alignments to references
			this->alignmentToReferences.free();
		}

	} while ( iterations - this->starting_iteration < this->maxIterations );

	// display reason for termination

	if ( iterations >= this->maxIterations ){
		cout << "Convergence reached because maximum number of iterations." << endl;
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: denoise( vector< nbfWedgedSubImage3D< Pixel > > & volumes, 
										    Array< Pixel, 3 > & alignments,
										    vector< nbfWedgedAverageImage3D< Pixel > > & denoised,
											Pixel h, Pixel th )
{
	// figure out weighting parameters before hand
	Array< Pixel, 2 > weights( volumes.size(), volumes.size() );
	weights = exp( - pow2( alignments( Range::all(), Range::all(), 0 ) / h ) );

	// limit contributions to percentage of self (if over acceptable limit)
	for ( int i = 0; i < weights.rows(); i++ ){

		weights( i, Range::all() ) = where( weights( i, Range::all() ) < .01, 0, weights( i, Range::all() ) );

		Pixel currentDecay = h;
		
		while ( sum( weights( i, Range::all() ) ) > 1 + th ){
			// compute corrected decay to make sure total contribution is below specified limit
			currentDecay -= .001;
			// recompute weights with new decay
			weights( i, Range::all() ) = exp( - pow2( alignments( i, Range::all(), 0 ) / currentDecay ) );
			weights( i, Range::all() ) = where( weights( i, Range::all() ) < .01, 0, weights( i, Range::all() ) );
		}
	}

	// build averages with computed weights
	nbfWedgedAverageImage3D< Pixel > av( this->input );

	denoised.clear();

	for ( int p = 0; p < volumes.size(); p++ ){
		Array< Pixel, 2 > alignmentToVolume( this->input.size(), 17 );
		alignmentToVolume( Range::all(), 0 ) = weights( p, Range::all() );
		alignmentToVolume( Range::all(), Range(1,toEnd) ) = this->alignments( p, Range::all(), Range(3,toEnd) );
		av.setAlignments( alignmentToVolume );
	}

	//for ( int p = 0; p < volumes.size(); p++ ){
	//	
	//	//cout << "Denoised volume " << p << ":";
	//	av.getVolumes().clear();
	//	av.getWeights().clear();

	//	for ( int i = 0; i < volumes.size(); i++ ){
	//		// add volume if non zero weight
	//		if ( weights(p,i) > 0 ){
	//			// retrieve alignment to reference
	//			double matrix[16];
	//			for ( int k = 0; k < 16; k++ ){
	//				matrix[k] = alignments( p, i, k + 3 );
	//			}

	//			if ( i != p ){
	//				vtkMatrix4x4 * mat = vtkMatrix4x4::New();
	//				mat->DeepCopy( matrix );
	//				volumes[ i ].setTransform( mat );
	//				mat->Delete();
	//			} else {
	//				volumes[i].setTransform( (vtkTransform*)NULL );
	//				av.useGivenWedge( av.getVolumes().size() );
	//			}

	//			av.getVolumes().push_back( volumes[i] );
	//			av.getWeights().push_back( weights(p,i) );
	//		}
	//	}
	//	//cout << endl;
	//	if ( av.getVolumes().size() > 0 ){
	//		denoised.push_back( av );
	//	} else {
	//		cerr << "ERROR - Average does not contain volumes. __FILE__ , __LINE__" << endl;
	//	}
	//}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: alignVolumesToReferences( bool reset, int alignment_mode )
{
	// initialize new alignment matrix of references to volumes
	this->alignmentToReferences.resize( this->references.size(), this->input.size(), 17, this->metric->getNumberOfCandidates() );
	this->alignmentToReferences = -1;

	if ( reset == false ){
		cout << "\tUsing most recent volume alignments for seeding current MRA round." << endl;
		//this->metric->getDistances( this->references, this->input, this->alignmentToReferences, 1 );
	} else {
		// reset transforms before computing distances
		for ( int i = 0; i < this->input.size(); i++ ){
			this->input[i].setTransform( (vtkTransform*)NULL );
		}
		cout << " Reseting current volume alignments for seeding current MRA round." << endl;
		//this->metric->getDistances( this->references, this->input, this->alignmentToReferences, 0 );
	}
	// full search
	this->metric->getDistances( this->references, this->input, this->alignmentToReferences, alignment_mode );
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: alignVolumesToExternalReferences( int alignment_mode )
{
	// initialize new alignment matrix of references to volumes
	this->alignmentToReferences.resize( this->externalReferences.size(), this->input.size(), 17, this->metric->getNumberOfCandidates() );
	this->alignmentToReferences = -1;

	if ( this->refinementOnly == true ){
		cout << "WARNING: Volume alignments are NOT properly set. In " << __FILE__ << ":" << __LINE__ << endl;
		this->metric->getDistances( this->externalReferences, this->input, this->alignmentToReferences, 1 );
	} else {
		// reset transforms before computing distances
		for ( int i = 0; i < this->input.size(); i++ ){
			this->input[i].setTransform( (vtkTransform*)NULL );
		}
		// full search
		this->metric->getDistances( this->externalReferences, this->input, this->alignmentToReferences, alignment_mode );
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: computeInitialReferences( Array< Pixel, 2 > & T, int iterations )
{
	// reset current references
	this->references.clear();

	// compute initial classes with hierarchical clustering
	nbfHierarchicalClustering< Pixel > hierarch;
	hierarch.setInput( this->input );
	hierarch.setMetric( this->metric, this->alignments );

	//vector< Pixel > sortedDistances;
	//for ( int i = 0; i < this->alignments.rows(); i++ ){
	//	for ( int j = 0; j < this->alignments.cols(); j++ ){
	//		sortedDistances.push_back( this->alignments( i, j, 0, 0 ) );
	//	}
	//}
	//sort( sortedDistances.begin(), sortedDistances.end() );
	//Pixel cutoff = sortedDistances[ floor( ( sortedDistances.size() - 1 ) / 2.0 ) ] * this->hierarchicalCutoff;

	//Pixel maxTh = max( this->alignments( Range::all(), Range::all(), 0, 0 ) );
	//Pixel minTh = min( where( this->alignments( Range::all(), Range::all(), 0, 0 ) == 0, numeric_limits< Pixel > :: max(), this->alignments( Range::all(), Range::all(), 0, 0 ) ) );
	//Pixel cutoff = ( maxTh - minTh ) * this->hierarchicalCutoff + minTh;

	//cout << "Using hierarchical cutoff = " << cutoff << " in range [" << minTh << "," << maxTh << "]." << endl;
	//hierarch.setMaxDistance( cutoff );
	hierarch.setMaxClusters( this->hierarchicalClasses );

	//hierarch.setMaxDistance( this->hierarchicalCutoff );
	hierarch.setMinOverlap( this->hierarchicalMinOverlap );
	hierarch.setMinElementNumber( this->minVolumeNumber );

	// assign already computed tree
	hierarch.currentTree.resize( T.shape() );
	hierarch.currentTree = T;

	hierarch.execute( this->classes );

	//// Post-processor
	//this->postProcessor( R, this->classes );

	T.resize( hierarch.currentTree.shape() );
	T = hierarch.currentTree;

	this->referencesIndexes = hierarch.classIndexes;
	
	// build class averages
	for ( int i = 0; i < this->classes.rows(); i++ ){
		nbfWedgedAverageImage3D< Pixel > average( this->input );

		// assign alignments
		Array< Pixel, 3 > referenceAlignments( this->input.size(), 17, 1 );
		referenceAlignments = 0;

		// set new weights
		referenceAlignments( Range::all(), 0, 0 ) = this->classes( i, Range::all(), 0 );

		for ( int j = 0; j < this->classes.cols(); j++ ){
			if ( this->classes( i, j, 0 ) > 0 ){
				// set alignment to best match
				vtkTransform * t = vtkTransform::New();
				this->input[j].getTransform(t);
				double matrix[16];
				vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
				t->Delete();
				for ( int e = 0; e < 16; e++ ){
					referenceAlignments( j, 1 + e, 0 ) = matrix[e];
				}
			}
		}

		average.setAlignments( referenceAlignments );

		this->references.push_back( average );
	}

	// build hierarchical tree if so required
	if ( this->hierarchicalToBuildEntireTree > 0 ){
		for ( int i = 0; i < hierarch.cumulativeClasses.size(); i++ ){
			stringstream currFile;
			if ( iterations < 10 ){
				currFile << this->fileHeader << "_iteration_00" << iterations << "_average_";
			} else if ( iterations < 100 ){
				currFile << this->fileHeader << "_iteration_0" << iterations << "_average_";
			} else {
				currFile << this->fileHeader << "_iteration_" << iterations << "_average_";
			}
			if ( hierarch.cumulativeIndexes[i][0] < 10 )
				currFile << "0";
			currFile << hierarch.cumulativeIndexes[i][0] << "_from_";
			if ( hierarch.cumulativeIndexes[i][1] < 10 )
				currFile << "0";
			currFile << hierarch.cumulativeIndexes[i][1] << "_";
			if ( hierarch.cumulativeIndexes[i][2] < 10 )
				currFile << "0";
			currFile << hierarch.cumulativeIndexes[i][2] << ".mrc";

			// build current reference

			nbfWedgedAverageImage3D< Pixel > average( this->input );
			Array< Pixel, 3 > referenceAlignments( this->input.size(), 17, 1 );
			referenceAlignments = 0;

			// set new weights
			for ( int j = 0; j < hierarch.cumulativeClasses[i].size(); j++ ){
				int currIndex = hierarch.cumulativeClasses[i][j];
				referenceAlignments( currIndex, 0, 0 ) = 1;
				// set alignment to the one stored in the input
				vtkTransform * t = vtkTransform::New();
				this->input[ currIndex ].getTransform(t);
				double matrix[16];
				vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
				t->Delete();
				for ( int e = 0; e < 16; e++ ){
					referenceAlignments( currIndex, 1 + e, 0 ) = matrix[e];
				}
			}
			
			average.setAlignments( referenceAlignments );
			this->metric->getImage( average );

			vtkImageData * averageVtk = vtkImageData::New();
			average.getImage( averageVtk );

			//// change geometry to make it compatible with mrc file
			//Array< double, 3 > A, B;
			//nbfVTKInterface::vtkToBlitzReference( averageVtk, A );

			//A.transposeSelf(thirdDim,secondDim,firstDim);
			//B.resize( A.shape() );
			//B = A.reverse(secondDim);
			//nbfVTKInterface::blitzToVtk(B,averageVtk);

			vtkImageCast * cast = vtkImageCast::New();
			cast->SetOutputScalarTypeToFloat();
			cast->SetInput( averageVtk );
			cast->Update();

			nbfMrcWriter mrcw;
			mrcw.setFileName( currFile.str().c_str() );
			if ( mrcw.write( cast->GetOutput() ) == false ){
				cerr << "Failed to write " << currFile.str() << endl;
			}
			cast->Delete();
			averageVtk->Delete();
		}
	}
}


template< class Pixel >
void nbfLoopClustering< Pixel > :: computeInitialReferencesKmeans( Array< Pixel, 3 > & R, int iterations )
{
	// compute initial classes with kmeans clustering
	KMterm  term(100, 0, 0, 0,   // run for 100 stages
		0.10, 0.10, 3,             // other typical parameter values 
		0.50, 10, 0.95);

	// DO NOT leave out first coordinate (we are substracting the mean now)
	Array< Pixel, 2 > globalA( R( Range :: all(), Range :: all(), 0 ) );
	//Array< Pixel, 2 > globalA( R( Range :: all(), Range(1,toEnd), 0 ) );

	//// FIRST ELIMINATE POINTS THAT ARE FAR AWAY FROM THE GLOBAL AVERAGE

	//// compute global average
	//Array< Pixel, 1 > globalAverage( globalA.cols() );
	//firstIndex i;
	//secondIndex j;
	//globalAverage = mean( globalA(j,i), j );

	//// compute cutoff distance to global average
	//Array< Pixel, 1 > distances( globalA.rows() );
	//vector< Pixel > sortedDistances;
	//for ( int i = 0; i < globalA.rows(); i++ ){
	//	distances(i) = sum( pow2( globalA( i, Range::all() ) - globalAverage ) );
	//	sortedDistances.push_back( distances(i) );
	//}
	//sort( sortedDistances.begin(), sortedDistances.end() );
	//Pixel maxDistance = sortedDistances[ floor( ( sortedDistances.size() - 1 ) * this->distanceTopCutoffPre ) ];

	////nbfMatlabWriter w;
	////w.setFileName("p.matlab");
	////w.write(distances);

	//// OVERRIDE
	//maxDistance = sortedDistances[ floor( sortedDistances.size() - 1.0 ) ];

	//// number of volumes submitted to classification
	//int effectiveVolumesToClassify = sum( where( distances <= maxDistance, 1, 0 ) );

	//if ( effectiveVolumesToClassify < distances.rows() ) {
	//	cout << "Eliminating " << distances.rows() - effectiveVolumesToClassify << " outlier volumes corresponding to " << 100*(1 - this->distanceTopCutoffPre) << "% of total." << endl;
	//}

	// outliers are set to 0
	secondIndex jx;
	Array< Pixel, 1 > distances( globalA.rows() );
	distances = sum( fabs(globalA), jx );

	//cout << globalA << endl;
	//cout << distances << endl;

	// number of volumes submitted to classification
	int effectiveVolumesToClassify = sum( where( distances != 0, 1, 0 ) );

	// extract points and store index information
	Array< Pixel, 2 > A( effectiveVolumesToClassify, globalA.cols() );
	Array< int, 1 > indexes( A.rows() );
	int count = 0;
	for ( int i = 0; i < globalA.rows(); i++ ){
		//if ( distances(i) <= maxDistance ){
		if ( distances(i) != 0 ){
			A( count, Range :: all() ) = globalA( i, Range :: all() );
			indexes(count++) = i;
		}
	}

	// KMEANS CLUSTERING

	// number of centers
	int	k = hierarchicalClasses;
	int	dim	= A.cols();		// dimension
	int	nPts = A.rows();		// number of data points

	if ( nPts < k ){
		cerr << "Number of classes is greater than number of volumes?. Skipping." << endl;;
		return;
	}

	KMdata dataPts( dim, nPts );			// allocate data storage

	for (int d = 0; d < dim; d++ ) {
		for ( int k = 0; k < nPts; k++ ) {
			dataPts.getPts()[k][d] = A(k,d);
		}
	}
	dataPts.buildKcTree();			// build filtering structure
	KMfilterCenters ctrs( k, dataPts );		// allocate centers

	// run the algorithm
	// KMlocalLloyds       kmAlg(ctrs, term);	// repeated Lloyd's
	// KMlocalSwap      kmAlg(ctrs, term);	// Swap heuristic
	// KMlocalEZ_Hybrid kmAlg(ctrs, term);	// EZ-Hybrid heuristic
	KMlocalHybrid    kmAlg(ctrs, term);	// Hybrid heuristic
	ctrs = kmAlg.execute();			// execute
	cout << "Number of stages: " << kmAlg.getTotalStages() << "\n";
	cout << "Average distortion: " << ctrs.getDist()/nPts << "\n";
	ctrs.print();				// print final centers

	KMctrIdx * closeCtr = new KMctrIdx[dataPts.getNPts()];
	double * sqDist = new double[dataPts.getNPts()];
	ctrs.getAssignments( closeCtr, sqDist );

	//*kmOut	<< "  (Cluster assignments:\n"
	//	<< "    Point  Center  Squared Dist\n"
	//	<< "    -----  ------  ------------\n";
	//for (int i = 0; i < dataPts.getNPts(); i++) {
	//	*kmOut	<< "   " << setw(5) << i
	//		<< "   " << setw(5) << closeCtr[i]
	//		<< "   " << setw(10) << sqDist[i]
	//		<< "\n";
	//}
	//*kmOut << "  )\n";

	// search for distance cutoff within each class
	vector< vector< double > > sortedSqDistances( k );
	for ( int i = 0; i < dataPts.getNPts(); i++ ){
		sortedSqDistances[ closeCtr[i] ].push_back( sqDist[i] );
	}
	vector< double > cutOffDistances;
	for ( int i = 0; i < k; i++ ){
		sort( sortedSqDistances[i].begin(), sortedSqDistances[i].end() );
		cutOffDistances.push_back( sortedSqDistances[i][ floor( ( sortedSqDistances[i].size() - 1 ) * this->distanceTopCutoff ) ] );
	}
	
	// transfer classification result to class attribute
	this->classes.resize( hierarchicalClasses, globalA.rows(), 1 );
	this->classes = 0;
	Array< int, 1 > outliers( hierarchicalClasses );
	outliers = 0;
	for ( int i = 0; i < dataPts.getNPts(); i++ ){
		if ( sqDist[i] <= cutOffDistances[ closeCtr[i] ] ){
			this->classes( closeCtr[i], indexes(i), 0 ) = 1;
		} else {
			outliers(closeCtr[i]) += 1;
		}
	}

	cout << "Number of high variance volumes eliminated (per class):\n" << outliers << endl;

	delete [] closeCtr;
	delete [] sqDist;

	// sort classes by decreasing number of volumes
	vector< int > classSizes;
	for ( int i = 0; i < this->classes.rows(); i++ ){
		int numVolumes = sum(  this->classes( i, Range::all(), 0 ) );
		classSizes.push_back( numVolumes );
	}
	sort( classSizes.begin(), classSizes.end() );
	reverse( classSizes.begin(), classSizes.end() );
	
	cout << "Classes sorted by decreasing number of members: " << endl;
	
	vector< bool > indexesBefore;
	for ( int i = 0; i < this->classes.rows(); i++ ){
		indexesBefore.push_back( true );
	}
	
	Array< Pixel, 3 > sortedClasses( this->classes.shape() );
	for ( int i = 0; i < this->classes.rows(); i++ ){
		int currentClass;
		// search for class with given number of volumes
		for ( int j = 0; j < this->classes.rows(); j++ ){
			if ( ( indexesBefore[j] == true ) && ( sum(  this->classes( j, Range::all(), 0 ) ) == classSizes[i] ) ){
				currentClass = j;
				indexesBefore[j] = false; // mark as already assigned
				break;
			}
		}
		sortedClasses( i, Range :: all(), Range :: all() ) = this->classes( currentClass, Range :: all(), Range :: all() );
	}
	this->classes = sortedClasses;	

	// reset current references
	this->references.clear();

	// build class averages
	if ( this->input.size() == this->classes.cols() ){
		for ( int i = 0; i < this->classes.rows(); i++ ){
			nbfWedgedAverageImage3D< Pixel > average( this->input );

			// assign alignments
			Array< Pixel, 3 > referenceAlignments( this->input.size(), 17, 1 );
			referenceAlignments = 0;

			// set new weights
			referenceAlignments( Range::all(), 0, 0 ) = this->classes( i, Range::all(), 0 );

			for ( int j = 0; j < this->classes.cols(); j++ ){
				if ( this->classes( i, j, 0 ) > 0 ){
					// set alignment to best match
					vtkTransform * t = vtkTransform::New();
					this->input[j].getTransform(t);
					double matrix[16];
					vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
					t->Delete();
					for ( int e = 0; e < 16; e++ ){
						referenceAlignments( j, 1 + e, 0 ) = matrix[e];
					}
				}
			}

			average.setAlignments( referenceAlignments );

			this->references.push_back( average );
		}
	}

	for ( int i = 0; i <  this->classes.rows(); i++ ){
		cout << "Class " << i << "(" << sum(  this->classes( i, Range::all(), 0 ) ) << "): " ;
		for ( int j = 0; j <  this->classes.cols(); j++ ){
			if (  this->classes( i, j, 0 ) > 0 ){
				cout << j << " " ;
			}
		}
		cout << endl;
	}

}


template< class Pixel >
void nbfLoopClustering< Pixel > :: computeNewReferences( int iteration )
{
	// reset current classification
	this->classes.resize( this->alignmentToReferences.rows(), this->alignmentToReferences.cols(), this->alignmentToReferences.extent(fourthDim) );

	// normalize distance ranges
	Array< Pixel, 2 > normalizedDistances( this->alignmentToReferences.rows(), this->alignmentToReferences.cols() );
	normalizedDistances = this->alignmentToReferences( Range::all(), Range::all(), 0, 0 );
	// Pixel xmin = min( normalizedDistances );
	// Pixel xmax = max( normalizedDistances );
	// cout << "xmin=" << xmin << ", xmax=" << xmax << endl;
	// normalizedDistances = ( normalizedDistances - xmin ) * xmax / ( xmax - xmin );

	// estimate decay as average of radius (defined as most distant element within a class)
	
	// sigmoid decay
	Pixel decay = .025;

	//Array< Pixel, 1 > R( dummy.rows() );
	this->referenceRadii.resize( normalizedDistances.rows() );

	firstIndex h; secondIndex k;
	Array< Pixel, 1 > A( normalizedDistances.cols() );
	A = min( normalizedDistances(k,h), k );
	for ( int i = 0; i < normalizedDistances.rows(); i++ ){
		
		// look for furthest element assigned to this class
		Pixel selfBased = max( where( A == normalizedDistances( i, Range::all() ), normalizedDistances( i, Range::all() ), - numeric_limits< Pixel > :: max() ) );
		// look for closest element assigned to another class
		Pixel otherBased = min( where( A == normalizedDistances( i, Range::all() ), numeric_limits< Pixel > :: max(), normalizedDistances( i, Range::all() ) ) );
		// compute radii as minima between self and other based distances
		this->referenceRadii(i) = min( selfBased, otherBased );

		// cout << "old radii(" << i << ") = " << this->referenceRadii(i) << endl;

		// compute CUTOFF as median of distances to reference
		vector< Pixel > medianVector;
		for ( int j = 0; j < normalizedDistances.cols(); j++ ){
			medianVector.push_back( normalizedDistances( i, j ) );
		}
		sort( medianVector.begin(), medianVector.end() );
		this->referenceRadii(i) = 1.0 * medianVector[ floor( medianVector.size() / 2.0 ) ];

		this->referenceRadii(i) = max( normalizedDistances( i, Range::all() ) ) + decay * log( .01 );

		cout << "Exponential decay (" << i << ") = " << this->referenceRadii(i) << endl;
	}

	// estimate class radii based on number of contributing volumes

	// compute weights from exponential decay (apply scaling?)
	for ( int i = 0; i < this->classes.rows(); i++ ){

		// set sigmoid function
		this->classes( i, Range::all(), Range::all() ) = ( 1 + 2.0 * exp( - ( this->alignmentToReferences( i, Range::all(), 0, Range::all() ) - this->referenceRadii(i) ) / decay ) ) /
			                               ( 1 + 1.0 * exp( - ( this->alignmentToReferences( i, Range::all(), 0, Range::all() ) - this->referenceRadii(i) ) / decay ) ) - 1.0;
		
		//this->classes( i, Range::all() ) = 0;
		//this->classes( i, i ) = 1;

		//// set gaussian width
		//this->classes( i, Range::all() ) = exp( - pow2( normalizedDistances( i, Range::all() ) / ( this->referenceRadii(i) / 1.0 ) ) );
		//
		//// truncate weights if too small wrt principal component
		//this->classes( i, Range::all() ) = where( this->classes( i, Range::all() ) > max( this->classes( i, Range::all() ) ) / 20.0, this->classes( i, Range::all() ), 0 );
		// truncate weights if too small
		this->classes( i, Range::all(), Range::all() ) = where( this->classes( i, Range::all(), Range::all() ) > .01, this->classes( i, Range::all(), Range::all() ), 0 );
	}
	// cout << "Weights after truncation = " << this->classes << endl;
	
	//// compute unnormalized radii for reference merging
	//Array< Pixel, 2 > dummy( this->alignmentToReferences( Range::all(), Range::all(), 0 ) );
	//A = min( dummy(k,h), k );
	//for ( int i = 0; i < dummy.rows(); i++ ){
	//	this->referenceRadii(i) = min( where( A == this->alignmentToReferences( i, Range::all(), 0 ), numeric_limits< Pixel > :: max(), this->alignmentToReferences( i, Range::all(), 0 ) ) );
	//}

	//cout << "Adjusted classes radii = " << this->referenceRadii << endl;

	// Rebuild references
	this->references.clear();
	for ( int i = 0; i < this->classes.rows(); i++ ){

		nbfWedgedAverageImage3D< Pixel > wedgedAverage( this->input );

		Array< Pixel, 3 > alignments2( this->input.size(), 17, this->alignmentToReferences.extent(fourthDim) );
		// set weights as first components
		alignments2( Range::all(), 0, Range::all() ) = this->classes( i, Range::all(), Range::all() );
		// set alignments
		alignments2( Range::all(), Range(1,toEnd), Range::all() ) = this->alignmentToReferences( i, Range::all(), Range(3,toEnd), Range::all() );
		// set to average
		wedgedAverage.setAlignments( alignments2 );

		this->references.push_back( wedgedAverage );
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: computeNewReferencesHierarchical( int iteration )
{
	// apply alignment of best match and eliminate outliers
	this->asignBestAlignmentAndClean();

	// save current alignments into file
	stringstream volumesFile;
	if ( iteration < 10 ){
		volumesFile << this->fileHeader << "_iteration_00" << iteration << "_volumes_tmp.txt";
	} else if ( iteration < 100 ){
		volumesFile << this->fileHeader << "_iteration_0" << iteration << "_volumes_tmp.txt";
	} else {
		volumesFile << this->fileHeader << "_iteration_" << iteration << "_volumes_tmp.txt";
	}
	nbfWedgedSubImage3D< Pixel > :: write( volumesFile.str(), this->input );

	nbfMatlabReader mreader;
	nbfMatlabWriter w;

	// compute lower dimensional representations
	Array< Pixel, 3 > R;
	this->metric->getRepresentations( volumesFile, R, 0, this->binFactorForClassification );

	stringstream repFourierFile;
	if ( iteration < 10 ){
		repFourierFile << this->fileHeader << "_iteration_00" << iteration << "_reducedFourierRepresentation.matlab";
	} else if ( iteration < 100 ){
		repFourierFile << this->fileHeader << "_iteration_0" << iteration << "_reducedFourierRepresentation.matlab";
	} else {
		repFourierFile << this->fileHeader << "_iteration_" << iteration << "_reducedFourierRepresentation.matlab";
	}
	w.setFileName( repFourierFile.str().c_str() );
	if ( w.write( R ) == false ){
		cerr << "ERROR: Failed to write " << repFourierFile.str() << endl;
	}

	// distances in lower dimensional space
	Array< Pixel, 2 > D;

	this->metric->getDistances( repFourierFile, D );

	// initialize distance alignments properly
	//this->alignments.resize( D.rows(), D.cols(), 19, 1 );
	// only use distance and overlaps
	this->alignments.resize( D.rows(), D.cols(), 2, 1 );
	
	// set distance matrix in first component
	this->alignments( Range::all(), Range::all(), 0, 0 ) = D;
	//this->alignments( Range::all(), Range::all(), 0, 0 ) = D( Range::all(), Range::all(), 0, 0 );

	// set overlaps and scalings to unit
	this->alignments( Range::all(), Range::all(), Range(1,toEnd), 0 ) = 1;

	stringstream treeFile;
	if ( iteration < 10 ){
		treeFile << this->fileHeader << "_iteration_00" << iteration << "_SVD_distancesTree.matlab";
	} else if ( iteration < 100 ){
		treeFile << this->fileHeader << "_iteration_0" << iteration << "_SVD_distancesTree.matlab";
	} else {
		treeFile << this->fileHeader << "_iteration_" << iteration << "_SVD_distancesTree.matlab";
	}

	Array< Pixel, 2 > T;
	
	mreader.setFileName( treeFile.str().c_str() );
	mreader.read( T );

	bool existingTree = false;

	if ( T.size() > 0 ){
		cout << "Using pre-computed distance tree from file:\n\t" << treeFile.str().c_str() << "\n\tSize = " << T.shape() << endl;
		existingTree = true;
	}

	cout << "Computing hierarchical clustering tree..." << endl;
	this->computeInitialReferences(T,iteration);
	
	if ( existingTree == false ){
		w.setFileName( treeFile.str().c_str() );
		if ( w.write( T ) == false ){
			cerr << "ERROR: Failed to write " << treeFile.str() << endl;
		}
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: asignBestAlignmentAndClean()
{
	// apply alignment of best match to input volumes
	vector< Pixel > distances;
	for ( int i = 0; i < this->alignmentToReferences.cols(); i++ ){
		TinyVector< int, 1 > best = minIndex( this->alignmentToReferences( Range::all(), i, 0, 0 ) );
		double matrix[16];
		int offset = 1;
		if ( this->alignmentToReferences.depth() > 1 + 16 ){
			offset = 2;
		} else if ( this->alignmentToReferences.depth() > 2 + 16 ){
			offset = 3;
		}
		for ( int k = 0; k < 16; k++ ){
			matrix[k] = this->alignmentToReferences( best[0], i, k + offset, 0 );
		}
		vtkMatrix4x4 * mat = vtkMatrix4x4::New();
		mat->DeepCopy( matrix );
		this->input[ i ].setTransform( mat );
		mat->Delete();
		distances.push_back( fabs( this->alignmentToReferences( best[0], i, 0, 0 ) ) ); 
	}

	// eliminate outliers: overwrite this->input with pruned list of volumes
	if ( ( distances.size() > 0 ) && ( this->distanceTopCutoffPre < 1 ) ){
		vector< Pixel > sortedDistances( distances.size() );
		sortedDistances = distances;
		sort( sortedDistances.begin(), sortedDistances.end() );
		Pixel cutoffDistance = sortedDistances[ this->distanceTopCutoffPre * ( sortedDistances.size() - 1 ) ];

		// transfer classification result to class attribute
		vector< nbfWedgedSubImage3D< Pixel > > cleanVolumeList;
		for ( int i = 0; i < this->input.size(); i++ ){
			if ( distances[i] < cutoffDistance ){
				cleanVolumeList.push_back( this->input[i] );
			}
		}

		cout << "Eliminating " << this->input.size() - cleanVolumeList.size() << " outlier volumes corresponding to " << 100*(1 - this->distanceTopCutoffPre) << "% of total." << endl;

		this->input.clear();
		this->input = cleanVolumeList;

		cleanVolumeList.clear();
	}
}


template< class Pixel >
void nbfLoopClustering< Pixel > :: computeNewReferencesMSA( int iteration )
{
	// apply alignment of best match and eliminate outliers
	this->asignBestAlignmentAndClean();

	if ( this->hierarchicalClasses > 1 ){

		// save current alignments into file
		stringstream volumesFile;
		if ( iteration < 10 ){
			volumesFile << this->fileHeader << "_iteration_00" << iteration << "_volumes_tmp.txt";
		} else if ( iteration < 100 ){
			volumesFile << this->fileHeader << "_iteration_0" << iteration << "_volumes_tmp.txt";
		} else {
			volumesFile << this->fileHeader << "_iteration_" << iteration << "_volumes_tmp.txt";
		}
		nbfWedgedSubImage3D< Pixel > :: write( volumesFile.str(), this->input );

		nbfMatlabReader mreader;
		nbfMatlabWriter w;

		// distances in lower dimensional space
		Array< Pixel, 2 > D;
		//stringstream distancesFileSVD;
		//distancesFileSVD << this->fileHeader << "_iteration_" << iteration << "_SVD_distances.matlab";
		//mreader.setFileName( distancesFileSVD.str().c_str() );
		//mreader.read( D );

		stringstream repSVDFile;
		if ( iteration < 10 ){
			repSVDFile << this->fileHeader << "_iteration_00" << iteration << "_LowDimensionalRepresentation.matlab";
		} else if ( iteration < 100 ){
			repSVDFile << this->fileHeader << "_iteration_0" << iteration << "_LowDimensionalRepresentation.matlab";
		} else {
			repSVDFile << this->fileHeader << "_iteration_" << iteration << "_LowDimensionalRepresentation.matlab";
		}
		Array< Pixel, 3 > R;
		mreader.setFileName( repSVDFile.str().c_str() );
		mreader.read( R );

		// is SVD representation available?
		if ( R.size() > 0 ){
			cout << "Using pre-computed lower dimensional coordinates from file:\n\t" << repSVDFile.str().c_str() << "\n\tSize = " << R.shape() << endl;
		} else {

			//stringstream repFourierFile;
			//repFourierFile << this->fileHeader << "_iteration_" << iteration << "_fourierRepresentation.matlab";

			//if ( R.size() > 0 ){
			//	cout << "Using low dimensional representation from file:\n\t" << repSVDFile.str().c_str() << "\n\tSize = " << R.shape() << endl;
			//} else {

			// need to compute SVD representation

			stringstream repFile;
			if ( iteration < 10 ){
				repFile << this->fileHeader << "_iteration_00" << iteration << "_reducedRepresentation.matlab";
			} else if ( iteration < 100 ){
				repFile << this->fileHeader << "_iteration_0" << iteration << "_reducedRepresentation.matlab";
			} else {
				repFile << this->fileHeader << "_iteration_" << iteration << "_reducedRepresentation.matlab";
			}

			mreader.setFileName( repFile.str().c_str() );
			mreader.read( R );

			if ( R.size() > 0 ){
				cout << "Using reduced representation from file:\n\t" << repFile.str().c_str() << "\n\tSize = " << R.shape() << endl;
			} else {

				cout << "Computing reduced representations..." << endl;
				if ( this->useRealRepresentation > 0 ){
					this->metric->getRepresentations( volumesFile, R, this->useRealRepresentation, this->binFactorForClassification, this->symmetryFactor );

					//// extract mean
					//Array< Pixel, 2 > sR( R( Range :: all(), Range :: all(), 0 ) );
					//Array< Pixel, 1 > globalAverage( sR.cols() );
					//firstIndex ix;
					//secondIndex jx;
					//globalAverage = mean( sR(jx,ix), jx );

					//// extract global average
					//for ( int i = 0; i < sR.rows(); i++ ){
					//	sR( i, Range :: all() ) -= globalAverage;
					//}

				} else {

					// data imputation

					stringstream imputedRepFile;
					if ( iteration < 10 ){
						imputedRepFile << this->fileHeader << "_iteration_00" << iteration << "_imputedFourierRepresentation.matlab";
					} else if ( iteration < 100 ){
						imputedRepFile << this->fileHeader << "_iteration_0" << iteration << "_imputedFourierRepresentation.matlab";
					} else {
						imputedRepFile << this->fileHeader << "_iteration_" << iteration << "_imputedFourierRepresentation.matlab";
					}

					Array< Pixel, 3 > Rimputed;

					mreader.setFileName( imputedRepFile.str().c_str() );
					mreader.read( Rimputed );

					if ( Rimputed.size() > 0 ){
						cout << "Using imputed representation from file:\n\t" << imputedRepFile.str().c_str() << "\n\tSize = " << Rimputed.shape() << endl;
					} else {

						stringstream repFourierFile;
						if ( iteration < 10 ){
							repFourierFile << this->fileHeader << "_iteration_00" << iteration << "_reducedFourierRepresentation.matlab";
						} else if ( iteration < 100 ){
							repFourierFile << this->fileHeader << "_iteration_0" << iteration << "_reducedFourierRepresentation.matlab";
						} else {
							repFourierFile << this->fileHeader << "_iteration_" << iteration << "_reducedFourierRepresentation.matlab";
						}

						// use distances taking the missing wedge into account
						mreader.setFileName( repFourierFile.str().c_str() );
						mreader.read( R );
						if ( R.size() > 0 ){
							cout << "Using pre-computed reciprocal space distances from file:\n\t" << repFourierFile.str().c_str() << "\n\tSize = " << R.shape() << endl;
						} else {
							cout << "Computing reciprocal space representations..." << endl;
							this->metric->getRepresentations( volumesFile, R, this->useRealRepresentation, this->binFactorForClassification );
							w.setFileName( repFourierFile.str().c_str() );
							if ( w.write( R ) == false ){
								cerr << "ERROR: Failed to write " << repFourierFile.str() << endl;
							}
						}

						if ( false ){
							stringstream distancesFile;
							if ( iteration < 10 ){
								distancesFile << this->fileHeader << "_iteration_00" << iteration << "_reciprocalSpaceDistances.matlab";
							} else if ( iteration < 100 ){
								distancesFile << this->fileHeader << "_iteration_0" << iteration << "_reciprocalSpaceDistances.matlab";
							} else {
								distancesFile << this->fileHeader << "_iteration_" << iteration << "_reciprocalSpaceDistances.matlab";
							}

							D.free();

							mreader.setFileName( distancesFile.str().c_str() );
							mreader.read( D );

							if ( D.size() > 0 ){
								cout << "Using pre-computed reciprocal space distances from file:\n\t" << distancesFile.str().c_str() << "\n\tSize = " << D.shape() << endl;
							} else {
								cout << "Computing reciprocal space distances..." << endl;
								this->metric->getDistances( repFourierFile, D );
								w.setFileName( distancesFile.str().c_str() );
								if ( w.write( D ) == false ){
									cerr << "ERROR: Failed to write " << distancesFile.str() << endl;
								}
							}

							// IMPUTE WITH K-NEAREST NEIGHBORS
							cout << "Imputing values in reciprocal space representation..." << endl;
							this->imputeInReciprocalSpace( R, D );

							// concatenate real and imaginary data
							Rimputed.resize( R.rows(), 2 * R.cols(), 1 );
							Rimputed( Range :: all(), Range( fromStart, R.ubound(secondDim) ), 0 ) = R( Range :: all(), Range :: all(), 0 );
							Rimputed( Range :: all(), Range( R.cols(), toEnd ), 0 ) = R( Range :: all(), Range :: all(), 1 );

							w.setFileName( imputedRepFile.str().c_str() );
							if ( w.write( Rimputed ) == false ){
								cerr << "Failed to write " << imputedRepFile.str() << endl;
							}
						}
					}
					R.resize( Rimputed.shape() );
					R = Rimputed;
				}
				w.setFileName( repFile.str().c_str() );
				// DOT NOT SAVE IN ORDER TO SAVE SPACE
				if ( w.write( R ) == false ){
					cerr << "Failed to write " << repFile.str() << endl;
				}
			}

			if ( R.size() > 0 ){
				cout << "Reducing dimensionality by SVD decomposition...\n" << R.shape() << endl;

				//w.setFileName("p.matlab");
				//w.write( R );

				//Array< Pixel, 3 > Rsmall( R.rows() / 2, R.cols(), R.depth() );
				//Rsmall = R( Range(fromStart,toEnd,2), Range::all(), Range::all() );
				this->doSVDDecomposition( R, 16 );

				//Array< Pixel, 3 > C( Rsmall( Range::all(), Range(1,toEnd), Range::all() ) );
				
				//w.write( Rsmall );

				//this->doSVDDecomposition( Rsmall, 17 );

				//w.write( Rsmall );

				w.setFileName( repSVDFile.str().c_str() );
				// DOT NOT SAVE IN ORDER TO SAVE SPACE
				if ( w.write( R ) == false ){
					cerr << "Failed to write " << repSVDFile.str() << endl;
				}
			}
		}

		//if ( R.size() > 0 ){
		//	cout << "Computing distances in lower dimensional space..." << endl;
		//	this->metric->getDistances( repSVDFile, D );
		//} else {
		//	this->metric->getDistances( repFourierFile, D );
		//}
		//w.setFileName( distancesFileSVD.str().c_str() );
		//if ( w.write( D ) == false ){
		//	cerr << "ERROR: Failed to write " << distancesFileSVD.str() << endl;
		//}

		this->computeInitialReferencesKmeans(R,iteration);

	} else {

		// transfer classification result to class attribute
		this->classes.resize( hierarchicalClasses, this->input.size(), 1 );
		this->classes = 1;

		// reset current references
		this->references.clear();

		// build class averages
		for ( int i = 0; i < this->classes.rows(); i++ ){
			nbfWedgedAverageImage3D< Pixel > average( this->input );

			// assign alignments
			Array< Pixel, 3 > referenceAlignments( this->input.size(), 17, 1 );
			referenceAlignments = 0;

			// set new weights
			referenceAlignments( Range::all(), 0, 0 ) = this->classes( i, Range::all(), 0 );

			for ( int j = 0; j < this->classes.cols(); j++ ){
				if ( this->classes( i, j, 0 ) > 0 ){
					// set alignment to best match
					vtkTransform * t = vtkTransform::New();
					this->input[j].getTransform(t);
					double matrix[16];
					vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
					t->Delete();
					for ( int e = 0; e < 16; e++ ){
						referenceAlignments( j, 1 + e, 0 ) = matrix[e];
					}
				}
			}

			average.setAlignments( referenceAlignments );

			this->references.push_back( average );
		}
	}

	//// initialize distance alignments properly
	////this->alignments.resize( D.rows(), D.cols(), 19, 1 );
	//// only use distance and overlaps
	//this->alignments.resize( D.rows(), D.cols(), 2, 1 );
	//
	//// set distance matrix in first component
	//this->alignments( Range::all(), Range::all(), 0, 0 ) = D;
	////this->alignments( Range::all(), Range::all(), 0, 0 ) = D( Range::all(), Range::all(), 0, 0 );

	//// set overlaps and scalings to unit
	//this->alignments( Range::all(), Range::all(), Range(1,toEnd), 0 ) = 1;

	//stringstream treeFile;
	//treeFile << this->fileHeader << "_iteration_" << iteration << "_SVD_distancesTree.matlab";

	//Array< Pixel, 2 > T;
	//
	//mreader.setFileName( treeFile.str().c_str() );
	//mreader.read( T );

	//bool existingTree = false;

	//if ( T.size() > 0 ){
	//	cout << "Using pre-computed distance tree from file:\n\t" << treeFile.str().c_str() << "\n\tSize = " << T.shape() << endl;
	//	existingTree = true;
	//}

	//cout << "Computing hierarchical clustering tree..." << endl;
	//this->computeInitialReferences(T,iteration);
	//
	//if ( existingTree == false ){
	//	w.setFileName( treeFile.str().c_str() );
	//	if ( w.write( T ) == false ){
	//		cerr << "ERROR: Failed to write " << treeFile.str() << endl;
	//	}
	//}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: computeNewReferencesSpectral( int iteration )
{
	// apply alignment of best match and eliminate outliers
	this->asignBestAlignmentAndClean();

	// save current alignments into file
	stringstream volumesFile;
	if ( iteration < 10 ){
		volumesFile << this->fileHeader << "_iteration_00" << iteration << "_volumes_tmp.txt";
	} else if ( iteration < 100 ){
		volumesFile << this->fileHeader << "_iteration_0" << iteration << "_volumes_tmp.txt";
	} else {
		volumesFile << this->fileHeader << "_iteration_" << iteration << "_volumes_tmp.txt";
	}
	nbfWedgedSubImage3D< Pixel > :: write( volumesFile.str(), this->input );

	nbfMatlabReader mreader;
	nbfMatlabWriter w;

	// where to store svd coordinates
	stringstream repSVDFile;
	if ( iteration < 10 ){
		repSVDFile << this->fileHeader << "_iteration_00" << iteration << "_SVD_representation.matlab";
	} else if ( iteration < 100 ){
		repSVDFile << this->fileHeader << "_iteration_0" << iteration << "_SVD_representation.matlab";
	} else {
		repSVDFile << this->fileHeader << "_iteration_" << iteration << "_SVD_representation.matlab";
	}
	mreader.setFileName( repSVDFile.str().c_str() );
	Array< Pixel, 2 > D;
	mreader.read( D );

	if ( D.size() > 0 ){
		cout << "Using pre-computed SVD coordinates from file:\n\t" << repSVDFile.str().c_str() << "\n\tSize = " << D.shape() << endl;
	} else {

		// where to store distance matrix
		stringstream distancesFile;
		if ( iteration < 10 ){
			distancesFile << this->fileHeader << "_iteration_00" << iteration << "_reciprocalSpaceDistances.matlab";
		} else if ( iteration < 100 ){
			distancesFile << this->fileHeader << "_iteration_0" << iteration << "_reciprocalSpaceDistances.matlab";
		} else {
			distancesFile << this->fileHeader << "_iteration_" << iteration << "_reciprocalSpaceDistances.matlab";
		}

		D.free();
		mreader.setFileName( distancesFile.str().c_str() );
		mreader.read( D );

		if ( D.size() > 0 ){
			cout << "Using pre-computed reciprocal space distances from file:\n\t" << distancesFile.str().c_str() << "\n\tSize = " << D.shape() << endl;
		} else {
			////////////////////////////
			// STORE ONLY TEMPORARILY //
			////////////////////////////
			stringstream reducedRepreresentationFile;
			reducedRepreresentationFile << "mpi_999.tmp";
			//if ( iteration < 10 ){
			//	reducedRepreresentationFile << this->fileHeader << "_iteration_00" << iteration << "_reduced_representation.matlab";
			//} else if ( iteration < 100 ){
			//	reducedRepreresentationFile << this->fileHeader << "_iteration_0" << iteration << "_reduced_representation.matlab";
			//} else {
			//	reducedRepreresentationFile << this->fileHeader << "_iteration_" << iteration << "_reduced_representation.matlab";
			//}

			Array< Pixel, 3 > R;
			mreader.setFileName( reducedRepreresentationFile.str().c_str() );
			mreader.read( R );

			if ( R.size() > 0 ){
				cout << "Using pre-computed reduced representations from file:\n\t" << reducedRepreresentationFile.str().c_str() << "\n\tSize = " << R.shape() << endl;
			} else {
				cout << "Computing reduced representations..." << endl;
				this->metric->getRepresentations( volumesFile, R, 0, this->binFactorForClassification, this->symmetryFactor );

				// OVERRIDE MISSING WEDGE INFORMATION
				R( Range :: all(), Range :: all(), 2 ) = 1;

				/////// DO NORMALIZATION TO MIMIC MODULATION METRIC

				firstIndex i;
				secondIndex j;

				// left matrix metric
				Array< Pixel, 2 > Rreal(  R( Range :: all(), Range :: all(), 0 ) );
				Array< Pixel, 2 > Rimag(  R( Range :: all(), Range :: all(), 1 ) );
				Array< Pixel, 2 > Rwedge( R( Range :: all(), Range :: all(), 2 ) );
				Array< Pixel, 2 > RmoduloSqr( R.rows(), R.cols() );
				RmoduloSqr = ( pow2( Rreal ) + pow2( Rimag ) ) * Rwedge; 
				Array< Pixel, 1 > M( R.rows() );
				M = where( sum( Rwedge(i,j), j ) > 0, sqrt( sum( RmoduloSqr(i,j), j ) ) / sum( Rwedge(i,j), j ), sqrt( sum( RmoduloSqr(i,j), j ) ) );
				M = where( sqrt(M) > 0, sqrt(M), 1 );
				// M = 1.0 / sqrt( sum( pow2( A(i,j) ), j ) );

				// right matrix metric
				Array< Pixel, 1 > N( R.cols() );
				N = where( sum( Rwedge(j,i), j ) > 0, sqrt( sum( RmoduloSqr(j,i), j ) ) / sum( Rwedge(j,i), j ), sqrt( sum( RmoduloSqr(j,i), j ) ) );
				N = where( sqrt(N) > 0, sqrt(N), 1 );
				// N = 1.0 / sqrt( sum( pow2( A(j,i) ), j ) );

				for ( int i = 0; i < R.cols(); i++ ){
					Rreal( Range :: all(), i ) /= M;
					Rimag( Range :: all(), i ) /= M;
				}

				for ( int i = 0; i < R.rows(); i++ ){
					Rreal( i, Range :: all() ) /= N;
					Rimag( i, Range :: all() ) /= N;
				}

				w.setFileName( reducedRepreresentationFile.str().c_str() );
				if ( w.write( R ) == false ){
					cerr << "ERROR: Failed to write " << reducedRepreresentationFile.str() << endl;
				}
			}

			cout << "Computing reciprocal space distances..." << endl;
			this->metric->getDistances( reducedRepreresentationFile, D );
			w.setFileName( distancesFile.str().c_str() );
			//if ( w.write( D ) == false ){
			//	cerr << "ERROR: Failed to write " << distancesFile.str() << endl;
			//}
		}

		cout << "Reducing dimensionality by SVD decomposition...\n" << D.shape() << endl;
		Array< Pixel, 3 > R( D.rows(), D.cols(), 1 );
		R( Range :: all(), Range :: all(), 0 ) = D;
		this->doSpectralSVDDecomposition( R, this->hierarchicalClasses );
		
		// add dummy coordinate because the K-means ignores the first coordinate
		D.resize( R.rows(), R.cols() + 1 );
		D( Range :: all(), 0 ) = 0;
		D( Range :: all(), Range(1,toEnd) ) = R( Range :: all(), Range :: all(), 0 );

		//D.resize( R.rows(), R.cols() );
		//D = R( Range :: all(), Range :: all(), 0 );

		w.setFileName( repSVDFile.str().c_str() );
		//if ( w.write( D ) == false ){
		//	cerr << "Failed to write " << repSVDFile.str() << endl;
		//}
	}
	Array< Pixel, 3 > R( D.rows(), D.cols(), 1 );
	R( Range :: all(), Range :: all(), 0 ) = D;
	this->computeInitialReferencesKmeans( R, iteration );
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: buildInitialReference( nbfWedgedAverageImage3D< Pixel > & reference, int referenceNumber )
{
	vector< int > classIndexes;
	for ( int j = 0; j < this->classes.cols(); j++ ){
		if ( this->classes( referenceNumber, j ) > 0 ){
			classIndexes.push_back(j);
		}
	}

	// build two-particle initial averages (self + closest one)
	vector< nbfWedgedAverageImage3D< Pixel > > current;
	nbfWedgedAverageImage3D< Pixel > currentAverage( this->input );
	
	vector< nbfWedgedAverageImage3D< Pixel > > list1;
	vector< nbfWedgedSubImage3D< Pixel > > list2;

	// store pair-wise distance/alignment matrix of all volumes in class
	Array< Pixel, 3 > allConfigurations( classIndexes.size(), classIndexes.size(), 19 );

	// seed with each volume
	for ( int i = 0; i < classIndexes.size(); i++ ){
	
		// keep track of volumes in average
		vector< int > averagedVolumes;
		
		Array< Pixel, 2 > currentAlignments( this->input.size(), 17 );
		currentAlignments = 0;

		// add current particle with identity transform
		averagedVolumes.push_back( i );
		currentAlignments( classIndexes[i], 0 ) = 1;
		currentAlignments( classIndexes[i], Range(1,toEnd) ) = this->alignments( classIndexes[i], classIndexes[i], Range(3,toEnd) );

		// find closest particle in class
		int minIndex; 
		Pixel minDistance = numeric_limits< Pixel > :: max();
		for ( int j = 0; j < classIndexes.size(); j++ ){
			if ( ( i != j ) && ( this->alignments( classIndexes[i], classIndexes[j], 0 ) < minDistance ) ){
				minIndex = j; minDistance = this->alignments( classIndexes[i], classIndexes[j], 0 );
			}
		}

		// add closest particle aligned to first particle
		averagedVolumes.push_back( minIndex );
		currentAlignments( classIndexes[minIndex], 0 ) = 1;
		currentAlignments( classIndexes[minIndex], Range(1,toEnd) ) = this->alignments( classIndexes[i], classIndexes[minIndex], Range(3,toEnd) );
		
		// store transformations
		allConfigurations( i, i, Range::all() ) = this->alignments( classIndexes[i], classIndexes[i], Range::all() );
		allConfigurations( i, minIndex, Range::all() ) = this->alignments( classIndexes[i], classIndexes[minIndex], Range::all() );

		// set average alignments
		currentAverage.setAlignments( currentAlignments );

		// store indexes of computed distances
		vector< int > computedDistances;

		while ( averagedVolumes.size() < classIndexes.size() ){

			list1.clear(); list2.clear(); computedDistances.clear();

			// compute distances to remaining volumes
			for ( int j = 0; j < classIndexes.size(); j++ ){
				// see if not already in average
				vector< int > :: iterator iter = find( averagedVolumes.begin(), averagedVolumes.end(), j );
				if ( iter == averagedVolumes.end() ){
					this->input[ classIndexes[j] ].setTransform( (vtkTransform*)NULL );
					list1.push_back( currentAverage );
					list2.push_back( this->input[ classIndexes[j] ] );
					computedDistances.push_back(j);
				}
			}

			Array< Pixel, 2 > D( list1.size(), 19 );
			D = -1;
			
			vector< nbfWedgedImage3D< Pixel > * > plist1, plist2;
			vector< TinyVector< int, 2 > > positions;
			for ( int j = 0; j < list1.size(); j++ ){
				plist1.push_back( &(list1[j]) );
				plist2.push_back( &(list2[j]) );
				TinyVector< int, 2 > pos(j,j);
				positions.push_back( pos );
			}
			this->metric->getDistances( plist1, plist2, positions, D );
			
			// find closest to current average
			minDistance = numeric_limits< Pixel > :: max();
			for ( int j = 0; j < computedDistances.size(); j++ ){
				if ( D( j, 0 ) < minDistance ){
					minIndex = j; minDistance = D( j, 0 );
				}
			}

			// add to average
			averagedVolumes.push_back( computedDistances[ minIndex ] );
			currentAlignments( classIndexes[ computedDistances[minIndex] ], 0 ) = 1;
			currentAlignments( classIndexes[ computedDistances[minIndex] ], Range(1,toEnd) ) = D( minIndex, Range(3,toEnd) );

			// set new alignments
			currentAverage.setAlignments( currentAlignments );

			// store transformation
			allConfigurations( i, computedDistances[minIndex], Range::all() ) = D( minIndex, Range::all() );
		}
	}

	//cout << allConfigurations( Range::all(), Range::all(), 0 ) << endl;

	// evaluate energy for each configuration
	Array< Pixel, 4 > energies( classIndexes.size(), classIndexes.size(), classIndexes.size(), 19 );
	energies = 0;

	for ( int i = 0; i < allConfigurations.rows(); i++ ){
		// compute sum of all pairwise distances
		for ( int j = 0; j < allConfigurations.cols(); j++ ){
			for ( int k = j + 1; k < allConfigurations.cols(); k++ ){
				double matrix[16];
				for ( int m = 0; m < 16; m++ ){
					matrix[m] = allConfigurations( i, j, 3 + m );
				}
				vtkMatrix4x4 * mat = vtkMatrix4x4::New();
				mat->DeepCopy( matrix );

				this->input[ classIndexes[j] ].setTransform( mat );

				nbfFourierImageMetric< Pixel, 3 > fMetric( this->metric->imageFilter, this->metric->fourierFilter );

				fMetric.setInput1( reinterpret_cast< nbfWedgedImage3D< Pixel > * >( & ( this->input[ classIndexes[j] ] ) ) );

				vtkTransform * trans2 = vtkTransform::New();
				for ( int m = 0; m < 16; m++ ){
					matrix[m] = allConfigurations( i, k, 3 + m );
				}
				mat->DeepCopy( matrix );
				trans2->Concatenate( mat );
				this->input[ classIndexes[k] ].setTransform( (vtkTransform*)NULL );
				fMetric.setInput2( reinterpret_cast< nbfWedgedImage3D< Pixel > * >( & ( this->input[ classIndexes[k] ] ) ) );

				fMetric.executeFourierNewHalf( trans2 );
				trans2->Delete();

				mat->Delete();

				energies(i,j,k,0) = fMetric.getCorrelationPeak();
				energies(i,j,k,2) = fMetric.getWedgeOverlap();
			}
		}
		Array< Pixel, 3 > dummy( energies( i, Range::all(), Range::all(), Range::all() ) );
		this->metric->makeDistanceMatrix( dummy );
	}
	
	// detect best configuration

	nbfHierarchicalClustering< Pixel > hierarch;
	hierarch.setMaxDistance( this->hierarchicalClasses );
	hierarch.setMinOverlap( 0.0 );
	hierarch.setMinElementNumber( 2 );

	int topConfigurationSize = 0;
	Pixel topConfigurationEnergy = numeric_limits< Pixel > :: max();
	Array< Pixel, 2 > bestAlignments;

	vector< nbfWedgedSubImage3D< Pixel > > subInput;
	for ( int i = 0; i < classIndexes.size(); i++ ){
		subInput.push_back( this->input[ classIndexes[ i ] ] );
	}

	for ( int i = 0; i < classIndexes.size(); i++ ){
		hierarch.setInput( subInput );
		hierarch.setMetric( this->metric, energies( i, Range::all(), Range::all(), Range::all() ) );
		//cout << energies( i, Range::all(), Range::all(), 0 ) << endl;
		Array< Pixel, 2 > subClasses;
		hierarch.execute( subClasses );
		for ( int p = 0; p < subClasses.rows(); p++ ){
			// count number of elements in class
			int classSize = sum( subClasses( p, Range::all() ) );
			if ( classSize >= topConfigurationSize ){
				Pixel energy = 0;
				for ( int di = 0; di < subClasses.cols(); di++ ){
					for ( int dj = di + 1; dj < subClasses.cols(); dj++ ){
						if ( ( subClasses(p,di) == 1 ) && ( subClasses(p,dj) == 1 ) ){
							energy += energies( i, di, dj, 0 );
						}
					}
				}
				if ( ( energy < topConfigurationEnergy ) || ( classSize > topConfigurationSize ) ){
					topConfigurationEnergy = energy;
					topConfigurationSize = classSize;
					bestAlignments.resize( classSize, 19 );
					bestAlignments = allConfigurations( i, Range::all(), Range::all() );
					bestAlignments( Range::all(), 0 ) = subClasses( p, Range::all() );
				}
			}
		}
	}

	// assign best set of alignments
	Array< Pixel, 2 > referenceAlignments( this->input.size(), 17 );
	
	// set weights
	referenceAlignments = 0;
	for ( int i = 0; i < classIndexes.size(); i++ ){
		cout << "adding volume " << classIndexes[i] << " to average." << endl;
		referenceAlignments( classIndexes[i], 0 ) = 1.0; // bestAlignments( i, 0 );
		referenceAlignments( classIndexes[i], Range(1,toEnd) ) = bestAlignments( i, Range(3,toEnd) );
		//referenceAlignments( classIndexes[i], Range(1,toEnd) ) = this->alignments( classIndexes[0], classIndexes[i], Range(3,toEnd) );
	}

	//cout << referenceAlignments( Range::all(), 0 ) << endl;
	// set new alignments
	reference.setAlignments( referenceAlignments );

	// assign alignment to references
	this->alignmentToReferences( referenceNumber, Range::all(), 0 ) = referenceAlignments( Range::all(), 0 );
	this->alignmentToReferences( referenceNumber, Range::all(), Range(3,toEnd) ) = referenceAlignments( Range::all(), Range(1,toEnd) );

	//vtkImageData * data = vtkImageData::New();
	//reference.getImage( data );
	//vtkStructuredPointsWriter * vw1 = vtkStructuredPointsWriter::New();
	//vw1->SetInput( data );
	//vw1->SetFileName("averageRefined.vtk");
	//vw1->Write();
	//vw1->Delete();
	//data->Delete();
	//exit(0);
}


template< class Pixel >
void nbfLoopClustering< Pixel > :: buildInitialReferences()
{
	// store pair-wise distance/alignment matrix of all volumes in class
	vector< Array< Pixel, 2 > > allCurrentAlignments;
	
	// store volumes currently contributing to each average
	vector< vector< int > > allAveragedVolumes;

	// store indexes of already computed distances
	vector< vector< int > > allComputedDistances;

	// use a sparse tree class representation for convenience
	vector< vector< int > > allClassIndexes;

	// store particle within class which is a best candidate for seeding
	vector< int > seeds;

	// build class trees
	for ( int referenceNumber = 0; referenceNumber < this->classes.rows(); referenceNumber++ ){
		vector< int > classIndexes;
		for ( int j = 0; j < this->classes.cols(); j++ ){
			if ( this->classes( referenceNumber, j ) > 0 ){
				classIndexes.push_back(j);
			}
		}
		allClassIndexes.push_back( classIndexes );
		
		// build configuration volumes
		Array< Pixel, 2 > configurations( this->classes.cols(), 19 );

		// find element in class with the most closest neighbors
		int best = 0;
		Pixel bestDistance = numeric_limits< Pixel > :: max();
		for ( int i = 0; i < allClassIndexes[ referenceNumber ].size(); i++ ){
			Pixel currentD = 0;
			for ( int j = 0; j < allClassIndexes[ referenceNumber ].size(); j++ ){
				currentD += this->alignments( allClassIndexes[ referenceNumber ][i], allClassIndexes[ referenceNumber ][j], 0 );
			}
			if ( currentD < bestDistance ){
				bestDistance = currentD;
				best = i;
			}
		}
		seeds.push_back( best );
	}

	vector< nbfWedgedAverageImage3D< Pixel > > list1;
	vector< nbfWedgedSubImage3D< Pixel > > list2;

	vector< nbfWedgedAverageImage3D< Pixel > > current;
	nbfWedgedAverageImage3D< Pixel > currentAverage( this->input );

	//// compute distances of volumes within each class to initial 1-particle averages
	for ( int referenceNumber = 0; referenceNumber < this->classes.rows(); referenceNumber++ ){

		// keep track of volumes and distances for current reference
		vector< int > averagedVolumes;
		vector< int > computedDistances;

		int currentIndex = seeds[ referenceNumber ];

		Array< Pixel, 2 > currentAlignments( this->input.size(), 17 );
		currentAlignments = 0;

		// add current particle with identity transform
		int globalIndex = allClassIndexes[ referenceNumber ][ currentIndex ];
		averagedVolumes.push_back( globalIndex );
		currentAlignments( globalIndex, 0 ) = 1;
		currentAlignments( globalIndex, Range(1,toEnd) ) = this->alignments( globalIndex, globalIndex, Range(3,toEnd) );

		// set average alignments
		currentAverage.setAlignments( currentAlignments );

		// compute distances to remaining volumes
		for ( int j = 0; j < allClassIndexes[ referenceNumber ].size(); j++ ){
			if ( currentAlignments( allClassIndexes[ referenceNumber ][j], 0 ) == 0 ){ // if not already in average
				this->input[ allClassIndexes[ referenceNumber ][j] ].setTransform( (vtkTransform*)NULL ); // reset transformation
				list1.push_back( currentAverage );
				list2.push_back( this->input[ allClassIndexes[ referenceNumber ][ j ] ] );
				computedDistances.push_back( j );
			}
		}

		allComputedDistances.push_back( computedDistances );
		allCurrentAlignments.push_back( currentAlignments );
		allAveragedVolumes.push_back( averagedVolumes );
	}

	// add closest volume to current reference and re-compute distances to remaining volumes. repeat.
	do {
		// compute all distances with MPI
		Array< Pixel, 2 > D( list1.size(), 19 );
		D = -1;

		vector< nbfWedgedImage3D< Pixel > * > plist1, plist2;
		vector< TinyVector< int, 2 > > positions;
		for ( int j = 0; j < list1.size(); j++ ){
			plist1.push_back( &(list1[j]) );
			plist2.push_back( &(list2[j]) );
			TinyVector< int, 2 > pos(j,j);
			positions.push_back( pos );
		}

		// compute cheap distances (no refinement) in favour of speed
		this->metric->getDistances( plist1, plist2, positions, D, 2 );

		list1.clear(); list2.clear();

		// interpret results

		int count = 0;

		for ( int referenceNumber = 0; referenceNumber < this->classes.rows(); referenceNumber++ ){

			int currentIndex = seeds[ referenceNumber ];
			
			if ( allComputedDistances[ referenceNumber ].size() > 0 ){
				// find closest to current average
				int minLocalIndex, minGlobalIndex, minDistanceIndex;
				Pixel minDistance = numeric_limits< Pixel > :: max();
				for ( int j = 0; j < allComputedDistances[ referenceNumber ].size(); j++ ){
					if ( D( count, 0 ) < minDistance ){
						minLocalIndex = allComputedDistances[ referenceNumber ][j];
						minGlobalIndex = allClassIndexes[ referenceNumber ][ minLocalIndex ];
						minDistance = D( count, 0 );
						minDistanceIndex = count;
					}
					count++;
				}

				cout << "Add volume " << minGlobalIndex << " to cummulative average " << referenceNumber << ", D = " << minDistance << endl;

				// add to average
				allAveragedVolumes[ referenceNumber ].push_back( minGlobalIndex );
				allCurrentAlignments[ referenceNumber ]( minGlobalIndex, 0 ) = 1;
				allCurrentAlignments[ referenceNumber ]( minGlobalIndex, Range(1,toEnd) ) = D( minDistanceIndex, Range(3,toEnd) );

				// set new alignments
				currentAverage.setAlignments( allCurrentAlignments[ referenceNumber ] );

				allComputedDistances[ referenceNumber ].clear();

				// compute distances to remaining volumes
				for ( int j = 0; j < allClassIndexes[ referenceNumber ].size(); j++ ){
					// see if not already in average
					vector< int > :: iterator iter = find( allAveragedVolumes[ referenceNumber ].begin(), allAveragedVolumes[ referenceNumber ].end(), allClassIndexes[ referenceNumber ][j] );
					if ( iter == allAveragedVolumes[ referenceNumber ].end() ){
						this->input[ allClassIndexes[ referenceNumber ][j] ].setTransform( (vtkTransform*)NULL );
						list1.push_back( currentAverage );
						list2.push_back( this->input[ allClassIndexes[ referenceNumber ][j] ] );
						allComputedDistances[ referenceNumber ].push_back( j );
					}
				}
			}
		}
	} while ( list1.size() > 0 );

	// cout << "All configurations " << allConfigurations( Range::all(), Range::all(), 0 ) << endl;

	// given the above transformations, compute all pairwise distances within each class

	//vector< nbfWedgedSubImage3D< Pixel > > list3;
	//vector< nbfWedgedSubImage3D< Pixel > > list4;
	//vector< TinyVector< int, 2 > > positions;

	// for all references
	//for ( int referenceNumber = 0; referenceNumber < this->classes.rows(); referenceNumber++ ){

	//	// for each volume in this reference
	//	for ( int j = 0; j < allClassIndexes[ referenceNumber ].size(); j++ ){
	//		double matrix[16];
	//		for ( int m = 0; m < 16; m++ ){
	//			matrix[m] = allCurrentAlignments[ referenceNumber ]( allClassIndexes[ referenceNumber ][ j ], 1 + m );
	//		}
	//		vtkMatrix4x4 * mat = vtkMatrix4x4::New();
	//		mat->DeepCopy( matrix );

	//		this->input[ allClassIndexes[ referenceNumber ][j] ].setTransform( mat );

	//		list3.push_back( this->input[ allClassIndexes[referenceNumber][j] ] );

	//		for ( int k = j + 1; k < allClassIndexes[ referenceNumber ].size(); k++ ){

	//			for ( int m = 0; m < 16; m++ ){
	//				matrix[m] = allCurrentAlignments[ referenceNumber ]( allClassIndexes[ referenceNumber ][ k ], 1 + m );
	//			}
	//			mat->DeepCopy( matrix );

	//			this->input[ allClassIndexes[referenceNumber][k] ].setTransform( mat );

	//			list4.push_back( this->input[ allClassIndexes[referenceNumber][k] ] );
	//			TinyVector< int, 2 > pos( list3.size() - 1, list4.size() - 1 );
	//			positions.push_back( pos );
	//		}
	//		mat->Delete();
	//	}
	//}

	//vector< nbfWedgedImage3D< Pixel > * > plist1, plist2;
	//for ( int j = 0; j < list3.size(); j++ ){
	//	plist1.push_back( &(list3[j]) );
	//}
	//for ( int j = 0; j < list4.size(); j++ ){
	//	plist2.push_back( &(list4[j]) );
	//}
	//Array< Pixel, 2 > allDistances( positions.size(), 19 );
	//allDistances = -1;
	//this->metric->getDistances( plist1, plist2, positions, allDistances, 2 );

	//int count = 0;

	for ( int referenceNumber = 0; referenceNumber < this->classes.rows(); referenceNumber++ ){

		//// evaluate energy for each configuration
		//Array< Pixel, 3 > energies( allClassIndexes[ referenceNumber ].size(), allClassIndexes[ referenceNumber ].size(), 19 );
		//energies = 0;

		//// compute sum of all pairwise distances
		//for ( int j = 0; j < allClassIndexes[ referenceNumber ].size(); j++ ){
		//	for ( int k = j + 1; k < allClassIndexes[ referenceNumber ].size(); k++ ){
		//		energies(j,k,0) = allDistances(count,0);
		//		energies(j,k,2) = allDistances(count,2);
		//		count++;
		//	}
		//}
		//this->metric->makeDistanceMatrix( energies );
	
		//// detect best configuration

		//nbfHierarchicalClustering< Pixel > hierarch;
		//hierarch.setMaxDistance( this->hierarchicalCutoff );
		//hierarch.setMinOverlap( 0.0 );
		//hierarch.setMinElementNumber( 0 );

		//vector< nbfWedgedSubImage3D< Pixel > > subInput;
		//for ( int i = 0; i < allClassIndexes[ referenceNumber ].size(); i++ ){
		//	subInput.push_back( this->input[ allClassIndexes[ referenceNumber ][ i ] ] );
		//}

		//hierarch.setInput( subInput );
		//hierarch.setMetric( this->metric, energies );

		//Array< Pixel, 2 > subClasses;
		//hierarch.execute( subClasses );

		//// find most significant class
		//secondIndex j;
		//TinyVector< int, 1 > bestClass = maxIndex( sum( subClasses, j ) );

		// assign best set of alignments
		Array< Pixel, 2 > referenceAlignments( this->input.size(), 17 );

		// set new weights
		//allCurrentAlignments[ referenceNumber ]( Range::all(), 0 ) = 0;
		//for ( int i = 0; i < subClasses.cols(); i++ ){
		//	allCurrentAlignments[ referenceNumber ]( allClassIndexes[ referenceNumber ][i], 0 ) = subClasses( bestClass[0], i );
		//}

		// allCurrentAlignments[ referenceNumber ]( Range::all(), 0 ) = this->classes( referenceNumber, Range::all() );

		// set final alignments
		nbfWedgedAverageImage3D< Pixel > cAverage( this->input );
		cAverage.setAlignments( allCurrentAlignments[ referenceNumber ] );

		for ( int i = 0; i < allClassIndexes[ referenceNumber ].size(); i++ ){
			cout << "adding volume " << allClassIndexes[ referenceNumber ][i] << " to average with weight "  << allCurrentAlignments[ referenceNumber ]( allClassIndexes[ referenceNumber ][i], 0 ) << endl;
		}

		// assign alignment to references
		this->alignmentToReferences( referenceNumber, Range::all(), 0 ) = allCurrentAlignments[ referenceNumber ]( Range::all(), 0 );
		this->alignmentToReferences( referenceNumber, Range::all(), 1 ) = 1;
		this->alignmentToReferences( referenceNumber, Range::all(), Range(3,toEnd) ) = allCurrentAlignments[ referenceNumber ]( Range::all(), Range(1,toEnd) );

		this->references.push_back( cAverage );
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: bundleAlignment( vector< nbfWedgedSubImage3D< Pixel > > & volumes )
{
	// align all volumes to first in list

	for ( int iter = 0; iter < this->refinementIterations; iter++ ){
		// center all volumes to global average
		nbfWedgedAverageImage3D< Pixel > globalAverage( volumes );
		Array< Pixel, 3 > globalAlignments( volumes.size(), 17, 1 );
		globalAlignments = 1;
		double matrix[16];
		for ( int i = 0; i < volumes.size(); i++ ){
			vtkTransform * t = vtkTransform::New();
			volumes[i].getTransform(t);
			vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
			t->Delete();
			for ( int e = 0; e < 16; e++ ){
				globalAlignments( i, 1 + e, 0 ) = matrix[e];
			}
		}
		globalAverage.setAlignments( globalAlignments );
		this->metric->getImage( globalAverage );

		// reset original alignments
		for ( int i = 0; i < volumes.size(); i++ ){
			volumes[i].setTransform( (vtkTransform*)NULL );
		}

		vector< nbfWedgedAverageImage3D< Pixel > > reference;
		reference.push_back( globalAverage );
		Array< Pixel, 4 > ALS( 1, volumes.size(), 19, 1 );
		ALS = -1;
		
		// align with most accurate metric
		this->metric->getDistances( reference, volumes, ALS, 0 );

		// apply alignments
		for ( int i = 0; i < ALS.cols(); i++ ){
			for ( int j = 0; j < 16; j++ ){
				matrix[j] = ALS( 0, i, 3 + j, 0 );
			}
			vtkMatrix4x4 * matriz = vtkMatrix4x4::New();
			matriz->DeepCopy(matrix);
			volumes[i].setTransform( matriz );
			matriz->Delete();
		}
	}

	// rectify all transforms to first reference
	if ( volumes.size() > 0 ){
		for ( int i = 0; i < volumes.size(); i++ ){
			vtkTransform * t = vtkTransform :: New();
			volumes[0].getTransform(t);
			t->Inverse();
			vtkTransform * ct = vtkTransform::New();
			volumes[i].getTransform(ct);
			ct->Concatenate(t);
			volumes[i].setTransform(ct);
			ct->Delete();
			t->Delete();
		}
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: bundleAlignment( vector< nbfWedgedAverageImage3D< Pixel > > & akareferences, int alignment_mode )
{
	// save starting class memberships
	Array< Pixel, 2 > startingWeights( akareferences.size(), this->input.size() );
	startingWeights = 0;
	for ( int i = 0; i < akareferences.size(); i++ ){
		startingWeights( i, Range::all() ) = akareferences[i].weights( Range::all(), 0 );
	}

	for ( int iterations = 0; iterations < this->refinementIterations; iterations++ ){

		cout << "Bundle alignment iteration " << iterations + 1 << " of " << this->refinementIterations << endl;

		vector< nbfWedgedAverageImage3D< Pixel > > list1;
		vector< nbfWedgedSubImage3D< Pixel > > list2;
		
		vector< TinyVector< int, 2 > > indexesVector;
		vector< nbfWedgedImage3D< Pixel > * > list1References;
		vector< nbfWedgedImage3D< Pixel > * > list2References;

		for ( int i = 0; i < akareferences.size(); i++ ){
			//// clear object to reduce its size for MPI
			//akareferences[i].initialize();
			this->metric->getImage( akareferences[i] );
			list1.push_back( akareferences[i] );

			//// align all volumes to current reference
			//nbfWedgedAverageImage3D< Pixel > globalAverage( this->input );
			//globalAverage = akareferences[i];
			//globalAverage.initialize();
			//this->metric->getImage( globalAverage );
			//list1.push_back( globalAverage );

			for ( int j = 0; j < akareferences[i].weights.rows(); j++ ){
				if ( fabs( startingWeights( i, j ) ) > 0 ){
					list2.push_back( akareferences[i].getVolumesRO()[j] );
				}
			}
		}

		int counter = 0;
		for ( int i = 0; i < akareferences.size(); i++ ){
			list1References.push_back( & list1[i] );
			for ( int j = 0; j < akareferences[i].weights.rows(); j++ ){
				if ( fabs( startingWeights( i, j ) ) > 0 ){
					indexesVector.push_back( TinyVector< int, 2 >( i, counter ) );
					list2References.push_back( & list2[ counter ] );
					counter++;
				}
			}
		}

		//for ( int i = 0; i < list1.size(); i++ ){
		//	list1References.push_back( & list1[i] );
		//	for ( int j = 0; j < list2.size(); j++ ){
		//		indexesVector.push_back( TinyVector< int, 2 >(i,j) );
		//		list2References.push_back( & list2[j] );
		//	}
		//}

		// don't compute overlaps and scaling
		Array< Pixel, 3 > D;
		D.resize( indexesVector.size(), 17, this->metric->getNumberOfCandidates() );
		D = -1;

		// apply shifts only (=2)
		// compute rotation + shifts refinement (=1)
		//if ( ( iterations % 2 ) == 0 ){
		//	// adjust shifts only for even iteration number
		//	this->metric->getDistances( list1References, list2References, indexesVector, D, 2 );
		//} else {
		//	// adjust shifts and refine translation locally for odd iteration number
			this->metric->getDistances( list1References, list2References, indexesVector, D, alignment_mode );
		//}

		cerr << "Bundle alignment distances: " << sum(D( Range::all(), 0, 0 )) << endl;

		// apply alignments
		int count = 0;
		for ( int i = 0; i < akareferences.size(); i++ ){
			
			Array< Pixel, 3 > newAlignments( akareferences[i].getVolumes().size(), 17, 1 );
			newAlignments = 0;
			newAlignments( Range::all(), 0, 0 ) = akareferences[i].weights( Range::all(), 0 );

			// store distances of individual volumes to global average
			Array< Pixel, 1 > distancesToGlobalAverageArray( akareferences[i].getVolumes().size() );
			distancesToGlobalAverageArray = numeric_limits< Pixel > :: max();
			vector< Pixel > distancesToGlobalAverage;

			for ( int j = 0; j < akareferences[i].getVolumes().size(); j++ ){
				if ( fabs( startingWeights( i, j ) ) > 0 ){
					newAlignments( j, Range(1,toEnd), 0 ) = D( count, Range(1,toEnd), 0 );

					// retrieve previous alignment
					double matrix[16];
					for ( int k = 0; k < 16; k++ ){
						matrix[k] = akareferences[i].multipleAlignments(j,k,0);
					}
					vtkMatrix4x4 * mat1 = vtkMatrix4x4 :: New();
					mat1->DeepCopy( matrix );

					// retrieve correction
					for ( int k = 0; k < 16; k++ ){
						matrix[k] = D( count, k + 1, 0 );
					}
					vtkMatrix4x4 * mat2 = vtkMatrix4x4 :: New();
					mat2->DeepCopy( matrix );

					// store distance to compute weights later
					distancesToGlobalAverageArray( j ) = D( count, 0, 0 );
					distancesToGlobalAverage.push_back( D( count, 0, 0 ) );

					vtkMatrix4x4 * mat3 = vtkMatrix4x4 :: New();
					vtkMatrix4x4 ::Multiply4x4( mat1, mat2, mat3 );

					vtkMatrix4x4 :: DeepCopy( matrix, mat3 );
					for ( int k = 0; k < 16; k++ ){
						newAlignments( j, 1 + k, 0 ) = matrix[k];
					}

					akareferences[i].getVolumes()[j].setTransform(mat3);
					this->input[ j ].setTransform( mat3 );

					//// retrieve previous alignment
					//double matrix[16];
					//for ( int k = 0; k < 16; k++ ){
					//	matrix[k] = newAlignments(j,1+k,0);
					//}
					//vtkMatrix4x4 * mat1 = vtkMatrix4x4 :: New();
					//mat1->DeepCopy( matrix );
					//this->input[ j ].setTransform( mat1 );
					//mat1->Delete();

					//// check if transformed volume has still acceptable transformation
					//vtkTransform * test = vtkTransform :: New();
					//test->Concatenate( mat3 );
					//if ( this->metric->isTransformValid( test ) == false ){
					//	cerr << "ERROR: transformation should always be valid. File " << __FILE__ << ": " << __LINE__ << endl;
					//	cerr << "\tR=[" << test->GetOrientation()[0] << "," << test->GetOrientation()[1] << "," << test->GetOrientation()[2] << "],	t=[" << test->GetPosition()[0] << "," << test->GetPosition()[1] << "," << test->GetPosition()[2] << "]" << endl;
					//}
					//test->Delete();

					mat1->Delete();
					mat2->Delete();
					mat3->Delete();

					count++;
				}
			}

			// set new weights as sigmoid function
			//if ( iterations == this->refinementIterations - 1 ){
				// sort distances to global average to determine best top percentage
				sort( distancesToGlobalAverage.begin(), distancesToGlobalAverage.end() );
				Pixel distanceTh = distancesToGlobalAverage[ this->distanceTopCutoffPre * ( distancesToGlobalAverage.size() - 1 ) ];
				// distanceTh = ( distancesToGlobalAverage[0] + distancesToGlobalAverage[ distancesToGlobalAverage.size() - 1 ] ) * this->distanceTopCutoffPre;

				// sigmoid decay
				Pixel decay = .05;

				newAlignments( Range::all(), 0, 0 ) = where( distancesToGlobalAverageArray <  numeric_limits< Pixel > :: max(), 1, 0 ) * ( 
					( 1 + 2.0 * exp( - ( distancesToGlobalAverageArray - distanceTh ) / decay ) ) /
					( 1 + 1.0 * exp( - ( distancesToGlobalAverageArray - distanceTh ) / decay ) ) - 1.0 );
				// eliminate weights below 1e-2
				newAlignments( Range::all(), 0, 0 ) = where( newAlignments( Range::all(), 0, 0 ) > 1e-2, newAlignments( Range::all(), 0, 0 ), 0 );

				cout << "Distances to reference " << i << " = [" << distancesToGlobalAverage[0] << "," << distancesToGlobalAverage[distancesToGlobalAverage.size() - 1] << "], th = " << distanceTh << ". Left out " << distancesToGlobalAverage.size()-sum(where( newAlignments( Range::all(), 0, 0 ) > 0, 1, 0 )) << " of " << distancesToGlobalAverage.size() << endl; 
			//}

			akareferences[i].setAlignments( newAlignments );

			// update classes
			this->classes( i, Range::all(), 0 ) = newAlignments( Range::all(), 0, 0 );
		}
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: alignReferencesToStrongest( bool apply_rotational_symmetry )
{
	// align all references to the one with most volumes
	vector< nbfWedgedImage3D< Pixel > * > list1;
	vector< nbfWedgedImage3D< Pixel > * > list2;

	// find reference with most volumes (assume the first one is always used as reference)
	int bestReference = 0;
	Pixel maxSoFar = 0;

	//	for ( int i = 0; i < this->references.size(); i++ ){
	//		Pixel currMax = sum( this->references[i].weights( Range::all(), 0 ) );
	//		if ( currMax > maxSoFar ){
	//			bestReference = i;
	//			maxSoFar = currMax;
	//		}
	//	}

	// before aligning everybody to it, center volume in frame

	//// xy-plane centering with vertical cylinder
	//if ( this->metric->getRotationSearchRestriction() < 180 ){

	//	// center reference image using CCC with vertical cylinder
	//	nbfFourierImageMetric< Pixel, 3 > myMetric( this->metric->imageFilter, this->metric->fourierFilter );
	//	//myMetric.setTranslationSearchRestriction( this->metric->getTranslationSearchRestriction() );
	//	myMetric.setNumberOfCandidates( 1 );

	//	// alignment reference
	//	Array< Pixel, 3 > cylinder( this->input[0].getDimensions() );
	//	firstIndex i; secondIndex j;

	//	// using generic cylinder for alignment
	//	cylinder = sqrt( pow2( i - ( cylinder.rows() - 1 ) / 2.0 ) + pow2( j - ( cylinder.cols() - 1 ) / 2.0 ) );
	//	nbfWedgedSubImage3D< Pixel > reference;
	//	reference.setFixedImage( cylinder );

	//	myMetric.setInput1( &reference );

	//	myMetric.setInput2( &this->references[bestReference] );

	//	vtkTransform * t = vtkTransform :: New();
	//	myMetric.executeFourierNewHalf(t);

	//	// extract shifts
	//	double shifts[3];
	//	t->GetPosition(shifts);
	//	t->Delete();

	//	// only correct in-plane shifts (ignore z component)
	//	shifts[2] = 0.0;
	//	r->Translate( shifts );

	//// xyz-rotational centering
	//} else {
	if ( this->metric->getRotationSearchRestriction() == 180 ){ // GROEL - ONLY

		// we need to do rotational matching
		nbfProjectionRotationMetric3D< Pixel > myRotMetric( this->metric->imageFilter, this->metric->fourierFilter );
		myRotMetric.setNumberOfCandidatePeaksToSearch( 15 );
		myRotMetric.setRotationSearchRestriction( this->metric->getRotationSearchRestriction() );
		myRotMetric.setTranslationSearchRestriction( this->metric->getTranslationSearchRestriction() );
		myRotMetric.setToComputeOverlapNormalizedDistances( this->metric->getToComputeOverlapNormalizedDistances() );
		myRotMetric.setToUseMutualCorrelation( this->metric->getToUseMutualCorrelation() );
		myRotMetric.setMissingWedgeCompensation( this->metric->getMissingWedgeCompensation() );
		// myRotMetric.imageFilter->maskOn( this->metric->imageFilter->maskFile );
		myRotMetric.setNumberOfCandidates( 1 );

		// alignment reference
		Array< Pixel, 3 > cylinder( this->input[0].getDimensions() );
		firstIndex i; secondIndex j;

		// using (inverted) cylinder for alignment
		cylinder = sqrt( pow2( i - ( cylinder.rows() - 1 ) / 2.0 ) + pow2( j - ( cylinder.cols() - 1 ) / 2.0 ) );
		Pixel inner = 10;
		Pixel outer = 16;
		cylinder = where( (cylinder > inner ) && ( cylinder < outer ), 0, 1 );
		int height = floor( cylinder.depth() / 2.0 );
		int halfHeight = 16;
		cylinder( Range :: all(), Range :: all(), Range( fromStart, height - halfHeight ) ) = 1;
		cylinder( Range :: all(), Range :: all(), Range( height + halfHeight, toEnd ) ) = 1;

		//w.write(cylinder);

		nbfWedgedSubImage3D< Pixel > reference;
		reference.setFixedImage( cylinder );

		myRotMetric.setInput1( &reference );

		// set current best volume as second input
		myRotMetric.setInput2( &this->references[bestReference] );

		myRotMetric.execute();

		vtkTransform * r = vtkTransform :: New();
		r->SetMatrix( myRotMetric.getTransform(0)->GetMatrix() );

		// re-compute master reference with new alignments

		Array< Pixel, 3 > newAlignments( this->input.size(), 17, this->references[bestReference].weights.cols() );
		newAlignments = 0;
		newAlignments( Range::all(), 0, Range::all() ) = this->references[bestReference].weights;

		// retrieve alignment of current reference to common frame
		//vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
		//matrix2->DeepCopy( r->GetMatrix() );

		//cout << "Centering strongest reference transformation = " << *r->GetMatrix() << endl;

		// concatenate new transform to existing alignments
		for ( int j = 0; j < this->references[bestReference].getVolumes().size(); j++ ){

			if ( this->references[bestReference].weights(j,0) > 0 ){

				// retrieve original transformation
				double matrix[16];
				for ( int k = 0; k < 16; k++ ){
					matrix[k] = this->references[bestReference].multipleAlignments( j, k, 0 );
				}
				vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
				matrix1->DeepCopy( matrix );

				// compose the two
				vtkMatrix4x4 * matrix3 = vtkMatrix4x4 :: New();
				vtkMatrix4x4 :: Multiply4x4( matrix1, r->GetMatrix(), matrix3 );

				vtkMatrix4x4 :: DeepCopy( matrix, matrix3 );

				// update volume list alignment
				this->input[j].setTransform( matrix3 );

				matrix1->Delete();

				for ( int k = 0; k < 16; k++ ){
					newAlignments( j, 1 + k, 0 ) = matrix[k];
				}

				//// check if transformed volume has still acceptable transformation
				//vtkTransform * test = vtkTransform :: New();
				//test->Concatenate( matrix3 );
				//if ( this->metric->isTransformValid( test ) == false ){
				//	// volume is not longer valid, we remove from current average
				//	// newAlignments( j, 0, 0 ) = 0;
				//	cerr << "WARNING: Volume " << this->input[j].getFileName() << " does not satisfy constraints after MAIN reference re-centering. File " << __FILE__ << ": " << __LINE__ << endl;
				//}
				//test->Delete();

				matrix3->Delete();
			}
		}

		// set new alignments and re-compute average
		this->references[bestReference].setAlignments( newAlignments );
		this->metric->getImage( this->references[bestReference] );

		// cleanup
		//matrix2->Delete();
		r->Delete();

		// vertical centering

		nbfFourierImageMetric< Pixel, 3 > myMetric( this->metric->imageFilter, this->metric->fourierFilter );
		myMetric.setNumberOfCandidates( 1 );
		myMetric.setMissingWedgeCompensation( this->metric->getMissingWedgeCompensation() );

		// set current best volume as first input
		myMetric.setInput1( &this->references[bestReference] );

		// build its mirror image to use as second input
		nbfWedgedAverageImage3D< Pixel > mirrorImage( this->input );
		newAlignments.resize( this->input.size(), 17, this->metric->getNumberOfCandidates() );
		newAlignments( Range::all(), 0, Range::all() ) = this->references[bestReference].weights;
		newAlignments( Range::all(), Range(1,toEnd), Range::all() ) = this->references[bestReference].multipleAlignments;

		vtkTransform * concat = vtkTransform :: New();
		concat->Scale( 1, 1, -1 );

		// apply mirror transformation to average
		for ( int i = 0; i < this->input.size(); i++ ){
			if ( this->references[bestReference].weights(i,0) > 0 ){
				// retrieve original transformation
				double matrix[16];
				for ( int k = 0; k < 16; k++ ){
					matrix[k] = this->references[bestReference].multipleAlignments( i, k, 0 );
				}
				vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
				matrix1->DeepCopy( matrix );
				vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
				vtkMatrix4x4 :: Multiply4x4( matrix1, concat->GetMatrix(), matrix2 );
				vtkMatrix4x4 :: DeepCopy( matrix, matrix2 );
				for ( int k = 0; k < 16; k++ ){
					newAlignments( i, 1 + k, 0 ) = matrix[k];
				}
				matrix1->Delete();
				matrix2->Delete();
			}
		}
		mirrorImage.setAlignments( newAlignments );
		this->metric->getImage( mirrorImage );
		concat->Delete();

		// set second input
		myMetric.setInput2( &mirrorImage );

		vtkTransform * t = vtkTransform :: New();
		myMetric.executeFourierNewHalf(t);

		// extract shifts
		double shifts[3];
		t->GetPosition(shifts);
		t->Delete();

		// Keep only vertical displacement
		shifts[0] /= 2.0; shifts[1] /= 2.0; shifts[2] /= 2.0;

		vtkTransform * result1 = vtkTransform :: New();
		result1->Translate( shifts );

		// apply transformation to average
		for ( int i = 0; i < this->references[bestReference].getVolumes().size(); i++ ){
			if ( this->references[bestReference].weights(i,0) > 0 ){
				// retrieve original transformation
				double matrix[16];
				for ( int k = 0; k < 16; k++ ){
					matrix[k] = this->references[bestReference].multipleAlignments( i, k, 0 );
				}
				vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
				matrix1->DeepCopy( matrix );
				vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
				vtkMatrix4x4 :: Multiply4x4( matrix1, result1->GetMatrix(), matrix2 );
				vtkMatrix4x4 :: DeepCopy( matrix, matrix2 );
				for ( int k = 0; k < 16; k++ ){
					newAlignments( i, 1 + k, 0 ) = matrix[k];
				}
				matrix1->Delete();

				//// check if transformed volume has still acceptable transformation, otherwise remove volume from average
				//vtkTransform * test = vtkTransform :: New();
				//test->Concatenate( matrix2 );
				//if ( this->metric->isTransformValid( test ) == false ){
				//	// volume is not longer valid, we remove from current average
				//	// newAlignments( i, 0, 0 ) = 0;
				//	cerr << "WARNING: Volume " << this->input[i].getFileName() << " does not satisfy constraints after vertical re-centering. File " << __FILE__ << ": " << __LINE__ << endl;
				//}
				//test->Delete();

				matrix2->Delete();
			}
		}
		this->references[bestReference].setAlignments( newAlignments );
		this->metric->getImage( this->references[bestReference] );
		result1->Delete();
	}
 //else if ( apply_rotational_symmetry == true ){
	//	cout << " Aligning first reference to its rotationally symmetrized copy." << endl;
	//	// align best reference to symmetrized average

	//	// apply rotational symmetry
	//	Array< double, 3 > in, out, res;
	//	Array< bool, 3 > zone;
	//	vtkImageData * data = vtkImageData :: New();
	//	this->references[bestReference].getImage( data );
	//	nbfVTKInterface :: vtkToBlitzReference( data, in );
	//	res.resize( in.shape() );
	//	res = in;
	//	data->Delete();

	//	nbfCylindricalDomain3< double > cyl;
	//	TinyVector< int, 3 > center = in.shape() / 2;
	//	cyl.setCenter( center );
	//	cyl.setMaxRho( in.extent(firstDim) / 2 );
	//	cyl.setResRho( in.extent(firstDim) / 2 );
	//	cyl.setResTheta(180);
	//	cyl.cartesian2cylindrical( res, out, zone );

	//	// compute radial profile
	//	Array< double, 2 > profile( out.rows(), out.depth() );
	//	firstIndex i; secondIndex j; thirdIndex k;
	//	profile = sum( out(i,k,j), k );
	//	// the center axis does not get averaging, so we average with the closest ring outside.
	//	for ( int i = 0; i < out.cols(); i++ ){
	//		out( Range :: all(), i, Range :: all() ) = profile;
	//	}

	//	// back-transform to cartesian coordinates
	//	cyl.cylindrical2cartesian( out, res, zone );
	//	Array< Pixel, 3 > cylindricalAverage( res.shape() );
	//	cylindricalAverage = cast<Pixel>(res);

	//	// metric for doing the alignment
	//	nbfProjectionRotationMetric3D< Pixel > myRotMetric( this->metric->imageFilter, this->metric->fourierFilter );
	//	myRotMetric.setNumberOfCandidatePeaksToSearch( 25 );
	//	myRotMetric.setRotationSearchRestriction( this->metric->getRotationSearchRestriction() );
	//	myRotMetric.setTranslationSearchRestriction( this->metric->getTranslationSearchRestriction() );
	//	myRotMetric.setToComputeOverlapNormalizedDistances( this->metric->getToComputeOverlapNormalizedDistances() );
	//	myRotMetric.setToUseMutualCorrelation( this->metric->getToUseMutualCorrelation() );
	//	myRotMetric.setMissingWedgeCompensation( this->metric->getMissingWedgeCompensation() );
	//	myRotMetric.setNumberOfCandidates( 1 );

	//	// alignment reference
	//	nbfWedgedSubImage3D< Pixel > reference;
	//	reference.setFixedImage( cylindricalAverage );
	//	myRotMetric.setInput1( &reference );

	//	// set current best volume as second input
	//	myRotMetric.setInput2( &this->references[bestReference] );

	//	myRotMetric.execute();

	//	vtkTransform * r = vtkTransform :: New();
	//	r->SetMatrix( myRotMetric.getTransform(0)->GetMatrix() );

	//	// re-compute master reference with new alignments

	//	Array< Pixel, 3 > newAlignments( this->input.size(), 17, this->references[bestReference].weights.cols() );
	//	newAlignments = 0;
	//	newAlignments( Range::all(), 0, Range::all() ) = this->references[bestReference].weights;

	//	// concatenate new transform to existing alignments
	//	for ( int j = 0; j < this->references[bestReference].getVolumes().size(); j++ ){

	//		if ( this->references[bestReference].weights(j,0) > 0 ){

	//			// retrieve original transformation
	//			double matrix[16];
	//			for ( int k = 0; k < 16; k++ ){
	//				matrix[k] = this->references[bestReference].multipleAlignments( j, k, 0 );
	//			}
	//			vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
	//			matrix1->DeepCopy( matrix );

	//			// compose the two
	//			vtkMatrix4x4 * matrix3 = vtkMatrix4x4 :: New();
	//			vtkMatrix4x4 :: Multiply4x4( matrix1, r->GetMatrix(), matrix3 );

	//			vtkMatrix4x4 :: DeepCopy( matrix, matrix3 );

	//			// update volume list alignment
	//			this->input[j].setTransform( matrix3 );

	//			matrix1->Delete();

	//			for ( int k = 0; k < 16; k++ ){
	//				newAlignments( j, 1 + k, 0 ) = matrix[k];
	//			}

	//			matrix3->Delete();
	//		}
	//	}

	//	// set new alignments and re-compute average
	//	this->references[bestReference].setAlignments( newAlignments );
	//	this->metric->getImage( this->references[bestReference] );

	//	// cleanup
	//	r->Delete();
	//} else {
	//	// bundle align all references
	//	this->bundleAlignment( this->references );

	//	// align all to external reference to make membrane stand vertically
	//}

	// align remaining references to this one
	list1.push_back( &(this->references[bestReference]) );
	vector< TinyVector< int, 2 > > pos;
	int count = 0;
	for ( int i = 0; i < this->references.size(); i++ ){
		if ( i != bestReference ){
			list2.push_back( &(this->references[i]) );
			TinyVector< int, 2 > p(0,count);
			pos.push_back(p);
			count++;
		}
	}

	Array< Pixel, 3 > allAlignments( list2.size(), 17, 1 );
	allAlignments = -1;

	//// restrict search to in-plane rotation only
	//Pixel currentRotationRestriction =  this->metric->getRotationSearchRestriction();
	//if ( currentRotationRestriction < 180 ){
	//	this->metric->setRotationSearchRestriction( 0.0 );
	//}

	this->metric->getDistances( list1, list2, pos, allAlignments, 0 );
	
	//// restore current rotation search restriction
	//if ( currentRotationRestriction < 180 ){
	//	this->metric->setRotationSearchRestriction( currentRotationRestriction );
	//}

	// re-compute averages with new transformations

	count = 0;

	for ( int i = 0; i < this->references.size(); i++ ){

		if ( i != bestReference ){

			Array< Pixel, 3 > newAlignments( this->references[i].getVolumes().size(), 17, this->references[i].weights.cols() );
			newAlignments = 0;
			newAlignments( Range::all(), 0, Range::all() ) = this->references[i].weights;

			// retrieve alignment of current reference to common frame
			vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
			double matrix[16];

			for ( int k = 0; k < 16; k++ ){
				matrix[k] = allAlignments( count, 1 + k, 0 );
			}
			matrix2->DeepCopy( matrix );
			count++;

			// concatenate new transform to existing alignments
			for ( int j = 0; j < this->references[i].getVolumes().size(); j++ ){

				if ( this->references[i].weights(j,0) > 0 ){

					// retrieve original transformation
					for ( int k = 0; k < 16; k++ ){
						matrix[k] = this->references[i].multipleAlignments( j, k, 0 );
					}
					vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
					matrix1->DeepCopy( matrix );

					// compose the two
					vtkMatrix4x4 * matrix3 = vtkMatrix4x4 :: New();
					vtkMatrix4x4 :: Multiply4x4( matrix1, matrix2, matrix3 );

					vtkMatrix4x4 :: DeepCopy( matrix, matrix3 );

					// update volume list alignment
					this->input[j].setTransform( matrix3 );

					matrix1->Delete();

					for ( int k = 0; k < 16; k++ ){
						newAlignments( j, 1 + k, 0 ) = matrix[k];
					}

					//// check if transformed volume has still acceptable transformation
					//vtkTransform * test = vtkTransform :: New();
					//test->Concatenate( matrix3 );
					//if ( this->metric->isTransformValid( test ) == false ){
					//	// volume is not longer valid, we remove from current average
					//	// newAlignments( j, 0, 0 ) = 0;
					//	cerr << "WARNING: Volume " << this->input[j].getFileName() << " does not satisfy constraints after reference re-centering. File " << __FILE__ << ": " << __LINE__ << endl;
					//}
					//test->Delete();

					matrix3->Delete();
				}
			}
			matrix2->Delete();

			// set alignments and re-compute average
			this->references[i].setAlignments( newAlignments );
			//this->metric->getImage( this->references[i] );
		}
	}

	// condense class averages
	this->classes.resize( this->classesSelected.size(), this->classes.cols(), this->classes.depth() );
	Array< Pixel, 4 > newAlignments( this->input.size(), 17, this->references[0].weights.cols(), this->classesSelected.size() );
	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		newAlignments( Range :: all(), Range :: all(), Range :: all(), i ) = 0;
		for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
			int currentIndex = this->classesSelected[i][j];
			for ( int k = 0; k < this->input.size(); k++ ){
				if ( this->references[ currentIndex ].weights(k,0) > 0 ){
					newAlignments( k, 0, 0, i ) = this->references[ currentIndex ].weights(k,0);
					newAlignments( k, Range(1,toEnd), 0, i ) = this->references[ currentIndex ].multipleAlignments( k, Range :: all(), 0 );
				}
			}
		}
		this->classes( i, Range :: all(), 0 ) = newAlignments( Range :: all(), 0, 0, i );
	}
	this->references.resize( this->classesSelected.size() );
	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		Array< Pixel, 3 > A( newAlignments( Range :: all(), Range :: all(), Range :: all(), i ) );
		this->references[i].setAlignments( A );
	}

	//if ( this->combineAllClasses == true ){
	//	Array< Pixel, 3 > newAlignments( this->references[0].getVolumes().size(), 17, this->references[0].weights.cols() );
	//	newAlignments = 0;
	//	for ( int i = 0; i < this->references.size(); i++ ){
	//		for ( int j = 0; j < this->references[i].getVolumes().size(); j++ ){
	//			if ( this->references[i].weights(j,0) > 0 ){
	//				newAlignments( j, 0, 0 ) = this->references[i].weights(j,0);
	//				newAlignments( j, Range(1,toEnd), 0 ) = this->references[i].multipleAlignments( j, Range :: all(), 0 );
	//			}
	//		}
	//	}
	//	// get rid of all classes except one
	//	while ( this->references.size() > 1 ){
	//		this->references.pop_back();
	//	}

	//	// update classes
	//	this->classes.resize( this->references.size(), this->classes.cols(), this->classes.depth() );
	//	this->classes( 0, Range :: all(), 0 ) = newAlignments( Range :: all(), 0, 0 );

	//	this->references[0].setAlignments( newAlignments );
	//}

	// compute new class averages
	for ( int i = 0; i < this->references.size(); i++ ){
		this->metric->getImage( this->references[i] );
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: alignReferencesToStrongestBundle( bool apply_rotational_symmetry )
{
	// align all references to the selected one
	vector< nbfWedgedImage3D< Pixel > * > list1;
	vector< nbfWedgedImage3D< Pixel > * > list2;

	// always assume the first class is the one used as reference
	int bestReference = 0;
	
	// align all classes to selected reference
	list1.push_back( &(this->externalReferences[bestReference]) );
	vector< TinyVector< int, 2 > > pos;
	int count = 0;
	for ( int i = 0; i < this->externalReferences.size(); i++ ){
		if ( i != bestReference ){
			list2.push_back( &(this->externalReferences[i]) );
			TinyVector< int, 2 > p(0,count);
			pos.push_back(p);
			count++;
		}
	}

	Array< Pixel, 3 > allAlignments( list2.size(), 17, 1 );
	allAlignments = -1;

	this->metric->getDistances( list1, list2, pos, allAlignments, 0 );
	
	// re-compute averages with new transformations

	count = 0;

	Array< Pixel, 3 > newAlignments( this->references[bestReference].weights.rows(), 17, this->references[bestReference].weights.cols() );

	for ( int i = 0; i < this->references.size(); i++ ){

		if ( i != bestReference ){

			newAlignments = 0;
			newAlignments( Range::all(), 0, Range::all() ) = this->references[i].weights;

			// retrieve alignment of current reference to common frame
			vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
			double matrix[16];

			for ( int k = 0; k < 16; k++ ){
				matrix[k] = allAlignments( count, 1 + k, 0 );
			}
			matrix2->DeepCopy( matrix );
			count++;

			// concatenate new transform to existing alignments
			for ( int j = 0; j < this->references[i].weights.rows(); j++ ){

				if ( this->references[i].weights(j,0) > 0 ){

					// retrieve original transformation
					for ( int k = 0; k < 16; k++ ){
						matrix[k] = this->references[i].multipleAlignments( j, k, 0 );
					}
					vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
					matrix1->DeepCopy( matrix );

					// compose the two
					vtkMatrix4x4 * matrix3 = vtkMatrix4x4 :: New();
					vtkMatrix4x4 :: Multiply4x4( matrix1, matrix2, matrix3 );

					vtkMatrix4x4 :: DeepCopy( matrix, matrix3 );

					// update volume list alignment
					this->input[j].setTransform( matrix3 );

					matrix1->Delete();

					for ( int k = 0; k < 16; k++ ){
						newAlignments( j, 1 + k, 0 ) = matrix[k];
					}

					matrix3->Delete();
				}
			}
			matrix2->Delete();

			// set alignments and re-compute average
			this->references[i].setAlignments( newAlignments );
			//this->metric->getImage( this->references[i] );
		}
	}

	// global average
	nbfWedgedAverageImage3D< Pixel > globalAverage( this->input );
	
	// do multiple bundle iterations
	// this->refinementIterations
	for ( int iterations = 0; iterations < abs( this->refinementIterations ); iterations++ ){

		// compute new global average
		newAlignments = 0;
		for ( int i = 0; i < this->references.size(); i++ ){
			for ( int j = 0; j < this->references[i].weights.rows(); j++ ){
				if ( this->references[i].weights(j,0) > 0 ){
					// set current weight
					newAlignments( j, 0, Range :: all() ) = this->references[i].weights(j,0);
					// set current alignments
					newAlignments( j, Range(1,toEnd), Range :: all() ) = this->references[i].multipleAlignments( j, Range::all(), Range :: all() );
				}
			}
		}
		// set alignments and re-compute average
		globalAverage.setAlignments( newAlignments );
		this->metric->getImage( globalAverage );

		list1.clear();
		list1.push_back( &globalAverage );

		pos.clear();
		list2.clear();

		for ( int i = 0; i < this->references.size(); i++ ){
			this->metric->getImage( this->references[i] );
			list2.push_back( &(this->references[i]) );
			TinyVector< int, 2 > p(0,i);
			pos.push_back(p);
		}

		// align all to global average
		allAlignments.resize( list2.size(), 17, 1 );
		allAlignments = -1;

		// alignment with local refinement only
		this->metric->getDistances( list1, list2, pos, allAlignments, 0 );

		// cout << "Distances = " << allAlignments( Range :: all(), 0, Range :: all() ) << endl;
		cout << " Bundle score at iteration " << iterations << " = " << sum( allAlignments( Range :: all(), 0, Range :: all() ) ) << endl;
		// re-compute averages with new transformations

		for ( int i = 0; i < this->references.size(); i++ ){

			newAlignments = 0;
			newAlignments( Range::all(), 0, Range::all() ) = this->references[i].weights;

			// retrieve alignment of current reference to common frame
			vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
			double matrix[16];

			for ( int k = 0; k < 16; k++ ){
				matrix[k] = allAlignments( i, 1 + k, 0 );
			}
			matrix2->DeepCopy( matrix );
			
			// concatenate new transform to existing alignments
			for ( int j = 0; j < this->references[i].weights.rows(); j++ ){

				if ( this->references[i].weights(j,0) > 0 ){

					// retrieve original transformation
					for ( int k = 0; k < 16; k++ ){
						matrix[k] = this->references[i].multipleAlignments( j, k, 0 );
					}
					vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
					matrix1->DeepCopy( matrix );

					// compose the two
					vtkMatrix4x4 * matrix3 = vtkMatrix4x4 :: New();
					vtkMatrix4x4 :: Multiply4x4( matrix1, matrix2, matrix3 );

					vtkMatrix4x4 :: DeepCopy( matrix, matrix3 );

					// update volume list alignment
					this->input[j].setTransform( matrix3 );

					matrix1->Delete();

					for ( int k = 0; k < 16; k++ ){
						newAlignments( j, 1 + k, 0 ) = matrix[k];
					}

					matrix3->Delete();
				}
			}
			matrix2->Delete();

			// set alignments and re-compute average
			this->references[i].setAlignments( newAlignments );
		}
	}

	//==================
	// Center global density in x-y plane

	//// compute new global average
	//newAlignments = 0;
	//for ( int i = 0; i < this->references.size(); i++ ){
	//	for ( int j = 0; j < this->references[i].weights.rows(); j++ ){
	//		if ( this->references[i].weights(j,0) > 0 ){
	//			// set current weight
	//			newAlignments( j, 0, Range :: all() ) = this->references[i].weights(j,0);
	//			// set current alignments
	//			newAlignments( j, Range(1,toEnd), Range :: all() ) = this->references[i].multipleAlignments( j, Range::all(), Range :: all() );
	//		}
	//	}
	//}
	//// set alignments and re-compute average
	//globalAverage.setAlignments( newAlignments );
	//this->metric->getImage( globalAverage );

	nbfFourierImageMetric< Pixel, 3 > myMetric( this->metric->imageFilter, this->metric->fourierFilter );
	myMetric.setNumberOfCandidates( 1 );
	myMetric.setMissingWedgeCompensation( this->metric->getMissingWedgeCompensation() );
	myMetric.setTranslationSearchRestriction( this->metric->getTranslationSearchRestriction() );

	// set current best volume as first input
	myMetric.setInput1( &this->references[bestReference] );

	// build its mirror image to use as second input
	nbfWedgedAverageImage3D< Pixel > mirrorImage( this->input );
	////Array< Pixel, 3 > newAlignmentsCenter( this->input.size(), 17, this->metric->getNumberOfCandidates() );
	//newAlignments( Range::all(), 0, Range::all() ) = this->references[bestReference].weights;
	//newAlignments( Range::all(), Range(1,toEnd), Range::all() ) = this->references[bestReference].multipleAlignments;

	vtkTransform * concat = vtkTransform :: New();
	concat->Scale( -1, -1, 1 );

	// apply mirror transformation to average
	newAlignments = 0;
	newAlignments( Range::all(), 0, Range::all() ) = this->references[bestReference].weights;
	for ( int i = 0; i < this->references[bestReference].weights.rows(); i++ ){
		if ( this->references[bestReference].weights(i,0) > 0 ){
			// retrieve original transformation
			double matrix[16];
			for ( int k = 0; k < 16; k++ ){
				matrix[k] = this->references[bestReference].multipleAlignments( i, k, 0 );
			}
			vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
			matrix1->DeepCopy( matrix );
			vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
			vtkMatrix4x4 :: Multiply4x4( matrix1, concat->GetMatrix(), matrix2 );
			vtkMatrix4x4 :: DeepCopy( matrix, matrix2 );
			for ( int k = 0; k < 16; k++ ){
				newAlignments( i, 1 + k, 0 ) = matrix[k];
			}
			matrix1->Delete();
			matrix2->Delete();
		}
	}
	concat->Delete();
	mirrorImage.setAlignments( newAlignments );
	this->metric->getImage( mirrorImage );

	// set second input
	myMetric.setInput2( &mirrorImage );
	//myMetric.setInput2( &this->references[bestReference] );

	vtkTransform * t = vtkTransform :: New();
	myMetric.executeFourierNewHalf(t);

	// extract shifts
	double shifts[3];
	t->GetPosition(shifts);
	t->Delete();

	// Keep only x-y plane displacements
	shifts[0] /= 2.0; shifts[1] /= 2.0; shifts[2] = 0.0;

	vtkTransform * result1 = vtkTransform :: New();
	result1->Translate( shifts );

	// apply transformation to all averages
	for ( int r = 0; r < this->references.size(); r++ ){
		newAlignments( Range::all(), 0, Range::all() ) = this->references[r].weights;
		newAlignments( Range::all(), Range(1,toEnd), Range::all() ) = this->references[r].multipleAlignments;
		for ( int i = 0; i < this->references[r].weights.rows(); i++ ){
			if ( this->references[r].weights(i,0) > 0 ){
				// retrieve original transformation
				double matrix[16];
				for ( int k = 0; k < 16; k++ ){
					matrix[k] = this->references[r].multipleAlignments( i, k, 0 );
				}
				vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
				matrix1->DeepCopy( matrix );
				vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
				vtkMatrix4x4 :: Multiply4x4( matrix1, result1->GetMatrix(), matrix2 );
				vtkMatrix4x4 :: DeepCopy( matrix, matrix2 );
				for ( int k = 0; k < 16; k++ ){
					newAlignments( i, 1 + k, 0 ) = matrix[k];
				}
				matrix1->Delete();
				matrix2->Delete();
			}
		}
		this->references[r].setAlignments( newAlignments );
		// this->metric->getImage( this->references[r] );
	}
	result1->Delete();

	//==================

	// condense class averages
	this->classes.resize( this->classesSelected.size(), this->classes.cols(), this->classes.depth() );
	Array< Pixel, 4 > finalAlignments( this->input.size(), 17, this->references[0].weights.cols(), this->classesSelected.size() );
	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		finalAlignments( Range :: all(), Range :: all(), Range :: all(), i ) = 0;
		for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
			int currentIndex = this->classesSelected[i][j];
			for ( int k = 0; k < this->input.size(); k++ ){
				if ( this->references[ currentIndex ].weights(k,0) > 0 ){
					finalAlignments( k, 0, 0, i ) = this->references[ currentIndex ].weights(k,0);
					finalAlignments( k, Range(1,toEnd), 0, i ) = this->references[ currentIndex ].multipleAlignments( k, Range :: all(), 0 );
				}
			}
		}
		this->classes( i, Range :: all(), 0 ) = finalAlignments( Range :: all(), 0, 0, i );
	}
	this->references.resize( this->classesSelected.size() );
	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		Array< Pixel, 3 > A( finalAlignments( Range :: all(), Range :: all(), Range :: all(), i ) );
		this->references[i].setAlignments( A );
	}

	// compute new class averages
	for ( int i = 0; i < this->references.size(); i++ ){
		this->metric->getImage( this->references[i] );
	}
}

// Use already computed averages instead of recomputing each time
template< class Pixel >
void nbfLoopClustering< Pixel > :: alignReferencesToStrongestBundleExternal( bool apply_rotational_symmetry )
{
	// align all references to the selected one
	vector< nbfWedgedImage3D< Pixel > * > list1;
	vector< nbfWedgedImage3D< Pixel > * > list2;

	// always assume the first class is the one used as reference
	int bestReference = 0;
	
	// align all classes to selected reference
	list1.push_back( &(this->externalReferences[bestReference]) );
	vector< TinyVector< int, 2 > > pos;
	int count = 0;
	for ( int i = 0; i < this->externalReferences.size(); i++ ){
		list2.push_back( &(this->externalReferences[i]) );
		TinyVector< int, 2 > p(0,i);
		pos.push_back(p);
	}

	Array< Pixel, 3 > allAlignments( list2.size(), 17, 1 );
	allAlignments = -1;

	// align all references to first component in vector
	this->metric->getDistances( list1, list2, pos, allAlignments, 0 );
	
	// store global average
	nbfWedgedAverageImage3D< Pixel > globalAverage( this->externalReferences );
	
	Array< Pixel, 3 > newAlignments( this->externalReferences.size(), 17, 1 );
	
	// set weights = 1 for global average
	newAlignments( Range :: all(), 0, Range :: all() ) = 1;
	
	// bundle alignment
	for ( int iterations = 0; iterations < abs( this->refinementIterations ); iterations++ ){

		// update transformations
		for ( int i = 0; i < this->externalReferences.size(); i++ ){

			// correction transformation
			double matrix[16];
			for ( int k = 0; k < 16; k++ ){
				matrix[k] = allAlignments( i, 1 + k, 0 );
			}
			vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
			matrix2->DeepCopy( matrix );
			
			// update transformation
			this->externalReferences[i].setTransform( matrix2 );
			
			newAlignments( i, Range(1,toEnd), 0 ) = allAlignments( i, Range(1,toEnd), 0 );

			matrix2->Delete();
		}

		globalAverage.setAlignments( newAlignments );
		this->metric->getImage( globalAverage );

		list1.clear();
		list1.push_back( &globalAverage );

		pos.clear();
		list2.clear();

		for ( int i = 0; i < this->externalReferences.size(); i++ ){
			list2.push_back( &(this->externalReferences[i]) );
			TinyVector< int, 2 > p(0,i);
			pos.push_back(p);
		}

		// align all to global average
		allAlignments.resize( list2.size(), 17, 1 );
		allAlignments = -1;

		// alignment with local refinement only
		this->metric->getDistances( list1, list2, pos, allAlignments, 1 );

		cout << " Bundle score at iteration " << iterations << " = " << sum( allAlignments( Range :: all(), 0, Range :: all() ) ) << endl;
	}

	//// Center global density in x-y plane

	//nbfFourierImageMetric< Pixel, 3 > myMetric( this->metric->imageFilter, this->metric->fourierFilter );
	//myMetric.setNumberOfCandidates( 1 );
	//myMetric.setMissingWedgeCompensation( this->metric->getMissingWedgeCompensation() );
	//myMetric.setTranslationSearchRestriction( this->metric->getTranslationSearchRestriction() );

	// update global average and set as reference for alignment
	for ( int i = 0; i < this->externalReferences.size(); i++ ){

		// correction transformation
		double matrix[16];
		for ( int k = 0; k < 16; k++ ){
			matrix[k] = allAlignments( i, 1 + k, 0 );
		}
		vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
		matrix2->DeepCopy( matrix );

		//newAlignments( i, Range(1,toEnd), 0 ) = allAlignments( i, Range(1,toEnd), 0 );

		// update transformation
		this->externalReferences[i].setTransform( matrix2 );
		matrix2->Delete();
	}
	//globalAverage.setAlignments( newAlignments );
	//this->metric->getImage( globalAverage );

	//myMetric.setInput1( &globalAverage );

	//nbfMatlabWriter w;
	//w.setFileName("w.matlab");
	//vtkImageData * data = vtkImageData::New();
	//globalAverage.getImage(data);
	//Array< double, 3 > A;
	//nbfVTKInterface::vtkToBlitzReference(data,A);
	//w.write(A);

	//// build its mirror image to use as second input
	//nbfWedgedAverageImage3D< Pixel > mirrorImage( this->externalReferences );

	//vtkTransform * concat = vtkTransform :: New();
	//concat->Scale( -1, -1, 1 );

	//// apply mirror transformation to all references
	//for ( int i = 0; i < this->externalReferences.size(); i++ ){
	//	double matrix[16];
	//	vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
	//	vtkTransform * t = vtkTransform :: New();
	//	this->externalReferences[i].getTransform(t);
	//	vtkMatrix4x4 :: Multiply4x4( t->GetMatrix(), concat->GetMatrix(), matrix2 );
	//	t->Delete();
	//	vtkMatrix4x4 :: DeepCopy( matrix, matrix2 );
	//	for ( int k = 0; k < 16; k++ ){
	//		newAlignments( i, 1 + k, 0 ) = matrix[k];
	//	}
	//	matrix2->Delete();
	//}
	//concat->Delete();
	//mirrorImage.setAlignments( newAlignments );
	//this->metric->getImage( mirrorImage );

	//mirrorImage.getImage(data);
	//nbfVTKInterface::vtkToBlitzReference(data,A);
	//w.write(A);


	//// set second input
	//myMetric.setInput2( &mirrorImage );
	//
	//vtkTransform * t = vtkTransform :: New();
	//myMetric.executeFourierNewHalf(t);

	//// extract shifts
	//double shifts[3];
	//t->GetPosition(shifts);
	//t->Delete();

	//// Keep only x-y plane displacements
	//shifts[0] /= 2.0; shifts[1] /= 2.0; shifts[2] = 0.0;

	//vtkTransform * result1 = vtkTransform :: New();
	//result1->Translate( shifts );

	//// compute final transformations
	//for ( int i = 0; i < this->externalReferences.size(); i++ ){
	//	double matrix[16];
	//	vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
	//	vtkTransform * t = vtkTransform :: New();
	//	this->externalReferences[i].getTransform(t);
	//	vtkMatrix4x4 :: Multiply4x4( t->GetMatrix(), result1->GetMatrix(), matrix2 );
	//	t->Delete();
	//	this->externalReferences[i].setTransform( matrix2 );
	//	matrix2->Delete();
	//}
	//result1->Delete();

	//// update global transformations
	//newAlignments.resize( this->references[bestReference].weights.rows(), this->references[bestReference].weights.cols(), 1 );
	//newAlignments = 0;

	//for ( int r = 0; r < this->references.size(); r++ ){
	//	newAlignments( Range::all(), 0, Range::all() ) = this->references[r].weights;
	//	for ( int i = 0; i < this->references[r].weights.rows(); i++ ){
	//		if ( this->references[r].weights(i,0) > 0 ){
	//			// retrieve original transformation
	//			double matrix[16];
	//			for ( int k = 0; k < 16; k++ ){
	//				matrix[k] = this->references[r].multipleAlignments( i, k, 0 );
	//			}
	//			vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
	//			matrix1->DeepCopy( matrix );
	//			vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
	//			vtkTransform * t = vtkTransform :: New();
	//			this->externalReferences[r].getTransform(t);
	//			vtkMatrix4x4 :: Multiply4x4( matrix1, t->GetMatrix(), matrix2 );
	//			t->Delete();
	//			vtkMatrix4x4 :: DeepCopy( matrix, matrix2 );
	//			for ( int k = 0; k < 16; k++ ){
	//				newAlignments( i, 1 + k, 0 ) = matrix[k];
	//			}
	//			matrix1->Delete();
	//			matrix2->Delete();
	//		}
	//	}
	//	this->references[r].setAlignments( newAlignments );
	//}

	//==================

	// compute pre-computed references
	Array< Pixel, 4 > finalAlignmentsExternal( this->externalReferences.size(), 17, 1, this->classesSelected.size() );
	finalAlignmentsExternal = 0;
	
	// assign most current alignments
	for ( int i = 0; i < this->externalReferences.size(); i++ ){
		
		vtkTransform * t = vtkTransform :: New();
		this->externalReferences[i].getTransform(t);
		
		double matrix[16];
		vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
		
		// set new alignments
		for ( int k = 0; k < 16; k++ ){
			finalAlignmentsExternal( i, 1 + k, 0, Range :: all() ) = matrix[k];
		}

		t->Delete();
	}

	// merge classes
	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
			int currentIndex = this->classesSelected[i][j];
			finalAlignmentsExternal( currentIndex, 0, 0, i ) = - 1.0 / this->classesSelected[i].size();
		}
	}

	this->precomputedReferences.resize( this->classesSelected.size() );
	for ( int i = 0; i < this->precomputedReferences.size(); i++ ){
		this->precomputedReferences[i].getVolumes() = this->externalReferences;
		Array< Pixel, 3 > A( finalAlignmentsExternal( Range :: all(), Range :: all(), Range :: all(), i ) );
		this->precomputedReferences[i].setAlignments( A );

		this->metric->getImage( this->precomputedReferences[i] );
	}

	//////////////////////
	
	// condense class averages
	this->classes.resize( this->classesSelected.size(), this->classes.cols(), this->classes.depth() );
	Array< Pixel, 4 > finalAlignments( this->input.size(), 17, this->references[0].weights.cols(), this->classesSelected.size() );
	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		finalAlignments( Range :: all(), Range :: all(), Range :: all(), i ) = 0;
		for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
			int currentIndex = this->classesSelected[i][j];
			for ( int k = 0; k < this->input.size(); k++ ){
				if ( this->references[ currentIndex ].weights(k,0) > 0 ){
					finalAlignments( k, 0, 0, i ) = this->references[ currentIndex ].weights(k,0);

					// retrieve corresponding correction
					vtkTransform * t = vtkTransform :: New();
					this->externalReferences[ currentIndex ].getTransform(t);

					// get previous alignments
					double matrix[16];
					for ( int m = 0; m < 16; m++ ){
						matrix[m] = this->references[ currentIndex ].multipleAlignments( k, m, 0 );
					}

					vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
					matrix1->DeepCopy( matrix );
					vtkMatrix4x4 * matrix2 = vtkMatrix4x4 :: New();
					vtkMatrix4x4 :: Multiply4x4( matrix1, t->GetMatrix(), matrix2 );

					// asign most current alignemnt to volume list
					this->input[ k ].setTransform( matrix2 );

					vtkMatrix4x4 :: DeepCopy( matrix, matrix2 );
					for ( int m = 0; m < 16; m++ ){
						finalAlignments( k, 1 + m, 0, i ) = matrix[m];
					}

					matrix2->Delete();
					matrix1->Delete();
					t->Delete();
				}
			}
		}
		this->classes( i, Range :: all(), 0 ) = finalAlignments( Range :: all(), 0, 0, i );
	}
	this->references.resize( this->classesSelected.size() );
	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		Array< Pixel, 3 > A( finalAlignments( Range :: all(), Range :: all(), Range :: all(), i ) );
		this->references[i].setAlignments( A );
	}
}


// Use already computed averages instead of recomputing each time
template< class Pixel >
void nbfLoopClustering< Pixel > :: alignReferencesInBundles( bool apply_rotational_symmetry, int alignment_mode )
{
	// 1 - align classes to reference within each group
	
	vector< nbfWedgedImage3D< Pixel > * > list1;
	vector< nbfWedgedImage3D< Pixel > * > list2;

	vector< TinyVector< int, 2 > > pos;

	// store copy of references so we can apply symmetry leaving the original references unchanged
	vector< nbfWedgedAverageImage3D< Pixel > > tmpAverage( this->classesSelected.size() );

	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		// always assume the first class is the one used as reference
		tmpAverage[i] = this->references[ this->classesSelected[i][0] ];
		list1.push_back( &tmpAverage[i] );
		//list1.push_back( &( this->references[ this->classesSelected[i][0] ] ) );
		if ( this->symmetryFactor > 1 ){
			nbfWedgedAverageImage3D< Pixel > * volume = reinterpret_cast< nbfWedgedAverageImage3D< Pixel > * >( list1[ list1.size() - 1 ] );
			this->symmetrize( *volume, this->symmetryFactor, false );
		}
		for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
			list2.push_back( &(this->references[ this->classesSelected[i][j] ] ) );
			TinyVector< int, 2 > p( i, list2.size() - 1 );
			pos.push_back(p);
		}
	}

	Array< Pixel, 3 > allAlignments( pos.size(), 17, 1 );
	allAlignments = -1;

	// create identity transformation
	vtkMatrix4x4 * mat = vtkMatrix4x4 :: New();
	mat->Identity();
	double matel[16];
	vtkMatrix4x4 :: DeepCopy( matel, mat );
	Array< Pixel, 1 > identity( 16 );
	for ( int i = 0; i < 16; i++ ){
		identity( i ) = matel[i];
	}
	mat->Delete();
	
	// initialize alignments to identity
	for ( int i = 0; i < allAlignments.rows(); i++ ){
		allAlignments( i, Range(1,toEnd), 0 ) = identity;
	}

	// If this is the initial iteration, we use the first reference in each group (-i) as starting reference
	if ( this->starting_iteration == 1 ){

		// first element in list1 and list2 are the same, so no need to compute those entries
		int counter = 0;
		for ( int i = 0; i < this->classesSelected.size(); i++ ){
			for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
				if ( j == 0 ){
					// set distance to 0 (avoid distance computation)
					allAlignments( counter, 0, 0 ) = 0;
				}
				counter++;
			}
		}

		// compute remaining alignments
		this->metric->getDistances( list1, list2, pos, allAlignments, alignment_mode );
		// cout << "Distances = " << allAlignments( Range :: all(), 0, Range :: all() ) << endl;
		cout << "Starting bundle score = " << sum( allAlignments( Range :: all(), 0, Range :: all() ) ) << endl;
	}

	// at this point, class averages have ONLY been aligned to the selected reference within each class

	// 2 - now iterate within each group independently (within group bundle alignment)

	// store global average for each group
	vector< nbfWedgedAverageImage3D< Pixel > > currentAverages;
	
	// build averages of class averages from written mrc files which have inverted contrast to avoid re-computation of class averages
	nbfWedgedAverageImage3D< Pixel > currentAverage( this->precomputedReferences );
	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		currentAverages.push_back( currentAverage );
	}

	list1.clear();

	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		list1.push_back( &( currentAverages[i] ) );
	}

	Array< Pixel, 3 > currentAlignments( this->precomputedReferences.size(), 17, 1 );

	int count;

	// bundle alignment
	for ( int iterations = 0; iterations < this->refinementIterations; iterations++ ){

		// update global average within group with most recently computed alignments
		
		count = 0;
		for ( int i = 0; i < this->classesSelected.size(); i++ ){
			currentAlignments = 0;
			for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
				// use negative weights because mrc files have inverted contrast
				currentAlignments( this->classesSelected[i][j], 0, 0 ) = - 1.0 / this->classesSelected[i].size();
				currentAlignments( this->classesSelected[i][j], Range(1,toEnd), Range :: all() ) = allAlignments( count++, Range(1,toEnd), Range :: all() );
			}
			nbfWedgedAverageImage3D< Pixel > * volume = reinterpret_cast< nbfWedgedAverageImage3D< Pixel > * >( list1[i] );
			volume->setAlignments( currentAlignments );
			this->metric->getImage( *volume );
			
			// symmetrize if needed
			if ( this->symmetryFactor > 1 ){
				this->symmetrize( *volume, this->symmetryFactor, false );
			}

			// update average image and wedge
			if ( this->metric->getMissingWedgeCompensation() == true ){
				volume->updateAccumulatedWedgeImage();
				if ( this->metric->getId() == NBF_IMAGE_METRIC_PROJECTION ){
					TinyVector< int, 2 > size = volume->getVolumesRO()[0].wedge.sphericalWedge.shape();
					volume->updateSphericalWedgeImage( size, reinterpret_cast< nbfProjectionRotationMetric3D< Pixel > * >(this->metric) );
				}
			}
		}

		// align volumes in group to global average
		allAlignments = -1;

		// skip computation if single-volume groups
		count = 0;
		for ( int i = 0; i < this->classesSelected.size(); i++ ){
			for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
				// keep the reference of each group unchanged
				if ( this->classesSelected[i].size() == 1 ){
					allAlignments( count, 0, 0 ) = 0;
					allAlignments( count, Range(1,toEnd), 0 ) = identity;
				}
				count++;
			}
		}

		// compute alignments
		this->metric->getDistances( list1, list2, pos, allAlignments, alignment_mode );
		cout << " Iteration " << iterations << " bundle score = " << sum( allAlignments( Range :: all(), 0, Range :: all() ) ) << endl;
	}

	// assign latest set of alignments
	count = 0;
	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		currentAlignments( Range :: all(), 0, Range :: all() ) = 0;
		for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
			currentAlignments( this->classesSelected[i][j], 0, 0 ) = - 1.0 / this->classesSelected[i].size();
			currentAlignments( this->classesSelected[i][j], Range(1,toEnd), Range :: all() ) = allAlignments( count++, Range(1,toEnd), Range :: all() );
		}
		nbfWedgedAverageImage3D< Pixel > * volume = reinterpret_cast< nbfWedgedAverageImage3D< Pixel > * >( list1[i] );
		volume->setAlignments( currentAlignments );

		this->metric->getImage( *volume );

		if ( this->metric->getMissingWedgeCompensation() == true ){
			Array< Pixel, 3 > W( volume->getDimensions() );
			volume->getWedgeImage( W );
			if ( this->metric->getId() == NBF_IMAGE_METRIC_PROJECTION ){
				Array< Pixel, 2 > Ws( volume->getVolumesRO()[0].wedge.sphericalWedge.shape() );
				volume->getSphericalWedgeImage( Ws, (vtkTransform*)NULL, reinterpret_cast< nbfProjectionRotationMetric3D< Pixel > * >(this->metric) );
			}
		}
	}

	// align all groups to first reference (-1 index)
	list2.clear();
	list2.push_back( &tmpAverage[0] );
	//list2.push_back( &this->references[ this->classesSelected[0][0] ] );
	
	pos.clear();
	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		TinyVector< int, 2 > p( 0, i );
		pos.push_back(p);
	}
	Array< Pixel, 3 > refAlignments( list1.size(), 17, 1 );
	refAlignments = -1;

	// skip first element which aligns the reference to itself anyways
	//refAlignments( 0, 0, 0 ) = 0;
	//refAlignments( 0, Range(1,toEnd), 0 ) = identity;
	this->metric->getDistances( list2, list1, pos, refAlignments, alignment_mode );

	// compute and apply final corrections
	for ( int i = 0; i < this->classesSelected.size(); i++ ){
		
		// retrieve alignment of group to -1 reference
		double matrix[16];
		for ( int k = 0; k < 16; k++ ){
			matrix[k] = refAlignments( i, 1 + k, 0 );
		}
		vtkMatrix4x4 * m = vtkMatrix4x4 :: New();
		m->DeepCopy( matrix );

		nbfWedgedAverageImage3D< Pixel > * volume = reinterpret_cast< nbfWedgedAverageImage3D< Pixel > * >( list1[i] );

		for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
		
			// retrieve alignment of class within group average
			for ( int k = 0; k < 16; k++ ){
				matrix[k] = volume->multipleAlignments( this->classesSelected[i][j], k, 0 );
			}
			vtkMatrix4x4 * m1 = vtkMatrix4x4 :: New();
			m1->DeepCopy( matrix );

			// compose with incremental alignment to global reference (first reference in first group)

			vtkMatrix4x4 * r = vtkMatrix4x4 :: New();
			vtkMatrix4x4 :: Multiply4x4( m1, m, r );
			m1->Delete();
	
			// rotate reference permanently
			this->references[ this->classesSelected[i][j] ].rotate( r, reinterpret_cast< nbfProjectionRotationMetric3D<Pixel>* >(this->metric) );
			r->Delete();
		}
		m->Delete();
	}

	// save most current alignments to input volume list
	for ( int i = 0; i < this->input.size(); i++ ){
		for ( int j = 0; j < this->references.size(); j++ ){
			if ( abs( this->references[j].weights(i,0) ) > 0 ){		
				// retrieve alignment of class within group average
				double matrix[16];
				for ( int k = 0; k < 16; k++ ){
					matrix[k] = this->references[j].multipleAlignments( i, k, 0 );
				}
				vtkMatrix4x4 * m1 = vtkMatrix4x4 :: New();
				m1->DeepCopy( matrix );
				this->input[i].setTransform( m1 );
				m1->Delete();
			}
		}
	}

	//// assign correlation to latest reference
	//list2.clear();
	//list2.push_back( &tmpAverage[0] );
	//refAlignments = -1;
	//this->metric->getDistances( list2, list1, pos, refAlignments, 2 );
	//for ( int i = 0; i < this->classesSelected.size(); i++ ){
	//	for ( int j = 0; j < this->classesSelected[i].size(); j++ ){
	//		this->references[ this->classesSelected[i][j] ].setCutOffset( refAlignemnts() );
	//	}
	//}
}


template< class Pixel >
void nbfLoopClustering< Pixel > :: alignReferencesToCommonFrame()
{
	// build volume list with all reference volumes
	vector< nbfWedgedSubImage3D< Pixel > > list;
	for ( int i = 0; i < this->references.size(); i++ ){
		list.push_back( this->references[i].getVolumes()[0]);
	}
	nbfWedgedAverageImage3D< Pixel > allvolumes( list );
	Array< Pixel, 3 > refAlignments( list.size(), 17, 1 );
	refAlignments = 1;
	double matrix[16];
	vtkMatrix4x4 :: Identity( matrix );
	for ( int i = 0; i < 16; i++ ){
		refAlignments( Range::all(), 1 + i, 0 ) = matrix[i];
	}
	allvolumes.setAlignments( refAlignments );

	vector< nbfWedgedAverageImage3D< Pixel > > allVolumesList;
	allVolumesList.push_back( allvolumes );

	// restric angular search to in plane rotations strictly
	
	// we need to store the current restriction angle
	Pixel storedPreviousAngle = this->metric->getRotationSearchRestriction();
	if ( storedPreviousAngle < 180 ){
		this->metric->setRotationSearchRestriction(0.0);
	}

	//this->bundleAlignment( allVolumesList, 2 );

	// restore previous restriction angle
	if ( storedPreviousAngle < 180 ){
		this->metric->setRotationSearchRestriction( storedPreviousAngle );
	}

	// assign new images to references
	for ( int i = 0; i < this->references.size(); i++ ){
		vtkImageData * data = vtkImageData::New();
		this->references[i].getImage(data);
		Array< double, 3 > A;
		nbfVTKInterface::vtkToBlitzReference( data, A );
		this->references[i].setImage( A );
		data->Delete();
	}

	// initialize new alignment matrix of references to volumes
	this->alignmentBetweenReferences.resize( this->classes.rows(), this->classes.rows(), 17, this->metric->getNumberOfCandidates() );
	this->alignmentBetweenReferences = -1;

	//nbfHierarchicalClustering< Pixel > hierarch;

	//Pixel cutoff = sum( this->alignmentBetweenReferences( Range::all(), Range::all(), 0, 0 ) ) / this->alignmentBetweenReferences.rows() / this->alignmentBetweenReferences.cols();
	//cout << "Using cutoff = " << cutoff << endl;
	//hierarch.setMaxDistance( cutoff );

	////hierarch.setMaxDistance( this->hierarchicalCutoff );
	//
	//hierarch.setMinOverlap( 0.0 );
	//hierarch.setMinElementNumber( 2 );

	////hierarch.setInput( subInput );
	//hierarch.setMetric( this->metric, this->alignmentBetweenReferences );
	//Array< Pixel, 3 > subClasses;
	//hierarch.execute( subClasses );

	//vector< vector< int > > classIndexes;
	//for ( int i = 0; i < subClasses.rows(); i++ ){
	//	vector< int > elements;	
	//	for ( int j = 0; j < subClasses.cols(); j++ ){
	//		if ( subClasses(i,j,0) == 1 ){
	//			elements.push_back(j);
	//		}
	//	}
	//	classIndexes.push_back( elements );
	//}

	//// align volumes within class to first element in class
	//for ( int i = 0; i < classIndexes.size(); i++ ){
	//	// slip first element
	//	for ( int j = 1; j < classIndexes[i].size(); j++ ){
	//		double matrix[16];
	//		for ( int k = 0; k < 16; k++ ){
	//			matrix[k] = this->alignmentBetweenReferences( classIndexes[i][0], classIndexes[i][j], 3 + k, 0 );
	//		}
	//		vtkTransform * t = vtkTransform :: New();
	//		t->Concatenate( matrix );
	//		vtkImageData * data = vtkImageData::New();
	//		this->references[ classIndexes[i][j] ].getImage( data, t );
	//		t->Delete();
	//		Array< double, 3 > A;
	//		nbfVTKInterface::vtkToBlitzReference( data, A );
	//		this->references[ classIndexes[i][j] ].setImage( A );
	//		data->Delete();
	//	}
	//}

	// look for volume that aligns best to all the rest (i.e. sum of distances is minimum)
	Array< Pixel, 2 > distances( this->alignmentBetweenReferences( Range::all(), Range::all(), 0, 0 ) );
	secondIndex j;
	TinyVector< int, 1 > bestOverall = minIndex( sum( distances, j ) );

	// align volumes within class to first element in class
	for ( int i = 0; i < distances.rows(); i++ ){
		if ( i != bestOverall[0] ){
			double matrix[16];
			for ( int k = 0; k < 16; k++ ){
				matrix[k] = this->alignmentBetweenReferences( bestOverall[0], i, 1 + k, 0 );
			}
			vtkTransform * t = vtkTransform :: New();
			t->Concatenate( matrix );
			vtkImageData * data = vtkImageData::New();
			this->references[ i ].getImage( data, t );
			t->Delete();
			Array< double, 3 > A;
			nbfVTKInterface::vtkToBlitzReference( data, A );
			this->references[ i ].setImage( A );
			data->Delete();
		}
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: compareAndMergeReferences()
{
	cout << "Comparing references..." << endl;

	// initialize new alignment matrix of references to volumes
	this->alignmentBetweenReferences.resize( this->classes.rows(), this->classes.rows(), 19 );
	this->alignmentBetweenReferences = -1;

	//vector< nbfWedgedAverageImage3D< Pixel > > temp;
	//temp.push_back( this->references[0] );
	//temp.push_back( this->references[1] );
	//this->references.clear();
	//this->references = temp;
	//this->alignmentBetweenReferences.resize( this->references.size(), this->references.size(), 19 );
	//this->alignmentBetweenReferences = -1;

	for ( int i = 0; i < this->references.size(); i++ ){
		cout << "Reference " << i << " has " << sum( this->references[i].weights > 0 ) << " non-zero elements." << endl;
	}

	// get best transformation
	this->metric->getDistances( this->references, this->alignmentBetweenReferences );

	Array< int, 1 > currentClasses( this->classes.rows() );
	currentClasses = 1;

#if 1

	// evaluate inter-reference distance with volume-wise measurements
	vector< nbfWedgedSubImage3D< Pixel > > list1, list2;
	vector< TinyVector< int, 2 > > positions;

	// for all references
	for ( int i = 0; i < this->alignmentBetweenReferences.rows(); i++ ){
		
		for ( int j = i + 1; j < this->alignmentBetweenReferences.rows(); j++ ){
		
			double matrix[16];
			for ( int m = 0; m < 16; m++ ){
				matrix[m] = this->alignmentBetweenReferences( i, j, 3 + m );
			}
			vtkMatrix4x4 * mat = vtkMatrix4x4::New();
			mat->DeepCopy( matrix );

			// traverse first reference components
			for ( int k = 0; k < this->input.size(); k++ ){

				if ( this->references[i].weights(k) > 0 ){
				
					list1.push_back( this->references[i].getVolumes()[k] );

					// for all second reference components
					for ( int l = 0; l < this->input.size(); l++ ){

						if ( this->references[j].weights(l) > 0 ){

							vtkTransform * concatenate = vtkTransform::New();
							this->references[j].getVolumes()[l].getTransform( concatenate );
							// concatenate->PostMultiply();
							concatenate->Concatenate( mat );
							this->input[l].setTransform( concatenate );
							concatenate->Delete();

							list2.push_back( this->input[l] );

							TinyVector< int, 2 > currentPosition( list1.size() - 1, list2.size() - 1 );
							positions.push_back( currentPosition );
						}
					}
				}
			}
			mat->Delete();
		}
	}

	vector< nbfWedgedImage3D< Pixel > * > plist1, plist2;

	for ( int j = 0; j < list1.size(); j++ ){
		plist1.push_back( &(list1[j]) );
	}

	for ( int j = 0; j < list2.size(); j++ ){
		plist2.push_back( &(list2[j]) );
	}

	Array< Pixel, 2 > allDistances( positions.size(), 19 );
	allDistances = -1;
	this->metric->getDistances( plist1, plist2, positions, allDistances, 2 ); // use executeFourierNew

	int count = 0;
	
	// for all references
	for ( int i = 0; i < this->alignmentBetweenReferences.rows(); i++ ){
		for ( int j = i + 1; j < this->alignmentBetweenReferences.rows(); j++ ){
			Array< Pixel, 2 > distances( this->alignmentToReferences.cols(), this->alignmentToReferences.cols() );
			distances = 0;
			// traverse first reference components
			for ( int k = 0; k < this->input.size(); k++ ){
				if ( this->references[i].weights(k) > 0 ){
					// for all second reference components
					for ( int l = 0; l < this->input.size(); l++ ){
						if ( this->references[j].weights(l) > 0 ){
							distances(k,l) = allDistances(count,0);
							count++;
						}
					}
				}
			}
			Array< Pixel, 1 > subDistances( this->input.size() );
			subDistances = 0;
			for ( int k = 0; k < this->input.size(); k++ ){
				subDistances( k ) = sum( this->references[j].weights * distances( k, Range::all() ) );
			}
			this->alignmentBetweenReferences( i, j, 0 ) = sum( this->references[i].weights * subDistances )	/ sum( this->references[i].weights ) / sum( this->references[j].weights );
			this->alignmentBetweenReferences( j, i, 0 ) = this->alignmentBetweenReferences( i, j, 0 );
			// if distance between references is smaller that sum of radii -> merge them into one
			if ( this->alignmentBetweenReferences( i, j, 0 ) < this->hierarchicalClasses ){
			//if ( this->alignmentBetweenReferences( i, j, 0 ) < this->referenceRadii(i) + this->referenceRadii(j) ){
			//if ( this->alignmentBetweenReferences( i, j, 0 ) < this->minDistanceBetweenReferences ){
				currentClasses(j) = 0;
			}
		}
	}

	cout << "Inter-reference distances = " << this->alignmentBetweenReferences( Range::all(), Range::all(), 0 ) << endl;

#else

	Array< Pixel, 2 > distances( this->classes.rows(), this->classes.rows() );
	distances = 0;

	// compute inter-reference distances as L1 norm of difference between weight vectors
	for ( int i = 0; i < this->references.size(); i++ ){
		for ( int j = i + 1; j < this->references.size(); j++ ){
			// if class still alive
			if ( currentClasses(j) > 0 ){
				// distances(i,j) = sum( abs( this->references[i].weights - this->references[j].weights ) );
				
				// normalize weight vectors
				Array< Pixel, 1 > w1( this->references[i].weights.shape() );
				w1 = this->references[i].weights / sqrt( sum( pow2(this->references[i].weights) ) );
				
				Array< Pixel, 1 > w2( this->references[j].weights.shape() );
				w2 = this->references[j].weights / sqrt( sum( pow2(this->references[j].weights) ) );
				
				// compute difference and magnitude of inner product		
				distances(i,j) = max( 0.0, abs( sum( w1 * w2 ) - 1.0 ) );
				
				// symmetrize
				distances(j,i) = distances(i,j);
				
				if ( distances(i,j) < this->minDistanceBetweenReferences ){
					currentClasses(j) = 0;
				}
			}
		}
	}

	cout << "Inter-reference distance = " << distances  << endl;

#endif

	Array< Pixel, 2 > classesAfterMerging( sum( currentClasses ), this->input.size() );
	vector< nbfWedgedAverageImage3D< Pixel > > referencesAfterMerging;

	int index = 0;
	for ( int i = 0; i < currentClasses.rows(); i++ ){
		if ( currentClasses(i) > 0 ){
			classesAfterMerging( index, Range::all() ) = this->classes( i, Range::all() );
			referencesAfterMerging.push_back( this->references[i] );
			index++;
		}
	}

	cout << "Classes after inter-reference merging = " << sum( currentClasses ) << " (of " << currentClasses.rows() << " before merging)" << endl;

	// assign new classes
	this->classes.resize( classesAfterMerging.shape() );
	this->classes = classesAfterMerging;

	// assign new references
	this->references.clear();
	this->references = referencesAfterMerging;
}



template< class Pixel >
void nbfLoopClustering< Pixel > :: saveState( stringstream & file, bool refinement )
{
	// save most current alignment of all raw volumes

	stringstream volumesFile;
	volumesFile << file.str() << "_volumes.txt";
	nbfWedgedSubImage3D< Pixel > :: write( volumesFile.str().c_str(), this->input );

	// save class memberships
	stringstream classesFile;
	classesFile << file.str() << "_classes.txt";
	vector< nbfWedgedSubImage3D< Pixel > > classVolumes;
	classVolumes = this->input;
	for ( int i = 0; i < classVolumes.size(); i++ ){
		classVolumes[i].setCutOffset(0);
	}
	for ( int i = 0; i < this->classes.rows(); i++ ){
		for ( int j = 0; j < this->classes.cols(); j++ ){
			if ( this->classes( i, j, 0 ) > 0 ){
				classVolumes[j].setCutOffset(i);
			}
		}
	}
	nbfWedgedSubImage3D< Pixel > :: write( classesFile.str().c_str(), classVolumes );

	// dump references into bin file
	stringstream objectData;
	objectData << this->references.size();
	for ( int i = 0; i < this->references.size(); i++ ){
		this->references[i].serialize( objectData );
	}
	stringstream processFile;
	processFile << file.str() << "_averages.bin";
	std :: ofstream dst( processFile.str().c_str(), std :: ofstream :: binary );
	dst << objectData.rdbuf();
	dst.close();

	// save references in human format
	vector< nbfWedgedSubImage3D< Pixel > > currentReferences;

	nbfWedgedSubImage3D< Pixel > nextVolume;
	TinyVector< int, 3 > size = this->input[0].getDimensions();
	nextVolume.setDimensions( size );
	TinyVector< Pixel, 3 > pos = size / 2.0;
	nextVolume.setPosition( pos );

	for ( int i = 0; i < this->references.size(); i++ ){

		stringstream fileName;

		if ( this->referencesIndexes.size() > 0 ){
			fileName << file.str() << "_average_" << this->hierarchicalClasses << "_" << this->referencesIndexes[i] << ".mrc";
		} else {
			if ( i < 10 ){
				fileName << file.str() << "_level_" << this->hierarchicalClasses << "_average_00" << i << ".mrc";
			} else if ( i < 100 ){
				fileName << file.str() << "_level_" << this->hierarchicalClasses << "_average_0" << i << ".mrc";
			} else {
				fileName << file.str() << "_level_" << this->hierarchicalClasses << "_average_" << i << ".mrc";
			}
		}
		
		// save image data (apply symmetry if neccesary)
		vtkImageData * averageVtk = vtkImageData::New();
		this->references[i].getImage( averageVtk );
		this->metric->imageFilter->symmetrize( averageVtk, this->symmetryFactor );
		Array< double, 3 > A;
		nbfVTKInterface :: vtkToBlitzReference( averageVtk, A );
		A *= -1;

		vtkImageCast * cast = vtkImageCast::New();
		cast->SetOutputScalarTypeToFloat();
		cast->SetInput( averageVtk );
		cast->Update();

		nbfMrcWriter mrcw;
		mrcw.setFileName( fileName.str().c_str() );
		if ( mrcw.write( cast->GetOutput() ) == false ){
			cerr << "ERROR: Failed to write " << fileName.str() << endl;
		}

		cast->Delete();
		averageVtk->Delete();

		nextVolume.setFileName( fileName.str().c_str() );
		currentReferences.push_back( nextVolume );
	}

	// save references file
	stringstream referencesFile;
	referencesFile << file.str() << "_averages.txt";
	nbfWedgedSubImage3D< Pixel > :: write( referencesFile.str().c_str(), currentReferences );

	if ( refinement == true ){

		// average selected references together
		if ( this->references.size() > 0 ){

			nbfWedgedAverageImage3D< Pixel > totalAverage;
			totalAverage = this->references[0];
			//Array< Pixel, 3 > newAlignments( this->references[0].weights.rows(), 17, blitz :: extrema :: max( 1, this->symmetryFactor ) );
			Array< Pixel, 3 > newAlignments( this->references[0].weights.rows(), 17, 1 );
			newAlignments = 0;
			for ( int i = 0; i < this->references[0].weights.rows(); i++ ){
				for ( int j = 0; j < this->references.size(); j++ ){
					if ( this->references[j].weights( i, 0 ) > 0 ){
						newAlignments( i, 0, Range::all() ) = this->references[j].weights( i, Range::all() );
						newAlignments( i, Range(1,toEnd), Range::all() ) = this->references[j].multipleAlignments( i, Range::all(), Range::all() );
					}
				}
			}
			totalAverage.setAlignments( newAlignments );
			this->metric->getImage( totalAverage );

			vtkImageData * averageVtk = vtkImageData::New();
			totalAverage.getImage( averageVtk );

			Array< double, 3 > A;
			nbfVTKInterface :: vtkToBlitzReference( averageVtk, A );
			A *= -1;

			vtkImageCast * cast = vtkImageCast::New();
			cast->SetOutputScalarTypeToFloat();
			cast->SetInput( averageVtk );
			cast->Update();

			stringstream fileName;
			fileName << file.str() << "_global_average.mrc";

			nbfMrcWriter mrcw;
			mrcw.setFileName( fileName.str().c_str() );
			if ( mrcw.write( cast->GetOutput() ) == false ){
				cerr << "ERROR: Failed to write " << fileName.str() << endl;
			}

			// save masked global average
			if ( false ){
				stringstream fileNameMasked;
				fileNameMasked << file.str() << "_global_average_masked.mrc";
				this->metric->imageFilter->execute( averageVtk );

				nbfVTKInterface :: vtkToBlitzReference( averageVtk, A );
				A *= -1;

				cast->Modified();
				cast->Update();
				mrcw.setFileName( fileNameMasked.str().c_str() );
				if ( mrcw.write( cast->GetOutput() ) == false ){
					cerr << "ERROR: Failed to write " << fileNameMasked.str() << endl;
				}
			}

			cast->Delete();
			averageVtk->Delete();

			// split reference in two for FSC calculations

			nbfWedgedAverageImage3D< Pixel > firstHalf, secondHalf;

			firstHalf = totalAverage;
			secondHalf = totalAverage;

			firstHalf.weights = 0;
			secondHalf.weights = 0;

			int counter = 0;

			// traverse all references
			for ( int j = 0; j < this->input.size(); j++ ){
				if ( totalAverage.weights( j, 0 ) > 0 ){
					if ( counter % 2 == 0 ){
						firstHalf.weights( j, Range::all() ) = totalAverage.weights( j, Range::all() );
						firstHalf.multipleAlignments( j, Range::all(), Range::all() ) = totalAverage.multipleAlignments( j, Range::all(), Range::all() );
						secondHalf.weights( j, Range::all() ) = 0;
					} else {
						firstHalf.weights( j, Range::all() ) = 0;
						secondHalf.weights( j, Range::all() ) = totalAverage.weights( j, Range::all() );
						secondHalf.multipleAlignments( j, Range::all(), Range::all() ) = totalAverage.multipleAlignments( j, Range::all(), Range::all() );
					}
					counter++;
				}
			}

			// force state to out of date
			firstHalf.updateState();
			secondHalf.updateState();

			this->metric->getImage( firstHalf );
			this->metric->getImage( secondHalf );

			for ( int k = 1; k < 3; k++ ){

				vtkImageData * averageVtk = vtkImageData::New();

				if ( k == 1 ){
					firstHalf.getImage( averageVtk );
				} else {
					secondHalf.getImage( averageVtk );
				}

				Array< double, 3 > A;
				nbfVTKInterface :: vtkToBlitzReference( averageVtk, A );
				A *= -1;

				vtkImageCast * cast = vtkImageCast::New();
				cast->SetOutputScalarTypeToFloat();
				cast->SetInput( averageVtk );
				cast->Update();

				stringstream fileName;

				fileName << file.str() << "_fsc_" << k << ".mrc";

				nbfMrcWriter mrcw;
				mrcw.setFileName( fileName.str().c_str() );
				if ( mrcw.write( cast->GetOutput() ) == false ){
					cerr << "ERROR: Failed to write " << fileName.str() << endl;
				}
				cast->Delete();
				averageVtk->Delete();
			}
		}
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: loadState( stringstream & file )
{
	nbfMatlabReader r;
	if ( this->running_mode == NBF_LOOP_CLUSTERING_CLASS ){

		// load alignments to references
		stringstream alignmentsFile;
		alignmentsFile << file.str() << "_alignmentToReferences.matlab";
		r.setFileName( alignmentsFile.str().c_str() );
		r.read( this->alignmentToReferences );

		if ( this->alignmentToReferences.size() > 0 && this->alignmentToReferences.cols() != this->input.size() ){
			cerr << "ERROR - File " << alignmentsFile.str().c_str() << " size is inconsistent with volume file size." << endl;
		}
	} else {

		// store current geometry parameters
		TinyVector< int, 3 > currentVolumeDimensions = this->input[0].getDimensions();
		Pixel currentCutOffset = this->input[0].getCutOffset();

		// retrieve current volume alignments from latest run
		stringstream volumesFile;
		volumesFile << file.str() << "_volumes.txt";
		this->input.clear();
		nbfWedgedSubImage3D< Pixel > :: read( volumesFile.str().c_str(), this->input );

		//// extract current alignments and apply specified size and offset parameters
		//Array< Pixel, 3 > inputAlignments( this->input.size(), 17, 1 );
		//for ( int i = 0; i < this->input.size(); i++ ){
		//	vtkTransform * t = vtkTransform :: New();
		//	this->input[i].getTransform(t);
		//	double matrix[16];
		//	vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
		//	t->Delete();
		//	for ( int k = 0; k < 16; k++ ){
		//		inputAlignments( i, 1 + k, 0 ) = matrix[k];
		//	}
		//	// set new geometry
		//	this->input[i].setDimensions( currentVolumeDimensions );
		//	this->input[i].setCutOffset( currentCutOffset );
		//}

		// retrieve references from file
		stringstream referencesBinaryFile;
		referencesBinaryFile << file.str() << "_averages.bin";
		std :: ifstream src( referencesBinaryFile.str().c_str(), std :: ofstream :: binary );
		stringstream inputString1;
		src >> inputString1.rdbuf();
		src.close();
		int numberOfReferences;
		inputString1 >> numberOfReferences;
		this->references.clear();
		for ( int i = 0; i < numberOfReferences; i++ ){
            nbfWedgedAverageImage3D< Pixel > vol1;
			vol1.unserialize( inputString1 );
			this->references.push_back( vol1 );
		}
	
		if ( ( numberOfReferences > 0 ) && ( this->references[0].weights.rows() == this->input.size() ) ){
			// everything is OK
		} else {
			cerr << "ERROR: volume list and reference memberships have different sizes: " << this->references[0].weights.rows() << " != " << this->input.size() << ". In " << __FILE__ << " : " << __LINE__ << endl;
			return;
		}

		// retrieve user class selection
		stringstream referencesFile;
		referencesFile << file.str() << "_averages.txt";
		//vector< nbfWedgedSubImage3D< Pixel > > currentReferences;
		nbfWedgedSubImage3D< Pixel > :: read( referencesFile.str().c_str(), this->precomputedReferences );

			//// reset memberships to ground truth
			//cout << "WARNING: OVERRIDING CLASS MEMBERSHIPS WITH GROUND TRUTH" << endl;
			//nbfMrcReader w;
			//w.setFileName( this->precomputedReferences[0].getFileName().c_str() );
			//vtkImageData * data = vtkImageData :: New();
			//w.read( data );
			//Array< double, 3 > Avg;
			//nbfVTKInterface::vtkToBlitz( data, Avg );
			//Avg *= -1;
			//this->references[0].setAverageImage( Avg );
			//this->references[0].weights( Range(fromStart,4843), Range :: all() ) = 1;
			//this->references[0].weights( Range(4844,toEnd), Range :: all() ) = 0;

			//w.setFileName( this->precomputedReferences[1].getFileName().c_str() );
			//w.read( data );
			//nbfVTKInterface::vtkToBlitz( data, Avg );
			//Avg *= -1;
			//this->references[1].setAverageImage( Avg );
			//this->references[1].weights( Range(fromStart,4843), Range :: all() ) = 0;
			//this->references[1].weights( Range(4844,toEnd), Range :: all() ) = 1;
			//data->Delete();
			//// done reseting memberships to ground truth

		vector< int > indexes( this->precomputedReferences.size() );
		int discarded = 0;
		for ( int i = 0; i < this->precomputedReferences.size(); i++ ){
			indexes[i] = this->precomputedReferences[i].getCutOffset();
			if ( abs( indexes[i] ) == 0 ){
				discarded++;
			}
			this->precomputedReferences[i].setDimensions( currentVolumeDimensions );
			this->precomputedReferences[i].setCutOffset( 0 );

			if ( this->references[i].wedgeImage.size() > 0 ){
				this->precomputedReferences[i].wedge.setImage( this->references[i].wedgeImage );
			}
			if ( this->references[i].sphericalWedgeImage.size() > 0 ){
				this->precomputedReferences[i].wedge.setSphericalImage( this->references[i].sphericalWedgeImage );
				//nbfMatlabWriter w;
				//w.setFileName("p.matlab");
				//w.write(this->precomputedReferences[i].wedge.sphericalWedge);
			}
		}

		// store new combination of classes
		this->classesSelected.clear();

		// if no selection was made, keep all classes and use biggest as reference
		if ( discarded == this->precomputedReferences.size() ){

			// find class with most volumes
			Array< int, 1 > numberOfVolumesPerClass( this->precomputedReferences.size() );
			for ( int i = 0; i < this->precomputedReferences.size(); i++ ){
				numberOfVolumesPerClass(i) = sum( this->references[i].weights );
			}
			TinyVector< int, 1 > useAsReference = maxIndex( numberOfVolumesPerClass );

			this->classesSelected.resize( this->precomputedReferences.size() );
			int count = 1;
			for ( int i = 0; i < this->precomputedReferences.size(); i++ ){
				if ( i == useAsReference[0] ){
					this->classesSelected[0].push_back(i);
				} else {
					this->classesSelected[count].push_back(i);
					count++;
				}
			}

		} else {
			// get rid of un-selected classes
			int reducedNumberOfClasses = 0;
			vector< int > :: iterator iter = indexes.begin();
			typename vector< nbfWedgedAverageImage3D< Pixel > > :: iterator iterR = this->references.begin();
			typename vector< nbfWedgedSubImage3D< Pixel > > :: iterator iterP = this->precomputedReferences.begin();
			while ( iter != indexes.end() ){
				if ( abs( *iter ) == 0 ){
					this->references.erase( iterR );
					this->precomputedReferences.erase( iterP );
					indexes.erase( iter );
				} else {
					if ( abs(*iter) > reducedNumberOfClasses ){
						reducedNumberOfClasses = abs(*iter);
					}
					iter++;
					iterR++;
					iterP++;
				}
			}

			this->classesSelected.resize( reducedNumberOfClasses );
			for ( int i = 0; i < indexes.size(); i++ ){
				int currentIndex = indexes[i];
				if ( currentIndex > 0 ){
					this->classesSelected[ currentIndex - 1 ].push_back(i);
				} else {
					this->classesSelected[ abs(currentIndex) - 1 ].insert( this->classesSelected[ abs(currentIndex) - 1 ].begin(), i );
				}
			}
		}

		//// swap reference to first slot
		//swap( this->references[0], this->references[this->classesSelected[0][0]] );

		// extract class information
		this->classes.resize( this->references.size(), this->references[0].weights.rows(), this->references[0].weights.cols() );
		for ( int i = 0; i < this->references.size(); i++ ){
			this->classes( i, Range :: all(), Range :: all() ) = this->references[i].weights;
		}
	}
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: symmetrize( nbfWedgedAverageImage3D< Pixel > & reference, Pixel fold, bool align )
{
	//// QUICK AND DIRTY

	//if ( fold > 1 ){
	//	Array< double, 3 > result;

	//	vtkImageData * data = vtkImageData :: New();

	//	// assume symmetry axis is Z direction
	//	for ( int i = 0; i < fold; i++ ){
	//		vtkTransform * t = vtkTransform :: New();
	//		t->RotateZ( i * 360.0 / fold );
	//		reference.getImage( data, t );
	//		t->Delete();
	//		Array< double, 3 > current;
	//		nbfVTKInterface :: vtkToBlitzReference( data, current );
	//		if ( result.size() == 0 ){
	//			result.resize( current.shape() );
	//			result = 0;
	//		}
	//		result += current;
	//	}

	//	data->Delete();

	//	result /= fold;
	//	reference.setAverageImage( result );
	//}

	// WELL DONE
	
	vtkImageData * data = vtkImageData :: New();
	reference.getImage( data );	
	Array< double, 3 > newAverage;
	nbfVTKInterface::vtkToBlitz( data, newAverage );
	newAverage = 0;

	if ( fold > 1 ){
		Array< Pixel, 3 > newAlignments( reference.weights.rows(), 17, fold );

		for ( int f = 0; f < fold; f++ ){

			newAlignments( Range::all(), 0, f ) = reference.weights( Range::all(), 0 );

			vtkTransform * t = vtkTransform :: New();
			t->RotateZ( f * 360.0 / fold );

			// align current fold to initial reference
			if ( ( f > 0 ) && ( align == true ) ){
				nbfFourierImageMetric< Pixel, 3 > fMetric( this->metric->imageFilter, this->metric->fourierFilter );
				fMetric.setRotationSearchRestriction( this->metric->getRotationSearchRestriction() );
				fMetric.setTranslationSearchRestriction( this->metric->getTranslationSearchRestriction() );
				fMetric.setToComputeOverlapNormalizedDistances( this->metric->getToComputeOverlapNormalizedDistances() );
				fMetric.setToUseMutualCorrelation( this->metric->getToUseMutualCorrelation() );
				//currentMetric->imageFilter->maskOn( this->imageFilter->maskFile );
				fMetric.setInput1( &reference );
				fMetric.setInput2( &reference );
				fMetric.executeFourierNewHalf( t );
			}

			for ( int i = 0; i < reference.weights.rows(); i++ ){

				// retrieve current transformation
				vtkMatrix4x4 * mat1 = vtkMatrix4x4 :: New();

				double mat[16];
				for ( int m = 0; m < 16; m++ ){
					mat[m] = reference.multipleAlignments( i, m, 0 );
				}
				mat1->DeepCopy( mat );

				vtkMatrix4x4 * mat3 = vtkMatrix4x4 :: New();
				vtkMatrix4x4 :: Multiply4x4( mat1, t->GetMatrix(), mat3 );

				vtkMatrix4x4::DeepCopy( mat, mat3 );
				for ( int m = 0; m < 16; m++ ){
					newAlignments( i, 1 + m, f ) = mat[m];
				}
				mat3->Delete();
				mat1->Delete();
			}

			reference.getImage( data, t );	
			Array< Pixel, 3 > A;
			nbfVTKInterface::vtkToBlitz( data, A );
			newAverage += A;

			t->Delete();
		}

		newAverage /= fold;
		reference.setAlignments( newAlignments );
		if ( this->metric->getMissingWedgeCompensation() == false ){
			reference.setAverageImage( newAverage );
		} else {
			cerr << "WARNING - Symmetry imposition does not cause the missing-wedge to be re-computed." << endl;
			reference.setAverageImage( newAverage );
			reference.resetWedge();
			// this->metric->getImage( reference );
		}
	}

	if ( ( fold > 1 ) && ( align == true ) ){
		Array< double, 3 > result;

		// assume symmetry axis is Z direction
		for ( int i = 0; i < fold; i++ ){
			vtkTransform * t = vtkTransform :: New();
			t->RotateZ( i * 360.0 / fold );
			reference.getImage( data, t );
			t->Delete();
			Array< double, 3 > current;
			nbfVTKInterface :: vtkToBlitzReference( data, current );
			if ( result.size() == 0 ){
				result.resize( current.shape() );
				result = 0;
			}
			result += current;
		}

		result /= fold;
		reference.setAverageImage( result );
	}

	data->Delete();
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: doPCADecomposition( Array< Pixel, 3 > & fullR, int numFactors )
{
	Array< Pixel, 2 > R( fullR( Range :: all(), Range :: all(), 0 ) );

	firstIndex i;
	secondIndex j;

	// extract the mean from each vectorized image
	Array< Pixel, 1 > M( R.cols() );
	M = mean( R(j,i), j );
	for ( int in = 0; in < R.rows(); in++ ){
		R( in, Range::all() ) -= M;
	}

	// store covariance matrix A'*A:
	Array< Pixel, 2 > L( R.rows(), R.rows() );
	
	cout <<  "Computing covariance matrix..." << endl;

	// Fast computation of covariance matrix. 
	// This replaces the more expensive blitz operation: L = sum( R(i,k) * P(k,j), k );
	for ( int in = 0; in < R.rows(); in++ ){
		Array< Pixel, 1 > cu( R( in, Range::all() ) );
		L(in,in) = sum( cu * cu );
		for ( int jn = in + 1; jn < R.rows(); jn++ ){
			L(in,jn) = sum( cu * R( jn, Range::all() ) );
			L(jn,in) = L(in,jn);
		}
	}

	// Covariance matrix: provide C pointer interface for processing by VTK's Jacobi routine
	Pixel * L_data = L.data();
	Pixel ** cL = new Pixel *[ L.rows() ];
	for ( int i = 0; i < L.rows(); i++) {
		cL[i] = &L_data[i*L.cols()];
	}

	// store eigenvalues
	Array< Pixel, 1 > vals( L.rows() );

	// store eigenvectors and provide C pointer interface for processing by VTK's Jacobi routine
	Array< Pixel, 2 > vecs( L.shape() );
	Pixel * vecs_data = vecs.data();
	Pixel ** cvecs = new Pixel *[ L.rows() ];
	for ( int i = 0; i < L.rows(); i++) {
		cvecs[i] = &vecs_data[i*L.cols()];
	}

	cout <<  "Computing Jacobi diagonalization..." << endl;

	// get the eigenvalues and eigenvectors	
	vtkMath :: JacobiN( cL, L.rows(), vals.data(), cvecs );

	// convert the eigenvectors of A'*A into eigenvectors of A*A': Vectors = A * Vectors
	Array< Pixel, 2 > nvecs( R.cols(), numFactors );
	for ( int i = 0; i < R.cols(); i++ ){
		Array< Pixel, 1 > fR( R( Range :: all(), i ) );
		for ( int j = 0; j < numFactors; j++ ){
			nvecs(i,j) = sum( fR * vecs( Range::all(), j ) );
		}
	}

	// normalize eigenvalues
	vals /= ( R.rows() - 1.0 );

	// normalize eigenvectors to unit length
	for ( int in = 0; in < nvecs.cols(); in++ ){
		nvecs( Range :: all(), in ) /= sqrt( sum( pow2( nvecs( Range :: all(), in ) ) ) );
	}

	// compute coordinates in lower dimensional space
	vecs.resize( R.rows(), numFactors );
	for ( int i = 0; i < R.rows(); i++ ){
		Array< Pixel, 1 > cR( R( i, Range :: all() ) );
		for ( int j = 0; j < nvecs.cols(); j++ ){
			vecs(i,j) = sum( cR * nvecs( Range :: all(), j ) );
		}
	}

	// copy to output
	fullR.resize( vecs.rows(), vecs.cols(), 1 );
	fullR( Range :: all(), Range :: all(), 0 ) = vecs;

	// cleanup
	delete [] cL;
	delete [] cvecs;
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: doSVDDecomposition( Array< Pixel, 3 > & fullR, int numFactors )
{
	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	TinyVector< int, 3 > fullRsize = fullR.shape();

	// CENTER DATASET

	// compute global average
	Array< Pixel, 2 > sR( fullR( Range :: all(), Range :: all(), 0 ) );
	Array< Pixel, 1 > globalAverage( sR.cols() );
	firstIndex ix;
	secondIndex jx;
	globalAverage = mean( sR(jx,ix), jx );

	Array< Pixel, 1 > distancesToGlobalAverage( sR.rows() );
	vector< Pixel > sortedDistancesToGlobalAverage( distancesToGlobalAverage.size() );

	// subtract global average and compute distance from center
	for ( int i = 0; i < sR.rows(); i++ ){
		sR( i, Range :: all() ) -= globalAverage;
		distancesToGlobalAverage(i) = sqrt( sum( pow2( sR( i, Range :: all() ) ) ) );
		sortedDistancesToGlobalAverage[i] = distancesToGlobalAverage(i);
	}

	cout << "Point cloud has been centered to the mean of the dataset." << endl;

	// sort by ascending distance from center
	sort( sortedDistancesToGlobalAverage.begin(), sortedDistancesToGlobalAverage.end() );
	
	// determine cutoff distance
	Pixel thresholdDistance = sortedDistancesToGlobalAverage[ floor( ( sortedDistancesToGlobalAverage.size() - 1 ) * this->distanceTopCutoffPre ) ];

	// find out what volumes are left
	vector< int > survivingVolumes;
	for ( int i = 0; i < distancesToGlobalAverage.size(); i++ ){
		if ( distancesToGlobalAverage(i) <= thresholdDistance ){
			survivingVolumes.push_back(i);
		}
	}

	if ( survivingVolumes.size() < distancesToGlobalAverage.size() ) {
		cout << "Eliminating " << distancesToGlobalAverage.size() - survivingVolumes.size() << " outlier volumes corresponding to " << 100*(1 - this->distanceTopCutoffPre) << "% of total." << endl;
	}

	fullR.transposeSelf( secondDim, firstDim, thirdDim );
	//Array< double, 2 > A( fullR.rows(), fullR.cols() );
	//A = cast< double >( fullR( Range :: all(), Range :: all(), 0 ) );

	Array< double, 2 > A( fullR.rows(), survivingVolumes.size() );
	for ( int i = 0; i < survivingVolumes.size(); i++ ){
		A( Range :: all(), i ) = cast< double>( fullR( Range :: all(), survivingVolumes[i], 0 ) );
	}

	//Array< double, 2 > A( fullR.rows(), fullR.cols() );
	//A = cast< double >( fullR );

	fullR.free();

	//// normalize to zero mean and unit variance
	//for ( int i = 0; i < A.cols(); i++ ){
	//	Array< double, 1 > lA( A( Range :: all(), i ) );
	//	double m = mean( lA );
	//	double var = sqrt( mean( pow2( lA - m ) ) );
	//	lA = ( lA - m ) / var;
	//	// lA = lA - m;
	//}

	//w.write(A);

	//cout << "doing GSVD " << endl;

	firstIndex i;
	secondIndex j;

	// left matrix metric
	Array< float, 1 > M( A.rows() );

	// right matrix metric
	Array< float, 1 > N( A.cols() );

	// double-stochastic iteration
	for ( int iter = 0; iter < 10; iter++ ){
	
		// left matrix metric
		M = 1.0 / sqrt( sum( pow2( A(i,j) ), j ) );

		// right matrix metric
		N = 1.0 / sqrt( sum( pow2( A(j,i) ), j ) );

		for ( int i = 0; i < A.cols(); i++ ){
			A( Range :: all(), i ) *= (M);
		}

		//w.write(A);

		for ( int i = 0; i < A.rows(); i++ ){
			A( i, Range :: all() ) *= (N);
		}

		//w.write(A);
	}

	// turn off SVDLIBC verbose

	// convert A to svdlibc
	double * A_data = A.data();
	double ** cA = new double *[ A.rows() ];
	for ( int i = 0; i < A.rows(); i++) {
		cA[i] = &A_data[i*A.cols()];
	}
	
	dmat Asvd;
	Asvd.rows = A.rows();
	Asvd.cols = A.cols();
	Asvd.value = cA;

	// convert dense to sparse representation
	DMat Ad = &Asvd;
	SMat As = NULL;
	As = svdConvertDtoS(Ad);
	
	//cout << "running SVD " << endl;

	// run svd
	SVDRec R = svdLAS2A( As, numFactors );

	Array< double, 1 > S( R->S, R->d );

	Array< double, 2 > V( R->Vt->value[0], shape( R->Vt->rows, R->Vt->cols ) );

	cout << "Eigenvectors = \n" << S << endl;

	// transform back
	for ( int i = 0; i < V.rows(); i++ ){
		V( i, Range::all() ) /= sqrt(N);
	}

	// multiply by eigenvalues to compute new representation
	for ( int i = 0; i < V.cols(); i++ ){
		V( Range( fromStart, S.ubound(firstDim) ), i ) *= S;
		if ( S.rows() < V.rows() ){
			V( Range( S.rows(), toEnd ), i ) = 0;
		}
	}
	//w.write(V);

	V.transposeSelf(secondDim,firstDim);

	//Array< double, 2 > U( R->Ut->value[0], shape( R->Ut->rows, R->Ut->cols ) );
	//for ( int i = 0; i < U.rows(); i++ ){
	//	U( i, Range::all() ) /= sqrt(M);
	//}
	//w.write(U);

	svdFreeSMat(As);

	// copy to output
	//fullR.resize( V.rows(), V.cols(), 1 );
	//fullR( Range::all(), Range::all(), 0 ) = cast< Pixel >(V);

	//fullR.resize( fullRsize[0], V.cols(), 1 );
	//fullR = 0;
	//for ( int i = 0; i < survivingVolumes.size(); i++ ){
	//	fullR( survivingVolumes[i], Range :: all(), 0 ) = cast< Pixel >( V( i, Range :: all() ) );
	//}

	fullR.resize( fullRsize[0], S.rows(), 1 );
	fullR = 0;
	for ( int i = 0; i < survivingVolumes.size(); i++ ){
		fullR( survivingVolumes[i], Range :: all(), 0 ) = cast< Pixel >( V( i, Range( fromStart, S.rows() - 1 ) ) );
	}

	svdFreeSVDRec(R);
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: doSpectralSVDDecomposition( Array< Pixel, 3 > & fullR, int numFactors )
{
	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");

	// fullR.transposeSelf( secondDim, firstDim, thirdDim );
	Array< double, 2 > A( fullR.rows(), fullR.cols() );
	A = cast< double >( fullR( Range :: all(), Range :: all(), 0 ) );

	//Array< double, 2 > A( fullR.rows(), fullR.cols() );
	//A = cast< double >( fullR );

	fullR.free();

	//// normalize to zero mean and unit variance
	//for ( int i = 0; i < A.cols(); i++ ){
	//	Array< double, 1 > lA( A( Range :: all(), i ) );
	//	double m = mean( lA );
	//	double var = sqrt( mean( pow2( lA - m ) ) );
	//	lA = ( lA - m ) / var;
	//	// lA = lA - m;
	//}

	//w.write(A);

	//cout << "doing GSVD " << endl;

	firstIndex i;
	secondIndex j;

	// ESTIMATE SIGMA FROM DATA (average distance to closest neighbor)
	// set diagonal to infty
	for ( int k = 0; k < A.rows(); k++ ){
		A(k,k) = numeric_limits< double > :: max();
	}
	// search for closest neighbor
	Array< double, 1 > closestNeighbor( A.rows() );
	closestNeighbor = min( A, j );
	Pixel beta = mean( closestNeighbor );
	
	// set diagonal back to zero
	for ( int k = 0; k < A.rows(); k++ ){
		A(k,k) = 0;
	}

	// Pixel beta = 1.0 / 3.0;
	A = exp( - A / 2.0 / beta );

	//// set diagonal to zero
	//for ( int k = 0; k < A.rows(); k++ ){
	//	A(k,k) = 0;
	//}

	// left matrix metric
	Array< float, 1 > M( A.rows() );
	M = 1.0 / sum( A(i,j), j );
	

	// L = D^-1/2 * A * D^-1/2
	for ( int i = 0; i < A.cols(); i++ ){
		A( Range :: all(), i ) *= sqrt(M);
	}

	for ( int i = 0; i < A.rows(); i++ ){
		A( i, Range :: all() ) *= sqrt(M);
	}

	// turn off SVDLIBC verbose
	SVDVerbosity = 2;

	// convert A to svdlibc
	double * A_data = A.data();
	double ** cA = new double *[ A.rows() ];
	for ( int i = 0; i < A.rows(); i++) {
		cA[i] = &A_data[i*A.cols()];
	}
	
	dmat Asvd;
	Asvd.rows = A.rows();
	Asvd.cols = A.cols();
	Asvd.value = cA;

	// convert dense to sparse representation
	DMat Ad = &Asvd;
	SMat As = NULL;
	As = svdConvertDtoS(Ad);
	
	//cout << "running SVD " << endl;

	// run svd
	SVDRec R = svdLAS2A( As, numFactors );

	Array< double, 1 > S( R->S, R->d );

	Array< double, 2 > V( R->Vt->value[0], shape( R->Vt->rows, R->Vt->cols ) );

	cout << "Eigenvectors = \n" << S << endl;

	//// multiply by eigenvalues to compute new representation
	//for ( int i = 0; i < V.cols(); i++ ){
	//	V( Range( fromStart, S.ubound(firstDim) ), i ) *= S;
	//	if ( S.rows() < V.rows() ){
	//		V( Range( S.rows(), toEnd ), i ) = 0;
	//	}
	//}

	// re-normalize
	for ( int k = 0; k < V.cols(); k++ ){
		V( Range::all(), k ) /= sqrt( sum( pow2( V( Range :: all(), k ) ) ) );
	}
	//w.write(V);

	V.transposeSelf(secondDim,firstDim);

	//Array< double, 2 > U( R->Ut->value[0], shape( R->Ut->rows, R->Ut->cols ) );
	//for ( int i = 0; i < U.rows(); i++ ){
	//	U( i, Range::all() ) /= sqrt(M);
	//}
	//w.write(U);

	svdFreeSMat(As);

	// copy to output
	fullR.resize( V.rows(), V.cols(), 1 );
	fullR( Range::all(), Range::all(), 0 ) = cast< Pixel >(V);
	//fullR.resize( V.rows(), V.cols() );
	//fullR = cast< Pixel >(V);

	svdFreeSVDRec(R);
}

template< class Pixel >
void nbfLoopClustering< Pixel > :: imputeInReciprocalSpace( Array< Pixel, 3 > & R, Array< Pixel, 2 > & D )
{
	// IMPUTE WITH K-NEAREST NEIGHBORS

	// create separate views for each data component
	Array< Pixel, 2 > rall( R( Range :: all(), Range :: all(), 0 ) );
	Array< Pixel, 2 > call( R( Range :: all(), Range :: all(), 1 ) );
	Array< Pixel, 2 > wall( R( Range :: all(), Range :: all(), 2 ) );
	
	// for each volume
	for ( int i = 0; i < D.rows(); i++ ){

		// create separate views for current volume components
		Array< Pixel, 1 > r( R( i, Range :: all(), 0 ) );
		Array< Pixel, 1 > c( R( i, Range :: all(), 1 ) );
		Array< Pixel, 1 > w( R( i, Range :: all(), 2 ) );
		
		firstIndex ix;
		secondIndex jx;

		Array< Pixel, 1 > d( D( i, Range :: all() ) );

		// compute exponential weights
		Array< Pixel, 1 > e( exp( - d ) );

		// normalize weights
		Array< Pixel, 1 > sumw( sum( wall( jx, ix ) * e(jx), jx ) );

		// impute new values where possible
		r = where( ( w == 0 ) && ( sumw > 0 ), sum( rall( jx, ix ) * wall( jx, ix ) * e(jx), jx ) / sumw, r );
		c = where( ( w == 0 ) && ( sumw > 0 ), sum( call( jx, ix ) * wall( jx, ix ) * e(jx), jx ) / sumw, c );
	}
}

//template< class Pixel >
//void nbfLoopClustering< Pixel > :: kmeans( Array< Pixel, 2 > & R, Array< Pixel, 3 > & classes )
//{
//	// eliminate global outliers
//
//	// compute global average
//	Array< Pixel, 2 > globalAverage( R.cols(), 2 );
//	firstIndex i;
//	secondIndex j;
//	Array< Pixel, 1 > wedges( R( Range :: all(), Range :: all(), 2 ) );
//	Array< Pixel, 1 > accumulatedWedges( R.cols() );
//	accumulatedWedges = sum( wedges(j,i), j );
//	Array< Pixel, 2 > real( R( Range :: all(), Range :: all(), 0 );
//	globalAverage( Range :: all(), 0 ) = sum( real(j,i), j ) / where( accumulatedWedges > 0, accumulatedWedges, 1 );
//	Array< Pixel, 2 > imag( R( Range :: all(), Range :: all(), 1 );
//	globalAverage( Range :: all(), 1 ) = sum( imag(j,i), j ) / where( accumulatedWedges > 0, accumulatedWedges, 1 );
//	
//	// compute cutoff distance to global average
//	Array< Pixel, 2 > distances( R.rows(), 1 );
//	this->metric->getDistances( representationFile, globalAverage, distances );
//
//	vector< Pixel > sortedDistances;
//	for ( int i = 0; i < distances.rows(); i++ ){
//		sortedDistances.push_back( distances(i) );
//	}
//	sort( sortedDistances.begin(), sortedDistances.end() );
//	Pixel maxDistance = sortedDistances[ floor( ( sortedDistances.size() - 1 ) * this->distanceTopCutoffPre ) ];
//
//	// number of volumes submitted to classification
//	int effectiveVolumesToClassify = sum( where( distances < maxDistance, 1, 0 ) );
//
//	cout << "Eliminating " << distances.rows() - effectiveVolumesToClassify << " outlier volumes corresponding to " << 1 - this->distanceTopCutoff << " percentage of total." << endl;
//
//	// extract points and store index information
//	Array< Pixel, 2 > A( effectiveVolumesToClassify, R.cols(), R.depth() );
//	Array< int, 1 > indexes( A.rows() );
//	int count = 0;
//	for ( int i = 0; i < R.rows(); i++ ){
//		if ( distances(i) < maxDistance ){
//			A( count, Range :: all(), Range :: all() ) = R( i, Range :: all(), Range :: all() );
//			indexes(count++) = i;
//		}
//	}
//	///////////////////////
//
//	int numberOfPasses = 10;
//	int numberOfIterations = 10;
//
//	Pixel bestScore = numeric_limits< Pixel > :: max();
//
//	for ( int pass = 0; i < numberOfPasses; i++ ){
//		
//		// select random centroids from given points
//		Array< Pixel, 3 > centroids( numberOfClasses, A.cols(), 3 );
//		for ( int h = 0; h < numberOfClasses; h++ ){
//			ranlib :: Uniform<float> noise;
//			noise.seed((unsigned int)time(0));
//			int index = floor( noise.random() * A.ubound(firstDim) );
//			centroids( h, Range :: all(), Range :: all() ) = A( index, Range :: all(), Range :: all() );
//		}
//		
//		for ( int j = 0; j < numerOfIterations; j++ ){
//		
//			// compute distances to centroids
//			Array< Pixel, 2 > distancesToCentroids( A.rows(), centroids.rows() );
//			this->metric->getDistances( newRepresentationFile, centroids, distancesToCentroids );
//
//			// assign to closest centroid
//			secondIndex jindex;
//			Array< Pixel, 1 > currentClasses( classes.rows() );
//			currentClasses = minIndex( distancesToCentroids, jindex );
//
//			// recompute centroids
//			for ( int k = 0; k , numberOfClasses; k++ ){
//				Array< Pixel, 1 > accumulatedWedges( A.cols() );
//				accumulatedWedges = sum( wedges(j,i), j );
//				Array< Pixel, 2 > real( R( Range :: all(), Range :: all(), 0 );
//				centroids( k, Range :: all(), 0 ) = sum( where( currentClasses == k, 1, 0 ) * real(j,i), j ) / where( accumulatedWedges > 0, accumulatedWedges, 1 );
//				Array< Pixel, 2 > imag( R( Range :: all(), Range :: all(), 1 );
//				centroids( k, Range :: all(), 1 ) = sum( where( currentClasses == k, 1, 0 ) * imag(j,i), j ) / where( accumulatedWedges > 0, accumulatedWedges, 1 );
//			}
//		}
//
//		// compute current kmeans score
//		Pixel currentScore = 0;
//		for (){
//			currentScore += ;
//		}
//
//		if ( currentScore < bestScore ){
//			this->classes = currentClasses;
//			bestScore = currentScore;
//		}
//	}
//}


//template< class Pixel >
//void nbfLoopClustering< Pixel > :: postProcessor( Array< Pixel, 2 > & R, Array< Pixel, 3 > & classes )
//{
//	// Post-Processor
//	if ( classes.size() > 0 ){
//
//		vector< vector< int > > finalClasses;
//		for ( int i = 0; i < classes.rows(); i++ ){
//			vector< int > elements;
//			for ( int j = 0; j < classes.cols(); j++ ){
//				if ( classes( i, j ) > 0 ){
//					elements.push_back(j);
//				}
//			}
//			finalClasses.push_back( elements );
//		}
//
//		Array< Pixel, 1 > centroid( R.cols() );
//		Array< Pixel, 2 > classCoords( R.shape() );
//
//		// evaluate current score
//		Array< Pixel, 1 > scores( classes.size() );
//		for ( int i = 0; i < classes.size(); i++ ){
//			// retrive current class memberships
//			Array< Pixel, 1 > currentClass( classes( i, Range::all(), 0 ) );
//			// retrieve coordinates of vectors in current class
//			for ( int k = 0; k < R.cols(); k++ ){
//				classCoords( Range :: all(), k ) = R( Range :: all(), k ) * currentClass;
//			}
//			// compute class centroid
//			centroid = mean( classCoords, 1 );
//			// compute ESS
//			scores(i) = sum( currentClass * pow2( classCoords - centroid ) );
//		}
//		cout << "Starting classification score " << sum(scores) << endl;
//
//		// do iterations until convergence (or until maxIters reached)
//		bool notFinished = true;
//		int maxIters = 0;
//		int currentIters = 0;
//		while ( ( notFinished == true ) && ( currentIters < maxIters ) ){
//
//			currentIters++;
//			notFinished = false;
//
//			// for all volumes
//			for ( int i = 0; i < this->alignments.rows(); i++ ){
//				Pixel currentScore = sum( scores );
//				TinyVector< int, 1 > cpc = maxIndex( classification( Range::all(), i, 0 ) );
//				int currentPointClass = cpc[0];
//				if ( classification( currentPointClass, i, 0 ) == 0 ){
//					continue;
//				}
//				Pixel bestScore = currentScore;
//				int bestClass = currentPointClass;
//				Pixel bestScoreOffset = 0;
//				Pixel bestNewCurrentScoreOffset = 0;
//				// place point in each of the remaining classes
//				Pixel newCurrentClassScoreOffset = 0;
//				//for ( int c = 0; c < finalClasses[currentPointClass].size(); c++ ){
//				//	newCurrentClassScoreOffset += pow2( this->alignments( i, finalClasses[currentPointClass][c], 0 ) );
//				//}
//				// compute ESS if we remove point i from current assigned class
//				for ( int j = 0; j < finalClasses[currentPointClass].size(); j++ ){
//					if ( finalClasses[currentPointClass][j] != i ){
//						Pixel csum = 0;
//						for ( int k = 0; k < finalClasses[currentPointClass].size(); k++ ){
//							if ( finalClasses[currentPointClass][k] != i ){
//								csum += this->alignments( finalClasses[currentPointClass][j], finalClasses[currentPointClass][k], 0 );
//							}
//						}
//						newCurrentClassScoreOffset += csum * csum / (finalClasses[currentPointClass].size()-1) / (finalClasses[currentPointClass].size()-1);
//						//newCurrentClassScoreOffset += csum / (finalClasses[currentPointClass].size()-1);
//					}
//				}
//
//				for ( int j = 0; j < finalClasses.size(); j++ ){
//					// try shifting volume i to class j (if not already a member)
//					if ( j != currentPointClass ){
//						Pixel newOtherClassScoreOffset = 0;
//						//for ( int c = 0; c < finalClasses[j].size(); c++ ){
//						//	newOtherClassScoreOffset += pow2( this->alignments( i, finalClasses[j][c], 0 ) );
//						//}
//						// compute ESS if we remove point i from current assigned class
//						Pixel csum;
//						for ( int c = 0; c < finalClasses[j].size(); c++ ){
//							csum = 0;
//							for ( int k = 0; k < finalClasses[j].size(); k++ ){
//								csum += this->alignments( finalClasses[j][c], finalClasses[j][k], 0 );
//							}
//							csum += this->alignments( finalClasses[j][c], i, 0 );
//							newOtherClassScoreOffset += csum * csum / (finalClasses[j].size()+1) / (finalClasses[j].size()+1);
//							//newOtherClassScoreOffset += csum / (finalClasses[j].size()+1);
//						}
//						csum = 0;
//						for ( int k = 0; k < finalClasses[j].size(); k++ ){
//							csum += this->alignments( i, finalClasses[j][k], 0 );
//						}
//						newOtherClassScoreOffset += csum * csum / (finalClasses[j].size()+1) / (finalClasses[j].size()+1);
//						//newOtherClassScoreOffset += csum / (finalClasses[j].size()+1);
//
//						//Pixel tentativeScore = currentScore - newCurrentClassScoreOffset + newOtherClassScoreOffset;
//						Pixel tentativeScore = currentScore - scores(currentPointClass) - scores(j) + newCurrentClassScoreOffset + newOtherClassScoreOffset;
//						if ( tentativeScore < bestScore ){
//							bestScore = tentativeScore;
//							bestClass = j;
//							bestScoreOffset = newOtherClassScoreOffset;
//							bestNewCurrentScoreOffset = newCurrentClassScoreOffset;
//						}
//					}
//				}
//				// make the change permanent
//				if ( bestClass != currentPointClass ){
//					// remove from old class
//					vector< int > :: iterator result = find( finalClasses[currentPointClass].begin(), finalClasses[currentPointClass].end(), i );
//					finalClasses[currentPointClass].erase( result );
//					classification( currentPointClass, i, 0 ) = 0;
//					// add to new class
//					finalClasses[bestClass].push_back( i );
//					classification( bestClass, i, 0 ) = 1;
//					//scores( currentPointClass ) -= newCurrentClassScoreOffset;
//					//scores( bestClass ) += bestScoreOffset;
//					scores( currentPointClass ) = bestNewCurrentScoreOffset;
//					scores( bestClass ) = bestScoreOffset;
//					notFinished = true;
//				}
//			}
//#ifdef WIN32
//			// recompute scores to see the improvement
//			cout << "Current score = " << sum(scores) << endl;
//#endif
//		}
//		// recompute scores to see the improvement
//		cout << "Score at iteration " << currentIters << " = " << sum(scores) << endl;
//	}
//
//}
