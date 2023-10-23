#define NBF_PARALLEL_IMPLEMENTATION_MPI 1

#define NBF_AVERAGE_IN_RECIPROCAL_SPACE

#ifdef WIN32
//	 #define NBF_VERBOSE
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

#include <mxml.h>

#define PIXEL float

int main( int argc, char ** argv )
{
#if 0
	//nbfTimer t;
	//t.start();
	//int rows = 10000;
	//int cols = 10000;
	//double *matrix = new double[rows*cols];
	//double **T = new double *[rows];
	//for(int i = 0; i < rows; i++) {
	//	T[i] = &matrix[i*cols];
	//}

	//for ( int i = 0; i < rows; i++ ){
	//	for ( int j = 0; j < cols; j++ ){
	//		T[i][j] = i + j;
	//	}
	//}

	//double * ev = new double[rows];
	//double **evecMat = new double *[rows*cols];
	//for(int i = 0; i < rows; i++) {
	//	evecMat[i] = &matrix[i*cols];
	//}

	//vtkMath::JacobiN(T, rows, ev, evecMat);
	//t.stop();
	//cout << "t = " << t.elapsedSeconds() << endl;

	//return 0;

	Array< float, 3 > mask;
	nbfMatlabReader mread;
	mread.setFileName(argv[1]);
	mread.read(mask);

	nbfMrcWriter mrcw;
	mrcw.setFileName(argv[2]);
	mrcw.write(mask);
	return 0;

	vtkImageData * window = vtkImageData::New();
	vtkImageGaussianSmooth * filter = vtkImageGaussianSmooth::New();

	// filter down the window functions (to improve spherical harmonics representation)
	nbfVTKInterface::blitzToVtk( mask, window );
	filter->SetInput( window );
	filter->SetRadiusFactors(1.5,1.5,1.5);
	filter->Update();
	nbfVTKInterface::vtkToBlitz( filter->GetOutput(), mask );
	filter->Delete();
	window->Delete();

	

	nbfMatlabWriter mwriter;
	mwriter.setFileName("groel_mask_smooth_i1.blitz");
	mwriter.write(mask);
#endif

	//  Initialize MPI.
	MPI_Init ( &argc, &argv );

	if ( argc != 3 ){
		cout << "Usage: MPI_Classification config.xml list.txt" << endl;
		exit(0);
	}

	int my_id, num_procs;

	//  Get this processes's rank.
	MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

	//  Find out how many processes are available.
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	cout << "Process " << my_id << " is active.\n";

	// initialize XML file readout structures
	FILE *fp;
	mxml_node_t *tree;
	fp = fopen( argv[1], "r");
	tree = mxmlLoadFile(NULL, fp, MXML_NO_CALLBACK);
	fclose(fp);
	mxml_node_t * xmlnode;

	// the metric is the only thing all processes share

	// Setup image filter
	nbfImageFilter< PIXEL, 3 > imfilter;

	// use window
	xmlnode = mxmlFindElement(tree, tree, "use_image_window", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field use_image_window not present." << endl;
		return 1;
	}
	PIXEL use_image_window = atof( xmlnode->child->value.text.string );

	if ( use_image_window == 0 ){
		imfilter.windowOff();
	} else {
		PIXEL inner = max( 0.0, use_image_window - .05 );
		PIXEL outer = min( 1.0, use_image_window + .05 );
		imfilter.setMaskSize( inner, outer );
	}
	
	if ( use_image_window == 0 ){
		imfilter.windowOff();
	}

	// padding
	xmlnode = mxmlFindElement(tree, tree, "padding_size", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field padding_size not present." << endl;
		return 1;
	}
	int padding_size = max( 1, atoi( xmlnode->child->value.text.string ) );
	
	imfilter.paddingOn(padding_size);

	// median filter
	xmlnode = mxmlFindElement(tree, tree, "median_filter", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field median_filter not present." << endl;
		return 1;
	}
	int median_filter = atoi( xmlnode->child->value.text.string );
	if ( median_filter > 0 ){
		imfilter.medianFilterOn(median_filter);
	}

	// Setup Fourier filter
	nbfFourierFilter< PIXEL, 3 > fffilter;

	xmlnode = mxmlFindElement(tree, tree, "low_pass_cutoff", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field low_pass_cutoff not present." << endl;
		return 1;
	}
	PIXEL lp_cutoff = atof( xmlnode->child->value.text.string );
	xmlnode = mxmlFindElement(tree, tree, "low_pass_decay", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field low_pass_decay not present." << endl;
		return 1;
	}
	PIXEL lp_decay = atof( xmlnode->child->value.text.string );
	xmlnode = mxmlFindElement(tree, tree, "high_pass_cutoff", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field high_pass_cutoff not present." << endl;
		return 1;
	}
	PIXEL hp_cutoff = atof( xmlnode->child->value.text.string );
	xmlnode = mxmlFindElement(tree, tree, "high_pass_decay", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field high_pass_decay not present." << endl;
		return 1;
	}
	PIXEL hp_decay = atof( xmlnode->child->value.text.string );

	if ( hp_cutoff > lp_cutoff ){
		cout << "ERROR - Check the bandpass filter range (" << hp_cutoff << ">" << lp_cutoff << ")" << endl;
		return 1;
	}
	fffilter.bandPassOn( hp_cutoff, lp_cutoff, hp_decay, lp_decay );

	// Setup metric configuration
	nbfProjectionRotationMetric3D< PIXEL > metric( &imfilter, &fffilter );

	xmlnode = mxmlFindElement(tree, tree, "number_of_candidate_peaks_to_search", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field number_of_candidate_peaks_to_search not present." << endl;
		return 1;
	}
	int candidates_to_search = atof( xmlnode->child->value.text.string );
	metric.setNumberOfCandidatePeaksToSearch( candidates_to_search );

	xmlnode = mxmlFindElement(tree, tree, "rotation_search_restriction", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field rotation_search_restriction not present." << endl;
		return 1;
	}
	PIXEL rotation_search_restriction = atof( xmlnode->child->value.text.string );
	metric.setRotationSearchRestriction( rotation_search_restriction );

	xmlnode = mxmlFindElement(tree, tree, "translation_search_restriction", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field translation_search_restriction not present." << endl;
		return 1;
	}
	PIXEL translation_search_restriction = atof( xmlnode->child->value.text.string );
	metric.setTranslationSearchRestriction( translation_search_restriction );

	xmlnode = mxmlFindElement(tree, tree, "number_of_alignment_candidates", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field number_of_alignment_candidates not present." << endl;
		return 1;
	}
	int number_of_alignment_candidates = atoi( xmlnode->child->value.text.string );
	metric.setNumberOfCandidates( number_of_alignment_candidates );
	
	if ( my_id == 0 ){ // Master process

		cout << "MPI - Master process:\n";
		cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
		cout << "  The number of processes is " << num_procs << "\n";

		xmlnode = mxmlFindElement(tree, tree, "volumes", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field volumes not present." << endl;
			return 1;
		}
		char * spikesFileName = xmlnode->child->value.text.string;

		xmlnode = mxmlFindElement(tree, tree, "pattern", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field pattern not present." << endl;
			return 1;
		}
		char * pattern = xmlnode->child->value.text.string;

		xmlnode = mxmlFindElement(tree, tree, "starting_iteration", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field starting_iteration not present." << endl;
			return 1;
		}
		int pattern_iteration = atoi( xmlnode->child->value.text.string );

		xmlnode = mxmlFindElement(tree, tree, "initial_classes", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field initial_classes not present." << endl;
			return 1;
		}
		char * initialClasses = NULL;
		if ( xmlnode->child != NULL ){
			initialClasses = xmlnode->child->value.text.string;
		}

		xmlnode = mxmlFindElement(tree, tree, "hierarchical_classes", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field hierarchical_classes not present." << endl;
			return 1;
		}
		PIXEL hierarchicalThreshold = atof( xmlnode->child->value.text.string );

		xmlnode = mxmlFindElement(tree, tree, "build_hierarchical_tree", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field build_hierarchical_tree not present." << endl;
			return 1;
		}
		PIXEL buildHierarchicalTree = atof( xmlnode->child->value.text.string );

		xmlnode = mxmlFindElement(tree, tree, "loop_iterations", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field loop_iterations not present." << endl;
			return 1;
		}
		int maxIter = atoi( xmlnode->child->value.text.string );

		// umbral usado como criterio para juntar clases - 
		xmlnode = mxmlFindElement(tree, tree, "min_elements_per_class", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field min_elements_per_class not present." << endl;
			return 1;
		}
		int minElementsPerClass = atoi( xmlnode->child->value.text.string );

		vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
		nbfWedgedSubImage3D< PIXEL > :: read( spikesFileName, volumeList );

		cout << "\nTotal " << volumeList.size() << " volumes to process." << endl;

		// size to cutout volumes
		xmlnode = mxmlFindElement(tree, tree, "volume_size", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field volume_size not present." << endl;
			return 1;
		}
		int volumeSize = atoi( xmlnode->child->value.text.string );

		// set volume size and cutting offset to center of subvolume
		for ( int i = 0; i < volumeList.size(); i++ ){
			TinyVector< int, 3 > tsize( volumeSize, volumeSize, volumeSize );
			volumeList[i].setDimensions( tsize );
			volumeList[i].setCutOffset( volumeSize / 2.0 );
		}

		// size to cutout volumes
		xmlnode = mxmlFindElement(tree, tree, "use_symmetrization", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field use_symmetrization not present." << endl;
			return 1;
		}
		int symmetry = atoi( xmlnode->child->value.text.string );

		// size to cutout volumes
		xmlnode = mxmlFindElement(tree, tree, "use_only_top_percentage", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field use_only_top_percentage not present." << endl;
			return 1;
		}
		PIXEL top_percentage = atof( xmlnode->child->value.text.string );

		vector< nbfWedgedSubImage3D< PIXEL > > newVolumeList;
		nbfWedgedSubImage3D< PIXEL > :: read( argv[2], newVolumeList );

		//vtkTransform * t = vtkTransform :: New();
		//t->Translate(0,0,5);
		//newVolumeList[0].setTransform(t);
		//vtkImageData * data = vtkImageData :: New();
		//newVolumeList[0].getImage(data,t);
		//t->Delete();

		//Array< double, 3 > A, B;
		//nbfVTKInterface::vtkToBlitzReference( data, A );

		//A.transposeSelf(thirdDim,secondDim,firstDim);
		//B.resize( A.shape() );
		//B = A.reverse(secondDim);
		//nbfVTKInterface::blitzToVtk(B,data);

		//vtkImageCast * cast = vtkImageCast::New();
		//cast->SetOutputScalarTypeToFloat();
		//cast->SetInput( data );
		//cast->Update();

		//nbfMrcWriter mrcw;
		//mrcw.setFileName( "tmp.mrc" );
		//mrcw.write( cast->GetOutput() );
		//cast->Delete();

		//data->Delete();
		//return 0;

		vector< TinyVector< int, 2 > > positions;
		vector< nbfWedgedImage3D< PIXEL > * > plist1, plist2;

		plist1.push_back( &(newVolumeList[0]) );

		for ( int j = 1; j < newVolumeList.size(); j++ ){
			plist2.push_back( &(newVolumeList[j]) );

			TinyVector< int, 2 > currentPosition( plist1.size() - 1, plist2.size() - 1 );
			positions.push_back( currentPosition );
		}

		Array< PIXEL, 3 > allDistances( plist2.size(), 19, metric.getNumberOfCandidates() );
		allDistances = -1;
		metric.getDistances( plist1, plist2, positions, allDistances, 0 );
		
		nbfMatlabWriter w;
		w.setFileName("p.matlab");
		w.write( allDistances( Range::all(), 0, 0 ) );

		cout << "Distances to XRAY = " << allDistances( Range::all(), 0, 0 ) << endl;

		for ( int j = 1; j < newVolumeList.size(); j++ ){
			// apply alignments and save volume

			double matrix[16];
			for ( int k = 0; k < 16; k++ ){
				matrix[k] = allDistances( j - 1, k + 3, 0 );
			}
			vtkMatrix4x4 * mat = vtkMatrix4x4 :: New();
			mat->DeepCopy( matrix );
			vtkTransform * t = vtkTransform :: New();
			t->SetMatrix( mat );
			mat->Delete();
			//newVolumeList[j].setTransform(t);
			vtkImageData * data = vtkImageData :: New();
			newVolumeList[j].getImage(data,&imfilter,t);
			t->Delete();

			Array< double, 3 > A, B;
			nbfVTKInterface::vtkToBlitzReference( data, A );

			A.transposeSelf(thirdDim,secondDim,firstDim);
			B.resize( A.shape() );
			B = A.reverse(secondDim);
			nbfVTKInterface::blitzToVtk(B,data);

			vtkImageCast * cast = vtkImageCast::New();
			cast->SetOutputScalarTypeToFloat();
			cast->SetInput( data );
			cast->Update();

			nbfMrcWriter mrcw;
			stringstream file;
			file << newVolumeList[j].getFileName() << "_ali.mrc";
			mrcw.setFileName( file.str().c_str() );
			mrcw.write( cast->GetOutput() );
			cast->Delete();

			data->Delete();

		}

		// end MPI
		metric.finalizeMPI();
	} else {
		metric.slaveMPI();
	}
	MPI_Finalize();

	return 0;
}