#define NBF_PARALLEL_IMPLEMENTATION_MPI 1

#define NBF_AVERAGE_IN_RECIPROCAL_SPACE

//#ifdef WIN32
//	 #define NBF_VERBOSE
//#endif

#define BZ_GENERATE_GLOBAL_INSTANCES

#include "mpi.h"

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

#include <blitz/bzdebug.h>        // Test suite globals

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

#include <nbfCylindricalDomain3.h>

#include <mxml.h>

extern "C" {
#include <svdlib.h>
}

#include <nbfTimer.h>

#define PIXEL float

int main( int argc, char ** argv )
{	
	/*
	double matrix[16];
	
	matrix[0]=0.8925;
	matrix[1]=-0.3848;
	matrix[2]=0.2249;
	matrix[3]=0;
	matrix[4]=0.3642;
	matrix[5]=0.9224;
	matrix[6]=0.1283;
	matrix[7]=0;
	matrix[8]=-0.2568;
	matrix[9]=-0.03298;
	matrix[10]=0.9659;
	matrix[11]=0;
	matrix[12]=0;
	matrix[13]=0;
	matrix[14]=0;
	matrix[15]=1;

	//matrix[0]=0.6601;
	//matrix[1]=0.1874;
	//matrix[2]=0.7275;
	//matrix[3]=0;
	//matrix[4]=-0.2894;
	//matrix[5]=0.9571;
	//matrix[6]=0.01607;
	//matrix[7]=0;
	//matrix[8]=-.6932;
	//matrix[9]=-0.2211;
	//matrix[10]=0.686;
	//matrix[11]=0;
	//matrix[12]=0;
	//matrix[13]=0;
	//matrix[14]=0;
	//matrix[15]=1;
	vtkMatrix4x4 * rot = vtkMatrix4x4::New();
	rot->DeepCopy(matrix);
	cout << *rot << endl;

	vtkMatrix4x4 * trans = vtkMatrix4x4::New();
	matrix[0]=1;matrix[1]=matrix[2]=0; matrix[3]=16.37;
	matrix[4]=0;matrix[5]=1;matrix[6]=0; matrix[7]=12.79;
	matrix[8]=matrix[9]=0;matrix[10]=1; matrix[11]=-8.101;
	matrix[12]=matrix[13]=matrix[14]=0;	matrix[15]=1;
	//matrix[0]=1;matrix[1]=matrix[2]=0; matrix[3]=-237.4;
	//matrix[4]=0;matrix[5]=1;matrix[6]=0; matrix[7]=-45.84;
	//matrix[8]=matrix[9]=0;matrix[10]=1; matrix[11]=-4.3497;
	//matrix[12]=matrix[13]=matrix[14]=0;	matrix[15]=1;
	trans->DeepCopy(matrix);
	cout << *trans << endl;

	vtkMatrix4x4 * comp = vtkMatrix4x4::New();
	vtkMatrix4x4::Multiply4x4(trans,rot,comp);
	cout << *comp << endl;

	vtkTransform * t = vtkTransform :: New();
	t->SetMatrix( comp );					
	double angles[3];
	t->GetOrientation(angles);
	cout << angles[0] << ", " << angles[1] << ", " << angles[2] << endl;
*/

    /*
    Array< double, 3 > A;
    vtkImageData * data = vtkImageData :: New();
    if ( false ){
        nbfMrcReader reader;
        reader.setFileName( "Tomo_BaL_20140208_20140614_global_average_symmetrized.mrc" );
        reader.read( data );
    } else { 
		vector< nbfWedgedSubImage3D< PIXEL > > volumes;
		nbfWedgedSubImage3D< PIXEL > :: read( "average.txt", volumes );
        volumes[0].getImage( data );
    }
    nbfVTKInterface :: vtkToBlitz( data, A );        

	nbfMatlabWriter w;
	w.setFileName("p.matlab");
    // w.write( A );
    exit(0);
    */

	int my_id, num_procs;

// b2 settings
#if 1

    int    namelen;
    char   processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Init ( &argc, &argv );
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Get_processor_name(processor_name,&namelen);
    // cout << "Process " << my_id << " is active.\n";

#else


	//  Initialize MPI.
	MPI_Init ( &argc, &argv );

	if ( argc != 2 ){
		cout << argv[0] << "\n";
		cout << " Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
		cout << "Usage: MPI_Classification config.xml" << endl;
		exit(0);
	}

	//  Get this processes's rank.
	MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

	//  Find out how many processes are available.
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

#endif

    // cout << "Process " << my_id << " is active.\n";

	// initialize XML file readout structures
	FILE *fp;
	mxml_node_t *tree;
	fp = fopen( argv[1], "r");

	if ( ( my_id == 0 ) && ( fp == NULL ) ){
		cout << "Failed to open configuration file " << argv[1] << endl;
		exit(0);
	}

	tree = mxmlLoadFile(NULL, fp, MXML_NO_CALLBACK);
	fclose(fp);
	mxml_node_t * xmlnode;

	// runing mode
	xmlnode = mxmlFindElement(tree, tree, "mode", NULL, NULL, MXML_DESCEND);
	if ( ( my_id == 0 ) && ( xmlnode == NULL ) ){
		cerr << "ERROR - Invalid configuration file. Field mode not present." << endl;
		return 1;
	}
	int running_mode = atoi( xmlnode->child->value.text.string );
	if ( ( my_id == 0 ) && ( abs( running_mode > 4 ) ) ){
		cerr << "ERROR - Invalid operation mode specified." << endl;
		return 1;
	}

	PIXEL use_image_window_x, use_image_window_y, use_image_window_z, use_image_window_sigma;

	PIXEL lp_cutoff, lp_decay, hp_cutoff, hp_decay;

	PIXEL bfactor;
		
	stringstream maskFile;

	PIXEL maskFileTh;

	PIXEL rotation_search_restriction = 180, translation_search_restriction;

	int symmetry;

	xmlnode = mxmlFindElement(tree, tree, "mra_raw_square_window_size", NULL, NULL, MXML_DESCEND);
	int squareMaskSize = 5;
	if ( xmlnode != NULL ){
		squareMaskSize = atoi( xmlnode->child->value.text.string );
	}

	switch ( abs( running_mode ) ){
			case NBF_LOOP_CLUSTERING_CLASS:
			case NBF_LOOP_CLUSTERING_CTF:

				// class symmetry factor
				xmlnode = mxmlFindElement(tree, tree, "class_use_symmetrization", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field class_use_symmetrization not present." << endl;
					return 1;
				}
				symmetry = atoi( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "class_image_window_x", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field class_image_window_x not present." << endl;
					return 1;
				}
				use_image_window_x = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "class_image_window_y", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field class_image_window_y not present." << endl;
					return 1;
				}
				use_image_window_y = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "class_image_window_z", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field class_image_window_z not present." << endl;
					return 1;
				}
				use_image_window_z = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "class_image_window_sigma", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field class_image_window_sigma not present." << endl;
					return 1;
				}
				use_image_window_sigma = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "class_file_image_window", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field class_file_image_window not present." << endl;
					return 1;
				}
				if ( xmlnode->child != NULL ){
					maskFile << xmlnode->child->value.text.string;
				}

				xmlnode = mxmlFindElement(tree, tree, "class_file_image_threshold", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field class_file_image_threshold not present." << endl;
					return 1;
				}
				maskFileTh = atof(xmlnode->child->value.text.string);

				xmlnode = mxmlFindElement(tree, tree, "class_low_pass_cutoff", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field class_low_pass_cutoff not present." << endl;
					return 1;
				}
				lp_cutoff = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "class_low_pass_decay", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field class_low_pass_decay not present." << endl;
					return 1;
				}
				lp_decay = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "class_high_pass_cutoff", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field class_high_pass_cutoff not present." << endl;
					return 1;
				}
				hp_cutoff = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "class_high_pass_decay", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field class_high_pass_decay not present." << endl;
					return 1;
				}
				hp_decay = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "class_bfactor", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "WARNING - Field class_bfactor not present in configuration file. Defaulting to 0." << endl;
					bfactor = 0;
				} else {
					bfactor = atof( xmlnode->child->value.text.string );
				}

				break;
			case NBF_LOOP_CLUSTERING_REFINE:

				// refine symmetry factor
				xmlnode = mxmlFindElement(tree, tree, "refine_use_symmetrization", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_use_symmetrization not present." << endl;
					return 1;
				}
				symmetry = atoi( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "refine_image_window_x", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_image_window_x not present." << endl;
					return 1;
				}
				use_image_window_x = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "refine_image_window_y", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_image_window_y not present." << endl;
					return 1;
				}
				use_image_window_y = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "refine_image_window_z", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_image_window_z not present." << endl;
					return 1;
				}
				use_image_window_z = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "refine_image_window_sigma", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_image_window_sigma not present." << endl;
					return 1;
				}
				use_image_window_sigma = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "refine_file_image_window", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_file_image_window not present." << endl;
					return 1;
				}
				if ( xmlnode->child != NULL ){
					maskFile << xmlnode->child->value.text.string;
				}

				xmlnode = mxmlFindElement(tree, tree, "refine_file_image_threshold", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_file_image_threshold not present." << endl;
					return 1;
				}
				maskFileTh = atof(xmlnode->child->value.text.string);

				xmlnode = mxmlFindElement(tree, tree, "refine_low_pass_cutoff", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_low_pass_cutoff not present." << endl;
					return 1;
				}
				lp_cutoff = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "refine_low_pass_decay", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_low_pass_decay not present." << endl;
					return 1;
				}
				lp_decay = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "refine_high_pass_cutoff", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_high_pass_cutoff not present." << endl;
					return 1;
				}
				hp_cutoff = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "refine_high_pass_decay", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_high_pass_decay not present." << endl;
					return 1;
				}
				hp_decay = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "refine_bfactor", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "WARNING - Field refine_bfactor not present in configuration file. Defaulting to 0." << endl;
					bfactor = 0;
				} else {
					bfactor = atof( xmlnode->child->value.text.string );
				}

				xmlnode = mxmlFindElement(tree, tree, "refine_out_of_plane_search_range", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_out_of_plane_search_range not present." << endl;
					return 1;
				}
				rotation_search_restriction = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "refine_shifts_tolerance", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field refine_shifts_tolerance not present." << endl;
					return 1;
				}
				translation_search_restriction = atof( xmlnode->child->value.text.string );

				break;
			case NBF_LOOP_CLUSTERING_CENTER:
			case NBF_LOOP_CLUSTERING_MRA:

				// mra symmetry factor
				xmlnode = mxmlFindElement(tree, tree, "mra_use_symmetrization", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_use_symmetrization not present." << endl;
					return 1;
				}
				symmetry = atoi( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "mra_image_window_x", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_image_window_x not present." << endl;
					return 1;
				}
				use_image_window_x = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "mra_image_window_y", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_image_window_y not present." << endl;
					return 1;
				}
				use_image_window_y = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "mra_image_window_z", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_image_window_z not present." << endl;
					return 1;
				}
				use_image_window_z = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "mra_image_window_sigma", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_image_window_sigma not present." << endl;
					return 1;
				}
				use_image_window_sigma = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "mra_file_image_window", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_file_image_window not present." << endl;
					return 1;
				}
				if ( xmlnode->child != NULL ){
					maskFile << xmlnode->child->value.text.string;
				}

				xmlnode = mxmlFindElement(tree, tree, "mra_file_image_threshold", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_file_image_threshold not present." << endl;
					return 1;
				}
				if ( xmlnode->child != NULL ){
					maskFileTh = atof(xmlnode->child->value.text.string);
				}

				xmlnode = mxmlFindElement(tree, tree, "mra_low_pass_cutoff", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_low_pass_cutoff not present." << endl;
					return 1;
				}
				lp_cutoff = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "mra_low_pass_decay", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_low_pass_decay not present." << endl;
					return 1;
				}
				lp_decay = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "mra_high_pass_cutoff", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_high_pass_cutoff not present." << endl;
					return 1;
				}
				hp_cutoff = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "mra_high_pass_decay", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_high_pass_decay not present." << endl;
					return 1;
				}
				hp_decay = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "mra_bfactor", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "WARNING - Field mra_bfactor not present in configuration file. Defaulting to 0." << endl;
					bfactor = 0;
				} else {
					bfactor = atof( xmlnode->child->value.text.string );
				}

				xmlnode = mxmlFindElement(tree, tree, "mra_out_of_plane_search_range", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_out_of_plane_search_range not present." << endl;
					return 1;
				}
				rotation_search_restriction = atof( xmlnode->child->value.text.string );

				xmlnode = mxmlFindElement(tree, tree, "mra_shifts_tolerance", NULL, NULL, MXML_DESCEND);
				if ( xmlnode == NULL ){
					cerr << "ERROR - Invalid configuration file. Field mra_shifts_tolerance not present." << endl;
					return 1;
				}
				translation_search_restriction = atof( xmlnode->child->value.text.string );

				break;
	}

	// the metric is the only thing all processes share

	// Setup image filter
	nbfImageFilter< PIXEL, 3 > imfilter;

	if ( use_image_window_x == 0 ){
		imfilter.maskOff();
	} else {
		// the last parameters sets to use cylyndrical mask
		imfilter.setMaskSize( use_image_window_x, use_image_window_y, use_image_window_z, use_image_window_sigma, use_image_window_z < 0 );
	}

	if ( maskFile.str().size() > 0 ){
		imfilter.setMaskFile( maskFile, maskFileTh );
		imfilter.fileMaskOn();
	}

	// padding
	xmlnode = mxmlFindElement(tree, tree, "padding_size", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field padding_size not present." << endl;
		return 1;
	}
	int padding_size = max( 1, atoi( xmlnode->child->value.text.string ) );
	
	imfilter.paddingOn( padding_size );

	// gaussian filter
	xmlnode = mxmlFindElement(tree, tree, "gaussian_filter", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field gaussian_filter not present." << endl;
		return 1;
	}
	int gaussian_filter = atoi( xmlnode->child->value.text.string );
	if ( gaussian_filter > 0 ){
		imfilter.gaussianFilterOn( gaussian_filter );
	}

	// size to cutout volumes
	xmlnode = mxmlFindElement(tree, tree, "volume_size", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field volume_size not present." << endl;
		return 1;
	}
	int volumeSize = atoi( xmlnode->child->value.text.string );

	// set square mask size for raw volumes
	squareMaskSize = floor( squareMaskSize + 0.0 );
	//squareMaskSize = floor( ( volumeSize - squareMaskSize ) / 2.0 );
	imfilter.squareMaskSize = squareMaskSize;

	// Setup Fourier filter
	nbfFourierFilter< PIXEL, 3 > fffilter;

	if ( hp_cutoff > lp_cutoff ){
		cout << "ERROR - Check the bandpass filter range (" << hp_cutoff << ">" << lp_cutoff << ")" << endl;
		return 1;
	}
	fffilter.bandPassOn( hp_cutoff, lp_cutoff, hp_decay, lp_decay );
	
	if ( bfactor != 0 ){
		fffilter.bfactorOn( bfactor );
	}

	// Setup metric configuration
	nbfProjectionRotationMetric3D< PIXEL > metric( &imfilter, &fffilter );

	xmlnode = mxmlFindElement(tree, tree, "number_of_candidate_peaks_to_search", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field number_of_candidate_peaks_to_search not present." << endl;
		return 1;
	}
	int candidates_to_search = atof( xmlnode->child->value.text.string );
	metric.setNumberOfCandidatePeaksToSearch( candidates_to_search );

	xmlnode = mxmlFindElement(tree, tree, "use_missing_wedge", NULL, NULL, MXML_DESCEND);
	bool use_missing_wedge = false;
	if ( xmlnode != NULL ){
		use_missing_wedge = atof( xmlnode->child->value.text.string );
	}
	metric.setMissingWedgeCompensation( use_missing_wedge );

	if ( abs( running_mode ) != NBF_LOOP_CLUSTERING_CLASS ){
		metric.setRotationSearchRestriction( rotation_search_restriction );
		metric.setTranslationSearchRestriction( translation_search_restriction );
	}

	xmlnode = mxmlFindElement(tree, tree, "number_of_alignment_candidates", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field number_of_alignment_candidates not present." << endl;
		return 1;
	}
	int number_of_alignment_candidates = atoi( xmlnode->child->value.text.string );
	metric.setNumberOfCandidates( number_of_alignment_candidates );
	
	// Metric configuration
	xmlnode = mxmlFindElement(tree, tree, "use_mutual_correlation", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field use_mutual_correlation not present." << endl;
		return 1;
	}
	int use_mutual_correlation = atoi( xmlnode->child->value.text.string );
	metric.setToUseMutualCorrelation( use_mutual_correlation > 0 );

	// Metric configuration
	xmlnode = mxmlFindElement(tree, tree, "use_wedge_overlap_normalization", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field use_wedge_overlap_normalization not present." << endl;
		return 1;
	}
	int use_wedge_overlap_normalization = atoi( xmlnode->child->value.text.string );
	metric.setToComputeOverlapNormalizedDistances( use_wedge_overlap_normalization > 0 );

	// save parameter values
	xmlnode = mxmlFindElement(tree, tree, "pattern", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field pattern not present." << endl;
		return 1;
	}
	char * pattern = xmlnode->child->value.text.string;

	xmlnode = mxmlFindElement(tree, tree, "current_iteration", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field current_iteration not present." << endl;
		return 1;
	}
	int pattern_iteration = atoi( xmlnode->child->value.text.string );

	std::ifstream src(argv[1]);
	stringstream parameterFile;
	if ( pattern_iteration < 10 ){
		parameterFile << pattern << "_iteration_00" << pattern_iteration << "_mode_" << abs(running_mode) << ".xml";
	} else if ( pattern_iteration < 100 ){
		parameterFile << pattern << "_iteration_0" << pattern_iteration << "_mode_" << abs(running_mode)  << ".xml";
	} else {
		parameterFile << pattern << "_iteration_" << pattern_iteration << "_mode_" << abs(running_mode)  << ".xml";
	}
	std::ofstream dst(parameterFile.str().c_str());
	dst << src.rdbuf();
	dst.close();

	if ( my_id == 0 ){ // Master process

		cout << argv[0] << "\n";
		cout << " Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
		cout << " The number of processes is " << num_procs << "\n";

		if ( pattern_iteration < 10 ){
			cout << " Executing run " << pattern << "_iteration_00" << pattern_iteration << "_mode_" << abs(running_mode) << endl;
		} else if ( pattern_iteration < 100 ){
			cout << " Executing run " << pattern << "_iteration_0" << pattern_iteration << "_mode_" << abs(running_mode) << endl;
		} else {
			cout << " Executing run " << pattern << "_iteration_" << pattern_iteration << "_mode_" << abs(running_mode) << endl;
		}

		xmlnode = mxmlFindElement(tree, tree, "volumes", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field volumes not present." << endl;
			return 1;
		}
		char * spikesFileName = xmlnode->child->value.text.string;

		xmlnode = mxmlFindElement(tree, tree, "mra_initial_classes", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field mra_initial_classes not present." << endl;
			return 1;
		}
		char * initialClasses = NULL;
		if ( xmlnode->child != NULL ){
			initialClasses = xmlnode->child->value.text.string;
		}

		// retrieve axis correction
		xmlnode = mxmlFindElement(tree, tree, "mra_axis_correction_x", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field mra_axis_correction_x not present." << endl;
			return 1;
		}
		PIXEL corr_x = atoi( xmlnode->child->value.text.string );

		xmlnode = mxmlFindElement(tree, tree, "mra_axis_correction_y", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field mra_axis_correction_y not present." << endl;
			return 1;
		}
		PIXEL corr_y = atoi( xmlnode->child->value.text.string );

		xmlnode = mxmlFindElement(tree, tree, "mra_axis_correction_z", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field mra_axis_correction_z not present." << endl;
			return 1;
		}
		PIXEL corr_z = atoi( xmlnode->child->value.text.string );

		xmlnode = mxmlFindElement(tree, tree, "class_number_of_classes", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field class_number_of_classes not present." << endl;
			return 1;
		}
		int hierarchicalClasses = atoi( xmlnode->child->value.text.string );

		int hierarchicalClassesHigh = hierarchicalClasses;
		xmlnode = mxmlFindElement(tree, tree, "class_number_of_classes_high", NULL, NULL, MXML_DESCEND);
		if ( xmlnode != NULL ){
			hierarchicalClassesHigh = atoi( xmlnode->child->value.text.string );
		}

		xmlnode = mxmlFindElement(tree, tree, "loop_iterations", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field loop_iterations not present." << endl;
			return 1;
		}
		int maxIter = atoi( xmlnode->child->value.text.string );

		xmlnode = mxmlFindElement(tree, tree, "apply_rotational_symmetry", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field apply_rotational_symmetry not present." << endl;
			return 1;
		}
		int apply_rotational_symmetry = atoi( xmlnode->child->value.text.string );

		xmlnode = mxmlFindElement(tree, tree, "compute_variance_map", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field compute_variance_map not present." << endl;
			return 1;
		}
		int compute_variance_map = atoi( xmlnode->child->value.text.string );

		vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
		nbfWedgedSubImage3D< PIXEL > :: read( spikesFileName, volumeList );

		if ( volumeList.size() == 0 ){
			cout << "Error reading sub-volumes from file " << spikesFileName << endl;
			return 1;
		} else {
			cout << "\nTotal " << volumeList.size() << " volumes to process." << endl;
		}

		// position to center volumes
		xmlnode = mxmlFindElement(tree, tree, "volume_cut_offset", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field volume_cut_offset not present." << endl;
			return 1;
		}
		int volumeCutOffset = atoi( xmlnode->child->value.text.string );

		// set volume size and cut offset along normal direction
		TinyVector< int, 3 > tsize( volumeSize, volumeSize, volumeSize );
		for ( int i = 0; i < volumeList.size(); i++ ){
			volumeList[i].setDimensions( tsize );
			volumeList[i].setCutOffset( volumeCutOffset );
		}

		// membership threshold by distance
		xmlnode = mxmlFindElement(tree, tree, "class_cutoff_percentage_pre", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field class_cutoff_percentage_pre not present." << endl;
			return 1;
		}
		PIXEL top_percentage_pre = atof( xmlnode->child->value.text.string );

		// membership threshold by distance
		xmlnode = mxmlFindElement(tree, tree, "class_cutoff_percentage", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field class_cutoff_percentage not present." << endl;
			return 1;
		}
		PIXEL top_percentage = atof( xmlnode->child->value.text.string );

		// membership threshold by distance
		xmlnode = mxmlFindElement(tree, tree, "class_cutoff_selection", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field class_cutoff_selection not present." << endl;
			return 1;
		}
		PIXEL top_percentage_selection = atof( xmlnode->child->value.text.string );

		// membership threshold by distance
		xmlnode = mxmlFindElement(tree, tree, "class_binning_factor", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field class_binning_factor not present." << endl;
			return 1;
		}
		PIXEL class_binning_factor = atof( xmlnode->child->value.text.string );

		// refinement iterations
		xmlnode = mxmlFindElement(tree, tree, "refinement_iterations", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field refinement_iterations not present." << endl;
			return 1;
		}
		int refinement = atoi( xmlnode->child->value.text.string );

		// refinement iterations
		xmlnode = mxmlFindElement(tree, tree, "refinement_fsc", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field refinement_fsc not present." << endl;
			return 1;
		}
		int fsc = atoi( xmlnode->child->value.text.string );

		// reset alignments before mra
		xmlnode = mxmlFindElement(tree, tree, "mra_reset_initial_alignments", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field mra_reset_initial_alignments not present." << endl;
			return 1;
		}
		int resetAlignments = atoi( xmlnode->child->value.text.string );

		// alignment mode
		xmlnode = mxmlFindElement(tree, tree, "alignment_mode", NULL, NULL, MXML_DESCEND);
		if ( xmlnode == NULL ){
			cerr << "ERROR - Invalid configuration file. Field alignment_mode not present." << endl;
			return 1;
		}
		int alignment_mode = atoi( xmlnode->child->value.text.string );

		Array< PIXEL, 4 > alignments;

		// Do clustering
		nbfLoopClustering< PIXEL > cluster;
		cluster.setFileHeader( pattern, pattern_iteration );
		cluster.setInput( volumeList );
		cluster.setMetric( &metric, alignments );
		cluster.setHierarchicalClasses( hierarchicalClasses, hierarchicalClassesHigh );
		cluster.setDistanceTopCutoff( top_percentage_pre, top_percentage, top_percentage_selection );
		cluster.setSymmetryFactor( symmetry );
		cluster.setInitialClasses( initialClasses );
		cluster.setRefinementIterations( refinement );
		cluster.setRunningMode( running_mode );
		cluster.binFactorForClassification = class_binning_factor;

		// old options
		cluster.useRealRepresentationOn(10);
		cluster.setHierarchicalToBuildEntireTree(1);
		cluster.setHierarchicalMinOverlap(0.0);

		if ( running_mode != 0 ){
			cluster.setMaxIterations( 1 );
			Array< PIXEL, 3 > classes;
			switch ( running_mode ){
				case 1:
					cluster.doClassAndSelect();
					break;
				case -1:
					cluster.doClass();
					break;
				case 2:
					cluster.doRefine( apply_rotational_symmetry > 0, fsc > 0, alignment_mode, volumeCutOffset );
					break;
				case 3:
					cluster.doMra( resetAlignments > 0, corr_x, corr_y, corr_z, alignment_mode );
					break;
				case 4:
					cluster.doCTFClass();
					break;
				case 5:
					cluster.do2DClass();
					break;
				case 6:
					cluster.do2DRefine();
					break;
				case 7:
					cluster.do2DMra();
					break;
			}
			//cluster.execute(classes);
		} else {

			// Pre-center data to rotationally symmetric global average

			alignments.resize( 1, volumeList.size(), 17, 1 );
			alignments = 1;

			for ( int iter = 0; iter <= maxIter; iter++ ){

				// store global average
				nbfWedgedAverageImage3D< PIXEL > globalAverage( volumeList );
				Array< PIXEL, 3 > globalAlignments( volumeList.size(), 17, 1 );
				globalAlignments = 1;

				// apply cutoff threshold
				if ( top_percentage < 1 ){
					// find threshold
					vector< PIXEL > distances;
					for ( int h = 0; h < alignments.cols(); h++ ){
						//distances.push_back( volumeList[h].getCutOffset() );
						distances.push_back( alignments( 0, h, 0, 0 ) );
					}
					sort( distances.begin(), distances.end() );
					PIXEL cutoffDistance = distances[ ( distances.size() - 1 ) * top_percentage ];
					cout << "Setting cutoff distance to " << cutoffDistance << endl;
					globalAlignments( Range :: all(), 0, 0 ) = where( alignments( 0, Range :: all(), 0, 0 ) <= cutoffDistance, 1, 0 );
					cout << "Current average has " << sum( globalAlignments( Range :: all(), 0, 0 ) ) << " volumes." << endl;
				}

				//if ( iter == 0 ){
				//	globalAlignments = 0;
				//	globalAlignments( 0, 0, 0 ) = 1;
				//}

				double matrix[16];
				for ( int i = 0; i < volumeList.size(); i++ ){
					vtkTransform * t = vtkTransform::New();
					volumeList[i].getTransform(t);
					vtkMatrix4x4 :: DeepCopy( matrix, t->GetMatrix() );
					t->Delete();
					for ( int e = 0; e < 16; e++ ){
						globalAlignments( i, 1 + e, 0 ) = matrix[e];
					}
				}
				globalAverage.setAlignments( globalAlignments );
				metric.getImage( globalAverage );

				// convert data to blitz array for manipulation
				vtkImageData * data = vtkImageData::New();
				globalAverage.getImage( data );

				Array< double, 3 > A;
				nbfVTKInterface :: vtkToBlitzReference( data, A );
				A *= -1;

				vtkImageCast * cast = vtkImageCast::New();
				cast->SetOutputScalarTypeToFloat();
				cast->SetInput( data );
				cast->Update();

				nbfMrcWriter mrcw;
				stringstream fname;
				if ( iter == 0 ){
					fname << pattern << "_global_average.mrc";
				} else {
					fname << pattern << "_global_average_" << iter << ".mrc";
				}
				mrcw.setFileName( fname.str().c_str() );
				mrcw.write( cast->GetOutput() );

				// compute and save image variance
				if ( compute_variance_map == true ){
					metric.getImageVariance( globalAverage );
					globalAverage.getImageVariance( data );

					cast->SetInput( data );
					cast->Update();

					stringstream fnameVariance;
					if ( iter == 0 ){
						fnameVariance << pattern << "_global_variance.mrc";
					} else {
						fnameVariance << pattern << "_global_variance_" << iter << ".mrc";
					}
					mrcw.setFileName( fnameVariance.str().c_str() );
					mrcw.write( cast->GetOutput() );

				}
				// put average back into 'data' for symmetrization
				globalAverage.getImage( data );

				// symmetrize if needed
				if ( ( apply_rotational_symmetry == true ) || ( symmetry > 1 ) ){

					if ( apply_rotational_symmetry == true ){
						// compute rotationally symmetric average

						// apply rotational symmetry
						Array< double, 3 > in, out, res;
						Array< bool, 3 > zone;
						nbfVTKInterface :: vtkToBlitzReference( data, in );
						res.resize( in.shape() );
						// res = in.transpose(thirdDim,secondDim,firstDim);
						res = in;

						nbfCylindricalDomain3< double > cyl;
						TinyVector< int, 3 > center = in.shape() / 2;
						cyl.setCenter( center );
						cyl.setMaxRho( in.extent(firstDim) / 2 );
						cyl.setResRho( in.extent(firstDim) / 2 );
						cyl.setResTheta(180);
						cyl.cartesian2cylindrical( res, out, zone );

						// compute radial profile
						Array< double, 2 > profile( out.rows(), out.depth() );
						firstIndex i; secondIndex j; thirdIndex k;
						profile = sum( out(i,k,j), k );
						// the center axis does not get averaging, so we average with the closest ring outside.
						// profile( 0, Range :: all() ) = ( profile( 0, Range :: all() ) + profile( 1, Range :: all() ) ) / 2.0;
						for ( int i = 0; i < out.cols(); i++ ){
							out( Range :: all(), i, Range :: all() ) = profile;
						}

						// back-transform to cartesian coordinates
						cyl.cylindrical2cartesian( out, res, zone );

						// set symmetrized image as current average
						globalAverage.setAverageImage( res );

						imfilter.execute( res );
						fffilter.executeHalf( res );
						nbfVTKInterface :: blitzToVtk( res, data );
					} else {
						cout << "Applying " << symmetry << "-fold symmetry." << endl;
						cluster.symmetrize( globalAverage, symmetry );
						globalAverage.getImage( data );

						imfilter.fileMaskOn();
						imfilter.execute( data );
						//Array< double, 3 > res;
						//nbfVTKInterface :: vtkToBlitzReference( data, res );
						//fffilter.executeHalf( res );
					}

					Array< double, 3 > A;
					nbfVTKInterface :: vtkToBlitzReference( data, A );
					A *= -1;

					cast->SetInput( data );
					cast->Update();

					stringstream fnamesym;
					if ( iter == 0 ){
						fnamesym << pattern << "_global_average_symmetrized.mrc";
					} else {
						fnamesym << pattern << "_global_average_" << iter << "_symmetrized.mrc";
					}
					mrcw.setFileName( fnamesym.str().c_str() );
					mrcw.write( cast->GetOutput() );
				}

				if ( iter < maxIter ){

					// reset original alignments
					if ( resetAlignments > 0 ){
						for ( int i = 0; i < volumeList.size(); i++ ){
							volumeList[i].setTransform( (vtkTransform*)NULL );
						}
					}

					vector< nbfWedgedAverageImage3D< PIXEL > > reference;
					reference.push_back( globalAverage );
					alignments.resize( 1, volumeList.size(), 17, 1 );
					alignments = -1;
					// align with shift only (=2)
					// align with rotation but no refinement (=3)
					//nbfTimer t;
					//t.start();
					metric.getDistances( reference, volumeList, alignments, alignment_mode );

					cout << "Current score = " << sum( alignments( Range :: all(), Range :: all(), 0, Range :: all() ) ) / sum( globalAlignments( Range :: all(), 0, 0 ) ) << endl;
					//cout << "Current distances = " << alignments( Range :: all(), Range :: all(), 0, Range :: all() ) << endl;
					//t.stop();
					//cout << t.elapsedSeconds() << endl;

					vector< PIXEL > distances;
					
					// apply alignments
					for ( int i = 0; i < alignments.cols(); i++ ){

						// retrieve original transformation
						for ( int j = 0; j < 16; j++ ){
							matrix[j] = globalAverage.multipleAlignments( i, j, 0 );
						}
						vtkMatrix4x4 * matrix1 = vtkMatrix4x4 :: New();
						matrix1->DeepCopy( matrix );
					
						// retrieve incremental alignment
						for ( int j = 0; j < 16; j++ ){
							matrix[j] = alignments( 0, i, 1 + j, 0 );
						}
						vtkMatrix4x4 * matriz = vtkMatrix4x4::New();
						matriz->DeepCopy( matrix );
						
						vtkMatrix4x4 * result = vtkMatrix4x4 :: New();
						vtkMatrix4x4::Multiply4x4( matrix1, matriz, result );
						
						if ( resetAlignments > 0 ){
							volumeList[i].setTransform( matriz );
						} else {
							volumeList[i].setTransform( result );
						}

						result->Delete();
						matriz->Delete();
						matrix1->Delete();

						// set current distance value to global reference
						volumeList[i].setCutOffset( alignments( 0, i, 0, 0 ) );
					
						distances.push_back( alignments( 0, i, 0, 0 ) ); 
					}

					stringstream prealignedVolumeFile;
					prealignedVolumeFile << pattern << "_volumes_pre_centered_" << iter + 1 << ".txt";
					nbfWedgedSubImage3D< PIXEL > :: write( prealignedVolumeFile.str().c_str(), volumeList );

					// eliminate outliers: overwrite this->input with pruned list of volumes
					if ( ( distances.size() > 0 ) && ( top_percentage < 1 ) ){
						vector< PIXEL > sortedDistances( distances.size() );
						sortedDistances = distances;
						sort( sortedDistances.begin(), sortedDistances.end() );
						PIXEL cutoffDistance = sortedDistances[ top_percentage * ( sortedDistances.size() - 1 ) ];

						// transfer classification result to class attribute
						vector< nbfWedgedSubImage3D< PIXEL > > cleanVolumeList;
						for ( int i = 0; i < volumeList.size(); i++ ){
							if ( distances[i] < cutoffDistance ){
								cleanVolumeList.push_back( volumeList[i] );
							}
						}

						stringstream cleanedVolumeFile;
						cleanedVolumeFile << pattern << "_volumes_pre_centered_clean_" << iter + 1 << ".txt";
						nbfWedgedSubImage3D< PIXEL > :: write( cleanedVolumeFile.str().c_str(), cleanVolumeList );

						cleanVolumeList.clear();
					}

					// revert cut offset
					for ( int k = 0; k < volumeList.size(); k++ ){
						volumeList[k].setCutOffset( volumeCutOffset );
					}

				}
				cast->Delete();
				data->Delete();
				}
		}

		// end MPI
		metric.finalizeMPI();

		cout << "Normal program termination." << endl;
	} else {
		metric.slaveMPI();
	}
	MPI_Finalize();
	return 0;
}
