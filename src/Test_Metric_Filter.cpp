#define NBF_PARALLEL_IMPLEMENTATION_MPI 1

#define NBF_AVERAGE_IN_RECIPROCAL_SPACE

#define BZ_GENERATE_GLOBAL_INSTANCES

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

#include <mxml.h>

#define PIXEL float

int main( int argc, char ** argv )
{
	if ( argc < 3 ){
		cout << "Usage: TestMetricFilter config.xml volume_index (filename)" << endl;
		cout << "       If volume_index = -1 you must also specify the mrc file to be filtered." << endl;
		exit(0);
	}

	cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";

	// initialize XML file readout structures
	FILE *fp;
	mxml_node_t *tree;
	fp = fopen( argv[1], "r");
	tree = mxmlLoadFile(NULL, fp, MXML_NO_CALLBACK);
	fclose(fp);
	mxml_node_t * xmlnode;

	// runing mode
	xmlnode = mxmlFindElement(tree, tree, "mode", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field mode not present." << endl;
		return 1;
	}
	int running_mode = atoi( xmlnode->child->value.text.string );
	if ( abs(running_mode) > 3 ){
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

	switch ( abs(running_mode) ){
			case NBF_LOOP_CLUSTERING_CLASS:

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
				maskFileTh = atof(xmlnode->child->value.text.string);

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
		// the last parameters sets to use cylyndrical mask. Assume cylindrical mask if size in Z == 0.
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

	int volIndex =  atoi( argv[2] );

	xmlnode = mxmlFindElement(tree, tree, "volumes", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field volumes not present." << endl;
		return 1;
	}
	char * spikesFileName = xmlnode->child->value.text.string;

	vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
	nbfWedgedSubImage3D< PIXEL > :: read( spikesFileName, volumeList );

	// size to cutout volumes
	xmlnode = mxmlFindElement(tree, tree, "volume_size", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field volume_size not present." << endl;
		return 1;
	}
	int volumeSize = atoi( xmlnode->child->value.text.string );

	// position to center volumes
	xmlnode = mxmlFindElement(tree, tree, "volume_cut_offset", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field volume_cut_offset not present." << endl;
		return 1;
	}
	int volumeCutOffset = atoi( xmlnode->child->value.text.string );

	// membership threshold by distance
	xmlnode = mxmlFindElement(tree, tree, "class_binning_factor", NULL, NULL, MXML_DESCEND);
	if ( xmlnode == NULL ){
		cerr << "ERROR - Invalid configuration file. Field class_binning_factor not present." << endl;
		return 1;
	}
	PIXEL class_binning_factor = atof( xmlnode->child->value.text.string );

	if ( volIndex == - 1 ){
		if ( argc < 4 ){
			cout << "ERROR - You must specify .mrc file to process if index is -1." << endl;
			return 0;
		} else {
			nbfWedgedSubImage3D< PIXEL > reference;
			reference.setFileName( argv[3] );
			nbfMrcReader mr;
			mr.setFileName( argv[3] );
			vtkImageData * data = vtkImageData :: New();
			mr.read(data);
			int dims[3];
			data->GetDimensions(dims);
			TinyVector< PIXEL, 3 > tdims( dims[0] / 2, dims[1] / 2, dims[2] / 2 );
			reference.setPosition( tdims );
			volumeList.push_back( reference );
			volIndex = volumeList.size() - 1;
		}
	}

	if ( ( volIndex >= 0 ) && ( volIndex <= volumeList.size() ) ){

		cout << "Processing volume " << volumeList[volIndex].getFileName() << endl;

		// set volume size and cutting offset to center of subvolume
		TinyVector< int, 3 > tsize( volumeSize, volumeSize, volumeSize );
		
		//TinyVector< PIXEL, 3 > mag(.5,.5,.5);
		//volumeList[volIndex].setMagnification( mag );
		
		//volumeList[volIndex].setNormal( TinyVector<PIXEL,3>(0,0,1) );

		volumeList[volIndex].setDimensions( tsize );
		// PIXEL p = volumeList[volIndex].getCutOffset();
		volumeList[volIndex].setCutOffset( volumeCutOffset );

		vtkImageData * data = vtkImageData::New();

		volumeList[ volIndex ].getImage( data );

		Array< PIXEL, 3 > A, B;
		nbfVTKInterface::vtkToBlitz( data, A );

		Array< double, 3 > result;
		if ( symmetry > 1 ){

			// assume symmetry axis is Z direction
			for ( int i = 0; i < symmetry; i++ ){
				vtkTransform * t = vtkTransform :: New();
				t->RotateZ( i * 360.0 / symmetry );
				volumeList[ volIndex ].getImage( data, t );
				t->Delete();
				Array< double, 3 > current;
				nbfVTKInterface :: vtkToBlitzReference( data, current );
				if ( result.size() == 0 ){
					result.resize( current.shape() );
					result = 0;
				}
				result += current;
			}

			result /= symmetry;

			A = cast< PIXEL >(result);
		}

		//// APPLY MAGNIFICATION FACTOR
		//if ( class_binning_factor != 1 ){
		//	nbfVTKInterface::blitzToVtk( A, data );
		//	vtkImageResample * resample = vtkImageResample::New();
		//	resample->SetDimensionality(3);
		//	resample->SetInput( data );
		//	resample->SetAxisMagnificationFactor( 0, class_binning_factor );
		//	resample->SetAxisMagnificationFactor( 1, class_binning_factor );
		//	resample->SetAxisMagnificationFactor( 2, class_binning_factor );
		//	resample->Update();
		//	data->DeepCopy( resample->GetOutput() );
		//	resample->Delete();
		//	nbfVTKInterface::vtkToBlitz( data, A );
		//}

		cout << A.shape() << endl;
		Array< double, 3 > C( A.shape() );
		C = cast<double>(A);
		nbfVTKInterface::blitzToVtk( C, data );
		fffilter.execute( data );
		nbfVTKInterface::vtkToBlitz( data, A );
		C = cast<double>(A);
		imfilter.execute( C );
		A = cast< PIXEL >(C);

		//// change geometry to make it compatible with mrc file
		//A.transposeSelf(thirdDim,secondDim,firstDim);
		//B.resize( A.shape() );
		//B = A.reverse(secondDim);
		//nbfVTKInterface::blitzToVtk( B, data );

		//vtkImageCast * cast = vtkImageCast::New();
		//cast->SetOutputScalarTypeToFloat();
		//cast->SetInput( data );
		//cast->Update();

		nbfMrcWriter mrcw;
		stringstream newname;
		newname << volumeList[ volIndex ].getFileName();
		newname << ".filtered.mrc";
		mrcw.setFileName( newname.str().c_str() );
		//mrcw.write( cast->GetOutput() );
		mrcw.write( A, true );
		//cast->Delete();

		//nbfMatlabWriter w;
		//w.setFileName("p.matlab");
		//w.write(A);

		data->Delete();
	} else {
		cout << "Invalid index into volume list. Valid range = [0," << volumeList.size() - 1 << "]." << endl;
	}
	return 0;
}