#define NBF_PARALLEL_IMPLEMENTATION_MPI 1
#include "mpi.h"

#define BZ_GENERATE_GLOBAL_INSTANCES

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkMath.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageResample.h>
#include <vtkImageMedian3D.h>
#include <vtkImageGaussianSmooth.h>

#include <bs/nbfBordStrategyMirror.h>
//#include <nbfDifferentials.h>

#include <em/nbfRadonStructure.h>
#include <em/nbfReconstruction3D.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <io/nbfImageReader.h>
#include <io/nbfMrcReader.h>
#include <io/nbfMrcWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfVTKInterface.h>
#include <io/nbfBlitzWriter.h>
#include <nbfTimer.h>

#include <mxml.h>

#define PIXEL float

#include <random/discrete-uniform.h>

int main( int argc, char ** argv )
{
	//  Initialize MPI.
	MPI_Init ( &argc, &argv );

	int my_id, num_procs;

	//  Get this processes's rank.
	MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

	//  Find out how many processes are available.
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	nbfReconstruction3D< PIXEL > rec3D;

	if ( my_id == 0 ){ // Master process

		cout << argv[0] << "\n";
		cout << " Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";

		//if ( argc != 12 ){
		//	std::cerr << "Missing Parameters " << endl;
		//	std::cerr << "Usage: " << argv[0] << endl;
		//	std::cerr << "  1. Input tilt series (Mrc file)" << endl;
		//	std::cerr << "  2. Input geometry (Xml file)" << endl;
		//	std::cerr << "  3. Algorithm (art=0, sirt=1, iterative=2, new=3) " << endl;
		//	std::cerr << "  4. Iterations " << endl;
		//	std::cerr << "  5. Output reconstruction (rec file)" << endl;
		//	return 1;
		//}

		//int left = atof( argv[6] ) + 1;
		//int right = atof( argv[7] ) + 1;
		//int top = atof( argv[8] );
		//int bottom = atof( argv[9] );
		//int zheight = - atof( argv[10] );
		//int thickness = atof( argv[11] );

		//if ( argc != 6 ){
		//	std::cerr << "Missing Parameters " << endl;
		//	std::cerr << "Usage: " << argv[0] << endl;
		//	std::cerr << "  1. Input tilt series (Mrc file)" << endl;
		//	std::cerr << "  2. Input geometry (Xml file)" << endl;
		//	std::cerr << "  3. Algorithm (art=0, sirt=1, iterative=2, new=3) " << endl;
		//	std::cerr << "  4. Iterations " << endl;
		//	std::cerr << "  5. Output reconstruction (rec file)" << endl;
		//	return 1;
		//}

		//FILE *fp;
		//mxml_node_t *tree;

		//fp = fopen( argv[2], "r");
		//tree = mxmlLoadFile(NULL, fp, MXML_NO_CALLBACK);
		//fclose(fp);

		//mxml_node_t *node = mxmlFindElement(tree, tree, "left", NULL, NULL, MXML_DESCEND);
		//int left = atof( node->child->value.text.string ) + 1;
		//node = mxmlFindElement(tree, tree, "right", NULL, NULL, MXML_DESCEND);
		//int right = atof( node->child->value.text.string ) + 1;
		//node = mxmlFindElement(tree, tree, "islice", NULL, NULL, MXML_DESCEND);
		//int top = atof( node->child->value.text.string );
		//node = mxmlFindElement(tree, tree, "jslice", NULL, NULL, MXML_DESCEND);
		//int bottom = atof( node->child->value.text.string );
		//node = mxmlFindElement(tree, tree, "zheight", NULL, NULL, MXML_DESCEND);
		//int zheight = - atof( node->child->value.text.string );
		//node = mxmlFindElement(tree, tree, "thickness", NULL, NULL, MXML_DESCEND);
		//int thickness = atof( node->child->value.text.string );

		// using IMOD convention
		if ( argc != 13 ){
			std::cerr << "Missing Parameters " << endl;
			std::cerr << "Usage: " << argv[0] << endl;
			std::cerr << "   1. Input tilt series (mrc file)" << endl;
			std::cerr << "   2. Input tlt file (containing angles)" << endl;
			std::cerr << "   3. FULLIMAGE[0]" << endl;
			std::cerr << "   4. SHIFT[0]" << endl;
			std::cerr << "   5. SHIFT[1]" << endl;
			std::cerr << "   6. SLICE[0]" << endl;
			std::cerr << "   7. SLICE[1]" << endl;
			std::cerr << "   8. THICKNESS" << endl;
			std::cerr << "   9. WIDTH" << endl;
			std::cerr << "  10. Algorithm (art=0, sirt=1, iterative=2, new=3) " << endl;
			std::cerr << "  11. Iterations " << endl;
			std::cerr << "  12. Output reconstruction (rec file)" << endl;
			return 1;
		}

		// -FULLIMAGE 2048,2048
		// -TILTFILE SIVmne_20080701_L3_7.rawtlt
	    // -SHIFT 148.0,48.0 -SLICE 856,1335 -THICKNESS 480 -WIDTH 480
		int fullimage = atoi( argv[3] );
		int shift[2];
		shift[0] = atoi( argv[4] );
		shift[1] = atoi( argv[5] );
		int slice[2];
		slice[0] = atoi( argv[6] );
		slice[1] = atoi( argv[7] );
		int thickness = atoi( argv[8] );
		int width = atoi( argv[9] );

		int top = fullimage - slice[1]; // slice[0];
		int bottom = fullimage - slice[0]; // slice[1];
		int left = - shift[0] + fullimage / 2.0 - width / 2.0; // shift[0] - width / 2;
		int right = - shift[0] + fullimage / 2.0 + width / 2.0; // shift[0] + width / 2;
		int zheight = shift[1]; // shift[1];

		//cout << " left = " << left << endl;
		//cout << "right = " << right << endl;

		int algorithm = atoi( argv[10] );
		int iterations = atoi( argv[11] );

		//// read ROI using Inspect3D conventions
		//int left = atoi(argv[2]);
		//int right = atoi(argv[3]);
		//int top = atoi(argv[4]) - 1;
		//int bottom = atoi(argv[5]) - 1;
		//int zheight = atoi(argv[6]); 
		//int thickness = atoi(argv[7]);
		//int algorithm = atoi(argv[8]);
		//int iterations = atoi(argv[9]);

		Array< short, 3 > tiltSeries;
		nbfMrcReader reader;
		reader.setFileName( argv[1] );

		vtkImageData * data = vtkImageData::New();
		Array< PIXEL, 1 > angles;

		Array< PIXEL, 1 > means;

		// read angles from tlt file
		if ( true ){
			vector< float > angles_list;
			// cout << argv[2] << endl;
			std :: ifstream src( argv[2], ifstream::in );
			stringstream inputString1;
			src >> inputString1.rdbuf();
			src.close();
			while ( inputString1.good() ){
				float numberOfReferences;
				inputString1 >> numberOfReferences;
				if ( inputString1.good() ){
					angles_list.push_back( numberOfReferences );
				}
			}
			// reader.read( data );
			reader.read( data, angles, means );
			angles.resize( angles_list.size() );
			for ( int k = 0; k < angles.rows(); k++ ){
				angles(k)=angles_list[k];
			}
		} else {
			// read angles from MRC header
			reader.read( data, angles, means );
		}
		cout << "Angles = " << angles << endl;

		if ( data->GetNumberOfPoints() == 0 ){
			cerr << "Error reading input file. Check that file is legal MRC of type 1.\n";
			return 1;
		}

		//nbfVTKInterface::vtkToBlitzReference( data, tiltSeries );
		//nbfMrcWriter mrcw;
		//mrcw.setFileName( "tiltseries.mrc" );
		//mrcw.write( tiltSeries );

		nbfMatlabWriter w;
		w.setFileName("p.matlab");

		if ( data->GetPointData()->GetScalars()->GetDataType() == VTK_FLOAT ){
			Array< float, 3 > realTiltSeries;
			nbfVTKInterface::vtkToBlitzReference( data, realTiltSeries );
			realTiltSeries -= min(realTiltSeries);
			realTiltSeries /= max(realTiltSeries);
			realTiltSeries = - numeric_limits< short >::max() + realTiltSeries * 2 * numeric_limits< short >::max();
			tiltSeries.resize( realTiltSeries.shape() );
			tiltSeries = cast< short >( realTiltSeries );
		} else {
			nbfVTKInterface::vtkToBlitzReference( data, tiltSeries );
		}

		if ( angles.rows() != tiltSeries.depth() ){
			cerr << "ERROR - The number of angles in tlt file and mrc file do not match." << endl;
			exit(0);
		}

		// SLICE AVERAGING
		//Array< float, 3 > A( tiltSeries.cols(), tiltSeries.depth(), 1 );
		//A = 0;
		//for ( int i = 82; i <= 84; i++ ){
		//	A( Range::all(), Range::all(), 0 ) += cast<float>( tiltSeries( i, Range::all(), Range::all() ) );
		//}
		//A = A / 3.0;
		//Array< short, 3 > output( A.shape() );
		//output = cast<short>( A );
		//nbfVTKInterface::blitzToVtk( output, data );

		//// eliminate hot pixels
		//for ( int i = 0; i < angles.size(); i++ ){
		//	Array< PIXEL, 2 > A( crop( Range::all(), Range::all(), i ) );
		//	A = where( A > numeric_limits<short>::min() + 1.5 * means(i), numeric_limits<short>::min() + means(i), A );
		//}

		int zeroTilt = 0;
		for ( int i = 0; i < angles.numElements(); i++ ){
			if ( fabs( angles(i) ) < fabs( angles( zeroTilt ) ) ){
				zeroTilt = i;
			}
		}

		angles *= vtkMath::DegreesToRadians();

		// BD5707dset03.ali.ali
		//left = 25; right = 865; top = 417; bottom = top + ( right - left ) + 1; // virus2
		//left = 0; right = 511; top = 110; bottom = 111; // virus2
		//left = 380; right = 468; top = 136 - 1; bottom = 244 - 1; // virus2

		//int zheight = 30; 
		//int thickness = 140;

		// 060315b-SIV_hpr_second.ali
		//left = 317; right = 836; top = 760; bottom = top + ( right - left ) + 1; // virus2
		//left = 617; right = 1136; top = 1125; bottom = top + ( right - left ) + 1; // virus1
		//left = 800; right = 1319; top = 736; bottom = top + ( right - left ) + 1; // virus3

		// Tcellsiv19.ali
		//left = 639; right = 950; top = 406; bottom = top + ( right - left ) + 1;
		//left = 639; right = 950; top = 406; bottom = 730;
		//left = 639; right = 950; top = 406; bottom = 407;
		//left = 639; right = 700; top = 300; bottom = top + ( right - left ) + 1;

		// 060205c-SIV_red0_0.ali
		// left = 891; right = 1334; top = 937; bottom = 1361;
		// left = 694;  right = 1230; top = 654; bottom = 1142;


		// 060205d-SIV_red0_0.ali
		// left = 1040; right = 1582; top = 792; bottom = 1196;
		// left = 694;  right = 1230; top = 654; bottom = 1142;

		// GroEL-Nov-06_7_final.ali
		// left = 1418; right = 1519; top = 916; bottom = 1017; thickness = 99; zheight = 11;


		PIXEL axis = ( tiltSeries.rows() - 1 ) / 2.0;
		PIXEL distanceToCenter = ( right + left ) / 2.0 - axis;
		int patchSize = 2 * sqrt( pow2( ( right - left + 1.0 ) / 2.0 ) + pow2( thickness / 2.0 ) ) + 2;
		PIXEL centerX = ( right + left ) / 2.0 - axis;
		PIXEL centerY = - zheight;
		Array< short, 3 > crop( patchSize, bottom - top + 1, tiltSeries.depth() );
		for ( int i = 0; i < angles.numElements(); i++ ){
			PIXEL angle = angles(i) - angles(zeroTilt);
			PIXEL offset = floor( cos( angle ) * centerX - sin(angle) * centerY );
			Range I( ceil( axis + offset - patchSize / 2 ), ceil( axis + offset - patchSize / 2 ) + crop.rows() - 1 );
			//cout << angles(i) * vtkMath::RadiansToDegrees() << " - " << I << endl;
			crop( Range::all(), Range::all(), i ) = tiltSeries( I, Range( top, bottom ), i );
		}

		nbfMrcWriter mrcw;
		//mrcw.setFileName( "crop.mrc" );
		vtkImageData * cropVtk = vtkImageData::New();
		Array< float, 3 > cropFloat( crop.shape() );
		cropFloat = cast<float>(crop);
		//nbfVTKInterface::blitzToVtk( cropFloat, cropVtk );
		//mrcw.write( cropVtk );

		cout << "Projection data croped to = " << crop.shape() << endl;

		tiltSeries.free();
		data->Delete();

		// normalize each projection
		cropFloat = cast<float>(crop);
		for ( int i = 0; i < crop.depth(); i++ ){
			vector< float > vImage( crop.rows() * crop.cols() );
			int c = 0;
			for ( int a = 0; a < crop.rows(); a++ ){
				for ( int b = 0; b < crop.cols(); b++ ){
					vImage[c++] = cropFloat(a,b,i);
				}
			}
			sort( vImage.begin(), vImage.end() );
			float median = vImage[ ceil( vImage.size() / 2.0 ) ];
			c = 0;
			for ( int a = 0; a < crop.rows(); a++ ){
				for ( int b = 0; b < crop.cols(); b++ ){
					vImage[c++] = pow2( cropFloat(a,b,i) - median );
				}
			}
			sort( vImage.begin(), vImage.end() );
			float variance = sqrt( vImage[ ceil( vImage.size() / 2.0 ) ] );
			cropFloat( Range :: all(), Range :: all(), i ) = ( cropFloat( Range :: all(), Range :: all(), i ) - median ) / variance;

			//vtkImageGaussianSmooth * smooth = vtkImageGaussianSmooth :: New();
			//vtkImageData * cropVtk = vtkImageData :: New();
			//nbfVTKInterface::blitzToVtk( cropFloat, cropVtk );
			//smooth->SetInput( cropVtk );
			//smooth->SetDimensionality(2);
			//smooth->SetStandardDeviations(1.0,1.0);
			//smooth->SetRadiusFactors(3.0,3.0);
			//smooth->Update();
			//nbfVTKInterface::vtkToBlitz( smooth->GetOutput(), cropFloat );
			//smooth->Delete();
			//cropVtk->Delete();
		}

		//mrcw.setFileName( "crop_normalized.mrc" );
		//nbfVTKInterface::blitzToVtk( cropFloat, cropVtk );
		//mrcw.write( cropVtk );
		cropVtk->Delete();

		cropFloat = - numeric_limits< short > :: max() + ( cropFloat - min(cropFloat) ) / ( max(cropFloat) - min(cropFloat) ) * numeric_limits< short > :: max();
		crop = cast< short >( cropFloat );

		//Array< PIXEL, 2 > slice( crop.rows(), crop.depth() );
		//Array< PIXEL, 2 > image( crop.rows(), crop.rows() );
		//image = 0;

		// store 3D reconstruction
		Array< PIXEL, 3 > reconstruction( crop.rows(), crop.cols(), crop.rows() );
		reconstruction = 0;

		rec3D.setProjections( crop );
		rec3D.setAngles( angles );
		rec3D.setImage( reconstruction );

		crop.free();

		cout << "Starting reconstruction" << endl;
		if ( algorithm < 2 ){
			rec3D.sirtMPI( iterations, algorithm );
		}
		else if ( algorithm == 2 ){
			rec3D.iterativeRefinement( iterations, 3 );
		} else {
			rec3D.progresiveIterativeRefinement(iterations);
			reconstruction = 0;
			rec3D.art( 1 );
		}

		cout << "done." << endl;

		//cout << "error = " << max( abs( recon - reconstruction ) ) << endl;

		//cout << "Reconstruction = [" << min(reconstruction) << "," << max(reconstruction) << "]" << endl;
		// cout << reconstruction.shape() << endl;

		//// set background to image mean
		//PIXEL mean = sum( reconstruction ) / sum( where( reconstruction != 0, 1, 0 ) );
		//reconstruction = where( reconstruction == 0, mean, reconstruction );

		//// normalize and cast output to MRC type 1
		//reconstruction -= min(reconstruction);
		//reconstruction = reconstruction / max(reconstruction) * 2 * numeric_limits<short>::max() - numeric_limits<short>::max();

		//cout << "Reconstruction = [" << min(reconstruction) << "," << max(reconstruction) << "]" << endl;

		int sizeX = right - left;

		int lowerbound = floor( ( reconstruction.rows() - sizeX ) / 2.0 );
		Range I( lowerbound, lowerbound + sizeX - 1 );
		// cout << I << endl;

		int lower = 0;
		if ( ( reconstruction.depth() - thickness ) / 2.0 - floor( ( reconstruction.depth() - thickness ) / 2.0 ) < .5 )
			lower = floor( ( reconstruction.depth() - thickness ) / 2.0 );
		else
			lower = ceil( ( reconstruction.depth() - thickness ) / 2.0 );
		Range K( lower, lower + thickness - 1 );

		// compatibilize output geometry to match Inpect3D's
		Array< float, 3 > subImage( reconstruction( I, Range::all(), K ) );
		//reconstruction.resize( subImage.shape() );
		//reconstruction = subImage.reverse(secondDim);
		//subImage = reconstruction.reverse(thirdDim);
				
		subImage -= min(subImage);
		subImage = subImage / max(subImage) * 2 * numeric_limits<short>::max() - numeric_limits<short>::max();

		// cout << "Reconstruction = [" << min(subImage) << "," << max(subImage) << "]" << endl;

		Array< short, 3 > output( subImage.shape() );
		output = cast< short >( subImage );

		//vtkImageData * output = vtkImageData::New();
		//nbfVTKInterface::blitzToVtk( subImage, output );
		//vtkImageCast * cast = vtkImageCast::New();
		//cast->SetInput( output );
		//cast->SetOutputScalarTypeToShort();
		//cast->Update();

		nbfMrcWriter mrc;
		mrc.setFileName( argv[12] );
		mrc.write( output );

		//mrc.setFileName( argv[5] );
		//mrc.write( cast->GetOutput() );

		//cast->Delete();
		//output->Delete();

		//w.setFileName("rec.matlab");
		//w.write(subImage);

		rec3D.finalizeMPI();

	} else {
		rec3D.slaveMPI();
	}

	MPI_Finalize();

	return 0;
}