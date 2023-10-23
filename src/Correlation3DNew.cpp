#include "mpi.h"
#define NBF_VERBOSE

#define BZ_GENERATE_GLOBAL_INSTANCES

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

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
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkXMLPImageDataWriter.h>
#include <vtkPoints.h>
#include <vtkPolyVertex.h>
#include <vtkProbeFilter.h>
#include <vtkContourFilter.h>
// #include <vtkImplicitModeller.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
//#include <io/nbfBlitzWriter.h>
//#include <io/nbfBlitzReader.h>
#include <io/nbfMrcWriter.h>
//#include <nbfVeselnessFilter.h>
//#include <bs/nbfBordStrategyMirror.h>

#include <em/nbfFourierImageMetricCore.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfExtractPointsAndNormals3D.h>
#include <em/nbfTemplateSearchEM.h>
#include <em/nbfCorrelationImageMetric.h>

#include <em/nbfCutSubVolumes.h>

//#include <vtkFastMarching.h>

#define PIXEL float

int main( int argc, char ** argv )
{
	cout << " " << argv[0] << " compiled on " << __DATE__ << " at " << __TIME__ << ".\n\n";

	if ( argc != 18 ){
		std::cerr << "Automatic spike picking by template matching on viral surface.\n" << endl;
		std::cerr << "Usage: " << argv[0] << endl;
		std::cerr << "  1. virus tomogram volume (excluding .rec extension)" << endl;
		std::cerr << "  2. virus membrane volume, possibly binned (mrc file) " << endl;
		std::cerr << "  3. level to extract surface from virus membrane volume" << endl;
		std::cerr << "  4. surface binning factor wrt to virus volume (1,2,3,etc) " << endl;
		std::cerr << "  5. lower tilt range" << endl;
		std::cerr << "  6. upper tilt range" << endl;
		std::cerr << "  7. dimension in which to restrict the template search,\n     e.g. to ignore bottom and top of virus: 0(X), 1(Y) or 2(Z)" << endl;
		std::cerr << "  8. lowest slice to search (in dimension specified above)" << endl;
		std::cerr << "  9. highest slice to search (in dimension specified above)" << endl;
		std::cerr << " 10. ellipsoid template size in X (0 if using template from file)" << endl;
		std::cerr << " 11. ellipsoid template size in Y (filename if using external template)" << endl;
		std::cerr << " 12. ellipsoid template size in Z (threshold for density if using external template, 0 if no thresholding desired)" << endl;
		std::cerr << " 13. minimum spacing between adjacent correlation peaks" << endl;
		std::cerr << " 14. correlation threshold (peaks with lower correlation\n     values than threshold are discarded)" << endl;
		//std::cerr << " 15. maximum number of extracted spikes" << endl;
		std::cerr << " 16. file with extracted spikes (txt file) " << endl;
		std::cerr << " 17. correlation surface for visualization (xml file) "<< endl;
		std::cerr << " 18. spike locations for visualization (xml file) "<< endl;
		return 0;
	}

	char * tomogramFileName = argv[1];
	char * surfaceFileName = argv[2];
	PIXEL surfaceLevelSet = atof( argv[3] );
	PIXEL surfaceMagnification = atof( argv[4] );

	PIXEL lower_tilt_range = atof( argv[5] );
	PIXEL upper_tilt_range = atof( argv[6] );

	int limitSearchDimension = atoi( argv[7] );
	int limitSearchLower = atoi( argv[8] );
	int limitSearchUpper = atoi( argv[9] );
	int templateSizeX = atoi( argv[10] ) + 4; // template padding
	int templateSizeY = atoi( argv[11] ) + 4;
	int templateSizeZ = atoi( argv[12] ) + 4;
	PIXEL cccPeakSpacing = atof( argv[13] );
	PIXEL cccPeakThreshold = atof( argv[14] );
	int cccSpikeCountThreshold = numeric_limits<int>::max();//atof( argv[15] );
	char * spikesGeometryFileName = argv[15];
	char * cccFileName = argv[16];
	char * spikesFileName = argv[17];

	// if size of template is 0, assume template is given as mrc file.
	char * templateFromFile;
	PIXEL templateFromFileThreshold = 0;
	if ( templateSizeX == 4 ){
		templateFromFile = argv[11];
		templateFromFileThreshold = atof( argv[12] );
	}

	//// read 3D image data
	//vtkImageData * myImage = vtkImageData::New();

	stringstream tomogramFileNameString;
	tomogramFileNameString << tomogramFileName << ".rec";
	nbfMrcReader reader;
	//reader.setFileName( tomogramFileNameString.str().c_str() );
	//reader.read( myImage );

	//Array< float, 3 > A;
	//nbfVTKInterface::vtkToBlitz( myImage, A );
	////cout << min(A) << "," << max(A) << endl;
	//A.reverseSelf(secondDim);
	////A( Range(72,82), Range(91,101), Range(52,62) ) = 0;
	//nbfVTKInterface::blitzToVtk( A, myImage );

	//nbfMrcWriter mw;
	//mw.setFileName( "test.mrc" );
	//mw.write( myImage );
	//reader.read( myImage );
	////return;

	//vtkXMLPImageDataWriter * pwrite = vtkXMLPImageDataWriter :: New();
	//pwrite->SetInput( myImage );
	//pwrite->SetFileName("image.xml");
	//pwrite->Write();
	//pwrite->Delete();

	// read 3D implicit surface
	vtkImageData * mySurface = vtkImageData::New();
	reader.setFileName( surfaceFileName );
	reader.read( mySurface );

	//vtkXMLPImageDataWriter * pwrite1 = vtkXMLPImageDataWriter :: New();
	//pwrite1->SetInput( mySurface );
	//pwrite1->SetFileName("image_seg.xml");
	//pwrite1->Write();
	//pwrite1->Delete();

#if 1

	// generate SPHERICAL template volume with specified size
	Array< PIXEL, 3 > myCropBlitz, stalk, head;

	Array< PIXEL, 3 > rawTemplate;
	if ( templateSizeX == 4 ){
		nbfMrcReader re;
		re.setFileName( templateFromFile );
		vtkImageData * data = vtkImageData :: New();
		re.read( data );
		nbfVTKInterface :: vtkToBlitz( data, rawTemplate );

		// binarize template using threshold
		if ( templateFromFileThreshold != 0 ){
			rawTemplate = where( rawTemplate < templateFromFileThreshold, 0, 1 );
		}

		myCropBlitz.resize( rawTemplate.shape() + 4 );
		myCropBlitz = mean( rawTemplate );
		myCropBlitz( Range(2,myCropBlitz.rows()-3), Range(2,myCropBlitz.cols()-3), Range(2,myCropBlitz.depth()-3)) = rawTemplate;

		//// enlarge template size by 1,5 times
		//myCropBlitz.resize( 1.5 * rawTemplate.shape() );

		//// set background and bottom of spike
		//myCropBlitz = min( rawTemplate );
		//myCropBlitz( Range :: all(), Range :: all(), Range(fromStart,1) ) = max( rawTemplate );

		//// set template image in centered-bottom part of enlarged image
		//Range I( floor(rawTemplate.rows() / 4.0), floor(rawTemplate.rows() / 4.0) + rawTemplate.rows() - 1 );
		//Range J( floor(rawTemplate.cols() / 4.0), floor(rawTemplate.cols() / 4.0) + rawTemplate.cols() - 1 );
		//Range K( 2, 2 + rawTemplate.depth() - 1 );
		//myCropBlitz( I, J, K ) = rawTemplate;

		data->Delete();
		myCropBlitz *= -1;
	} else {
		// build spike template

		// make spike head using a sphere
		myCropBlitz.resize( templateSizeX, templateSizeY, templateSizeZ );
		firstIndex i; secondIndex j; thirdIndex k;
		myCropBlitz = sqrt( pow2( 1.0 * i - templateSizeX / 2.0 ) / 1.85 + pow2( 1.0 * j - templateSizeY / 2.0 ) / 1.85 + pow2( 1.0 * k - templateSizeZ / 2.0 ) ) / .85;
		//myCropBlitz = sqrt( pow2( 1.0 * i - templateSizeX / 2.0 ) + pow2( 1.0 * j - templateSizeY / 2.0 ) );

		// add the stalk
		myCropBlitz( Range :: all(), Range :: all(), Range(0,10) ) = sqrt( pow2( 1.0 * i - templateSizeX / 2.0 ) + pow2( 1.0 * j - templateSizeY / 2.0 ) ) * 2.25;
		myCropBlitz = where( myCropBlitz < 20, myCropBlitz, max(myCropBlitz) );

		//myCropBlitz( Range :: all(), Range :: all(), Range( 35, toEnd) ) = max(myCropBlitz);


		for ( int i = 0; i < templateSizeZ / 2; i++ ){
			myCropBlitz( Range :: all(), Range :: all(), i ) = where( myCropBlitz( Range :: all(), Range :: all(), i ) > myCropBlitz( Range :: all(), Range :: all(), 0 ), myCropBlitz( Range :: all(), Range :: all(), 0 ), myCropBlitz( Range :: all(), Range :: all(), i ) );
		}

		// add the membrane
		// myCropBlitz( Range :: all(), Range :: all(), Range(0,1) ) = 0;

		// smooth the binary mask
		nbfImageFilter< PIXEL, 3 > imFilter;
		imFilter.gaussianFilterOn(1.75);
		vtkImageData * data = vtkImageData :: New();
		nbfVTKInterface :: blitzToVtk( myCropBlitz, data );
		imFilter.execute( data );
		nbfVTKInterface :: vtkToBlitz( data, myCropBlitz );
		data->Delete();

		myCropBlitz = myCropBlitz - mean(myCropBlitz);
		myCropBlitz = myCropBlitz / sqrt( sum( myCropBlitz * myCropBlitz ) );

		stalk.resize( myCropBlitz.shape() );
		stalk = max(myCropBlitz);
		stalk( Range :: all(), Range :: all(), Range(0,10) ) = myCropBlitz( Range :: all(), Range :: all(), Range(0,9) );

		head.resize( myCropBlitz.shape() );
		head = max(myCropBlitz);
		head( Range :: all(), Range :: all(), Range(11,toEnd) ) = myCropBlitz( Range :: all(), Range :: all(), Range(10,toEnd) );

		rawTemplate.resize( myCropBlitz.shape() - 4 );

		//// add blank on top
		//Array< float, 3 > enlarged( myCropBlitz.shape() + 20 );
		//enlarged = max(myCropBlitz);
		//Range I( 10, 10 + myCropBlitz.ubound(firstDim) );
		//enlarged( I, I, Range( fromStart, myCropBlitz.ubound(thirdDim) ) ) = myCropBlitz;
		//myCropBlitz.resize( enlarged.shape() );
		//myCropBlitz = enlarged;

		//myCropBlitz = sqrt( pow2( 1.0 * i - size[0] / 2.0 ) + pow2( 1.0 * j - size[1] / 2.0 ) );
		//myCropBlitz = where( myCropBlitz < 10, lowTh, highTh );

	}
#ifdef WIN32	
	nbfMatlabWriter w;
	w.setFileName("p.matlab");
	w.write(myCropBlitz);
	w.write(stalk);
	w.write(head);
#endif
	//nbfMrcWriter mrcw;
	//mrcw.setFileName("tempalte.mrc");
	//mrcw.write(myCropBlitz);

	//////////////////

	// generate list of points + normals from input surface for constrained template search
	vtkPolyData * points = vtkPolyData::New();
	nbfExtractPointsAndNormals3D< PIXEL > extract;
	extract.setSurface( mySurface, surfaceLevelSet );
	extract.setMagnification( surfaceMagnification );
	extract.execute( points );

	//Array< float, 3 > A;
	//nbfVTKInterface::vtkToBlitz( myImage, A );
	//float mi = min(A(Range::all(),Range::all(),59));

	// fix the Y-dimension reversal for coordinates of mrc files
	vtkPolyData * nbfpoints = vtkPolyData :: New();
	nbfpoints->DeepCopy( points );
	for ( int i = 0; i < points->GetNumberOfPoints(); i++ ){
		double point[3];
		nbfpoints->GetPoints()->GetPoint( i, point);
		//cout << "p=" << point[0] << ",\t" << point[1] << ",\t"<< point[2] << endl;
		//point[1] = myImage->GetDimensions()[1] - 1 - point[1];
		//nbfpoints->GetPoints()->SetPoint( i, point );
		//TinyVector< int, 3 > t(floor(point[0]), floor(point[1]), floor(point[2]));
		//A( t ) = mi;
	}

	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");
	//w.write(A);

	cout << "Number of surface points where CCC will be evaluated = " << points->GetNumberOfPoints() << endl;

	// build Fourier filter for computing ccc-based metric
	nbfFourierFilter< PIXEL, 3 > filter;
	//filter.bandPassOn(.02,.01,.5,.01);

	// define object for representing cropped volumes
	nbfWedgedSubImage3D< PIXEL > cutVolume;
	cutVolume.setFileName( tomogramFileNameString.str().c_str() );
	//cutVolume.setVolume( myImage );
	//cutVolume.setCutOffset( - ( myCropBlitz.depth() - 4 ) / 2.0 );
	cutVolume.setCutOffset( - rawTemplate.depth() / 2.0 );
	//cutVolume.setCutOffset( 0.0 );
	// cutVolume.setCutOffset( - templateSizeZ / 2.0 );
	//PIXEL wedge = 60;
	cutVolume.getWedge()->set( lower_tilt_range, upper_tilt_range );

	// build metric to be used by template search
	nbfImageFilter< PIXEL, 3 > imFilter;
	//imFilter.setMaskSize( ( myCropBlitz.rows() - 4 ) / 2.0, ( myCropBlitz.cols() - 4 ) / 2.0, ( myCropBlitz.depth() - 4 ) / 2.0, 2, false );
	//imFilter.gaussianFilterOn(4.0);
	imFilter.gaussianFilterOn(8.0);
	// imFilter.squareMaskSize = 4;

 	nbfCorrelationImageMetric< PIXEL, 3 > metric( &imFilter, &filter );
	//nbfFourierImageMetricCore< PIXEL, 3 > metric( &imFilter, &filter );
	metric.setTranslationSearchRestriction( 0 );
	
	// template search using given metric
	nbfTemplateSearchEM< PIXEL > tSearch( &metric );
	tSearch.setTemplate( myCropBlitz );
	tSearch.setTemplates( stalk, head );
	tSearch.setMinimumSeparation( cccPeakSpacing );
	tSearch.setCorrelationTh( cccPeakThreshold );

	// set range in full resolution size
	//tSearch.setHeightRange(thirdDim,122,238); // 060205c
	//tSearch.setHeightRange(thirdDim,124,220); // 060205b.0
	//tSearch.setHeightRange(thirdDim,0,101); // 060205b.1

	//// limit search around central slice
	//TinyVector< int, 3 > sizes = surfaceMagnification * reader.getDims();
	//limitSearchLower = sizes[limitSearchDimension] / 2 - limitSearchLower;
	//limitSearchUpper = sizes[limitSearchDimension] / 2 + limitSearchUpper;

	tSearch.setHeightRange( limitSearchDimension, limitSearchLower, limitSearchUpper );
	
	vector< int > spikes;
	vector< TinyVector< PIXEL, 3 > > spikes_normals;
	tSearch.execute( nbfpoints, &cutVolume, spikes, spikes_normals );

	cout << "Detected spikes = " << spikes.size() << endl;

	// transfer correlation value to original VTK structure
	for ( int i = 0; i < points->GetNumberOfPoints(); i++ ){
		points->GetPointData()->GetScalars()->SetTuple1( i, nbfpoints->GetPointData()->GetScalars()->GetTuple1(i) );
	}

	vector< float > originalCCC;
	for ( int i = 0; i < spikes.size(); i++ ){
		originalCCC.push_back( nbfpoints->GetPointData()->GetScalars()->GetTuple1( spikes[i] ) );
	}

	// sort spikes by correlation to template
	vector< int > sortedSpikes;
	for ( int i = 0; i < originalCCC.size(); i++ ){
		// search for current maxima
		float maximo = - numeric_limits< float > :: max();
		int cmaximo = 0;
		for ( int j = 0; j < originalCCC.size(); j++ ){
			if ( originalCCC[j] > maximo ){
				maximo = originalCCC[j];
				cmaximo = j;
			}
		}
		// ignore in future searches
		originalCCC[cmaximo] = - numeric_limits< float > :: max();
		sortedSpikes.push_back( spikes[cmaximo] );
		if ( sortedSpikes.size() == cccSpikeCountThreshold ){
			break;
		}
	}

	for ( int i = 0; i < sortedSpikes.size(); i++ ){
		cout << nbfpoints->GetPointData()->GetScalars()->GetTuple1( sortedSpikes[i] ) << endl;
	}

	cutVolume.setCutOffset(0);
	//cutVolume.setCutOffset( - rawTemplate.depth() / 2.0 );
	
	// build structure for saving
	vector< nbfWedgedSubImage3D< PIXEL > > volumes;
	volumes.resize( sortedSpikes.size() );
	for ( int i = 0; i < sortedSpikes.size(); i++ ){
		volumes[i] = cutVolume;
		
		double tmp[3];
		nbfpoints->GetPoint( sortedSpikes[i], tmp );
		TinyVector< PIXEL, 3 > p( tmp[0], tmp[1], tmp[2] );

		// retrieve current normal
		double normal[3];	
		nbfpoints->GetPointData()->GetNormals()->GetTuple( sortedSpikes[i], normal );
		TinyVector< PIXEL, 3 > n( normal[0], normal[1], normal[2] );

		p = floor( p - cutVolume.getCutOffset() * n + .5 );

		volumes[i].setPosition( p );
		//nbfpoints->GetPointData()->GetNormals()->GetTuple( sortedSpikes[i], tmp );
		//volumes[i].setNormal( TinyVector< PIXEL, 3 >( tmp[0], tmp[1], tmp[2] ) );
		volumes[i].setNormal( spikes_normals[ sortedSpikes[i] ] );
		volumes[i].setCutOffset(0);
		volumes[i].setCutOffset( nbfpoints->GetPointData()->GetScalars()->GetTuple1( sortedSpikes[i] ) );
	}

	//cout << "volumes = " << volumes.size() << endl;

	cout << "Saving spike geometry to file " << spikesGeometryFileName << " ... ";
	nbfWedgedSubImage3D< PIXEL > :: write( spikesGeometryFileName, volumes );
	cout << "done." << endl;

	/// Compute dintrinsic distance matrix between spikes.
	//Array< float, 2 > D( spikes.size(), spikes.size() );
	//D = 0;

	// save surface geometry with correlation values
	cout << "Saving CCC geometry to file " << cccFileName << " ... ";
	vtkXMLPPolyDataWriter * writer = vtkXMLPPolyDataWriter::New();
	writer->SetInput( points );
	//writer->SetFileTypeToBinary();
	writer->SetFileName( cccFileName );
	writer->Write();
	cout << "done." << endl;


	reader.setFileName( tomogramFileNameString.str().c_str() );
	TinyVector< int, 3 > sizes = reader.getDims();

	// write chimera model file
	
	stringstream chimeraFile;
	chimeraFile << tomogramFileName << "_auto.cmm";
	ofstream modelFile( chimeraFile.str().c_str(), ios::out );
	stringstream outputStream;
	outputStream << "<marker_sets>\n<marker_set name=\"spikes detected automatically\">\n";

	for ( int i = 0; i < sortedSpikes.size(); i++ ){
		double tmp[3];
		nbfpoints->GetPoint( sortedSpikes[i], tmp );

		double normal[3];	
		nbfpoints->GetPointData()->GetNormals()->GetTuple( sortedSpikes[i], normal );

		tmp[0] = floor( tmp[0] - cutVolume.getCutOffset() * normal[0] + .5 );
		tmp[1] = floor( tmp[1] - cutVolume.getCutOffset() * normal[1] + .5 );
		tmp[2] = floor( tmp[2] - cutVolume.getCutOffset() * normal[2] + .5 );

		outputStream << "<marker id=\"" << i + 1 << "\" x=\"" << tmp[0] / surfaceMagnification << "\" y=\"" << ( sizes[1] - 1 - tmp[1] ) / surfaceMagnification << "\" z=\"" << tmp[2] / surfaceMagnification << "\" r=\"1\" g=\"0.2\" b=\"1\" radius=\"2\"/>\n";
	}
	//cout << outputStream.str().c_str();
	outputStream << "</marker_set>\n</marker_sets>";
	modelFile << outputStream.str();
	modelFile.close();
	nbfpoints->Delete();


	//// set weights to uniform (w=1)
	//points->GetPointData()->GetScalars()->FillComponent(0,1);

	//for ( int i = 0; i < sortedSpikes.size(); i++ ){
	//	// Compute Distances
	//	vtkFastMarching * fm = vtkFastMarching :: New();
	//	fm->SetInput( points );
	//	fm->SetAlive( sortedSpikes[i] );
	//	fm->Update();
	//	for ( int j = i; j < sortedSpikes.size(); j++ ){
	//		D(i,j) = fm->GetOutput()->GetPointData()->GetScalars()->GetTuple1( sortedSpikes[j] );
	//	}
	//	fm->Delete();
	//}

	//// save result to matlab file
	//nbfMatlabWriter mw;
	//mw.setFileName( argc[5] );
	//mw.write(D);

	// save CCC function as surface in vtkPolyData
	vtkPolyData * poly = vtkPolyData::New();
	vtkPoints * pointsVtk = vtkPoints::New();
	vtkDoubleArray * normalsVtk = vtkDoubleArray::New();
	normalsVtk->SetNumberOfComponents(3);
	for ( int i = 0; i < sortedSpikes.size(); i++ ){
		pointsVtk->InsertNextPoint( points->GetPoint( sortedSpikes[i] ) );
		normalsVtk->InsertTuple( i, points->GetPointData()->GetNormals()->GetTuple( sortedSpikes[i] ) );
	}
	poly->SetPoints( pointsVtk );
	poly->GetPointData()->SetNormals( normalsVtk );

	cout << "Saving 3D positions in vtk file " << spikesFileName << " ... ";
	writer->SetInput( poly );
	writer->SetFileName( spikesFileName );
	//writer->SetFileTypeToBinary();
	writer->Write();
	cout << "done." << endl;

	writer->Delete();

	//return;
#elif 0

	nbfCutSubVolumes< float > cutter;
	cutter.setTemplate( myCropBlitz );
	cutter.setSurface( mySurface );
	cutter.setImage( myImage );
	cutter.setFullResolutionImage( myImageFullResolution );
	cutter.setTemplateFullResolution( myCropFullResolutionBlitz );
	cutter.setCorrelationTh(-1.0);
	cutter.setOffsetUnderSurface(0.0);
	cutter.setOffsetUnderSurfaceFullResolution(50.0);
	//cutter.setSurfaceLbound(-.5);
	//cutter.setSurfaceUbound(0);
	cutter.setSurfaceLbound(-0.05);
	cutter.setSurfaceUbound(0.05);
	cutter.setMinimumSeparation(20);
	//cutter.setHeightRange(thirdDim,55,125);
	cutter.setHeightRange(thirdDim,190/2,190/2+2);
	cutter.setNotToIgnoreValuesBelowSurface();

	vector< Array< short, 3 > > volumes;
	vector< vtkTransform * > transforms;
	cutter.execute( volumes, transforms );

	cout << volumes.size() << endl;
	
	//nbfMrcWriter mrwriter;
	//mrwriter.setFileName("test.mrc");
	//mrwriter.write( myImage );

#endif

	//vtkImageData * output = vtkImageData::New();

	//myImage->DeepCopy( cutVolume.getVolume() );

	//nbfMrcWriter writerm;
	//for ( int i = 1; i <= volumes.size(); i++ ){
	//	stringstream fileName;
	//	if ( i < 10 ){
	//		fileName << argv[8] << "_0" << i << ".mrc";
	//	}
	//	else{
	//		fileName << argv[8] << "_" << i << ".mrc";
	//	}
	//	writerm.setFileName( fileName.str().c_str() );
	//	//nbfVTKInterface converter;
	//	//cout << volumes[i-1].shape() << endl;
	//	//transforms[i]->Print(cout);
	//	//converter.blitzToVtk( volumes[i-1], output );
	//	vtkImageData * output = vtkImageData::New();
	//	volumes[i-1].setVolume( myImage);
	//	volumes[i-1].setDimensions( TinyVector< int, 3 >(56,56,56) );
	//	volumes[i-1].getSubVolume( output );
	//	//writer->SetInput( output );
	//	writerm.write( output );
	//	cout << "File: " << fileName.str() << " written." << endl;
	//	output->Delete();
	//}

	return 0;

	//// Store wedge limits and orientation
	///////////////////////////////////////
	//// Values are stored in a 2D blitz array.
	//// There is one row for each volume.
	//// Columns (8) are as follows:
	//// wedgeL wedgeU wedgeAngleX wedgeAngleY wedgeAngleZ positionX positionY positionZ

	//double lWedge = -59.99;
	//double uWedge = 59.99;
	//Array< double, 2 > save( transforms.size(), 8 );
	//save = 0;
	//for ( int i = 0; i < save.rows(); i++ ){
	//	save(i,0) = lWedge;
	//	save(i,1) = uWedge;
	//	double orient[3];
	//	transforms[i]->GetOrientation(orient);
	//	for ( int j = 0; j <= 2; j++ ){
	//		save(i,2+j) = orient[j];
	//	}
	//	double position[3];
	//	transforms[i]->GetPosition(position);
	//	for ( int j = 0; j <= 2; j++ ){
	//		save(i,5+j) = position[j];
	//	}
	//}
	//
	//nbfMatlabWriter bwriter;
	//bwriter.setFileName("geometry.matlab");
	//bwriter.write(save);

	//Array< float, 2 > K;
	//nbfBlitzReader breader;
	//breader.setFileName("geometry.blitz");
	//breader.read(K);
	//cout << K << endl;

	//writer->Delete();
	//myImage->Delete();
	//reader->Delete();
	//myImageFullResolution->Delete();
}
