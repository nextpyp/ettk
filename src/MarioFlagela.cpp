#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

//#include <MacrosFlujos.hh>
//#include <IO/FlujosIO.hh>
#include <nbfTimer.h>

#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <fm/nbfFastMarching3Dccc.h>

#include <bp/nbfFastGeodesicPath3D.h>
#include <bp/nbfGeodesicPath3D.h>

#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>
#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkImageCast.h>

#include <io/nbfVTKInterface.h>
#include <bs/nbfBordStrategyMirror.h>
#include <nbfVeselnessFilter.h>

#define PIXEL float
#define DIM 3

void main( int argv, char ** argc )
{
	// upsample and merge segmentation to original grid
#if 1
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	Array< short, DIM > input;

	// read original (full resolution) image
	//reader->SetFileName( argc[1] );
	//reader->Update();

	vtkImageData * image = vtkImageData::New();
	image->DeepCopy( reader->GetOutput() );

	// read outer membrane segmentation (binned and cropped)
	reader->SetFileName( argc[1] );
	reader->Update();

	vtkImageData * outer = vtkImageData::New();
	outer->DeepCopy( reader->GetOutput() );

	// upsample by factor = 2
	vtkImageResample * resample = vtkImageResample::New();
	resample->SetInput( outer );
	resample->SetDimensionality(3);
	resample->SetAxisMagnificationFactor(0,.5);
	resample->SetAxisMagnificationFactor(1,.5);
	resample->SetAxisMagnificationFactor(2,.5);
	resample->Update();
	outer->DeepCopy( resample->GetOutput() );

	// fix spacing
	outer->SetSpacing(1,1,1);

	nbfVTKInterface converter;

	//// convert original image to blitz
	//Array< short, 3 > imageBlitz;
	//converter.vtkToBlitz( image, imageBlitz );
	//cout << "imageBlitz = " << imageBlitz.shape() << endl;
	////cout << min(imageBlitz) << ", " << max(imageBlitz) << endl;

	//// convert upsampled segmentation to blitz
	//Array< PIXEL,3 > outerBlitz;
	//converter.vtkToBlitz( outer, outerBlitz );
	//cout << "outerBlitz = " << outerBlitz.shape() << endl;

	//// create storage for full size segmentation
	//Array< PIXEL,3 > outerFull( imageBlitz.shape() );
	//outerFull = 0;

	//// merge upsampled segmentation to enlarged full size segmentation
	//outerFull( Range(0,332), Range(0,198), Range(197,511) ) = outerBlitz;
	//cout << "outerFull = " << outerFull.shape() << endl;
	//
	//// mask input image so only background is visible (inft inside)
	//imageBlitz = where( outerFull < 7, imageBlitz, numeric_limits<short>::max() );

	//// convert back to blitz for writing result
	//vtkImageData * out = vtkImageData::New();
	//converter.blitzToVtk( imageBlitz, out );

	// write result into file
	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName( argc[2] );
	writer->SetInput( outer );
	writer->Write();

	return;
#elif 0
	
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName( argc[1] );
	reader->Update();

	// convert original image to blitz
	Array< PIXEL, 3 > V;

	nbfVTKInterface converter;
	converter.vtkToBlitz( reader->GetOutput(), V );

	BordStrategyMirrorSimple< PIXEL,3 > bsForV( V, 1 );
	bsForV.refresh();

	Array< PIXEL, 3 > W( V.shape() );
	for ( int i = 0; i < 50; i++ ){
		W = Laplacian3D(V);
		V = V + .1 * W;
		bsForV.refresh();
	}
	nbfVeselnessFilter< float, 3 > vesel( V );
	vesel.execute(W);

	PIXEL maxW = max(W);
	PIXEL minW = min(W);

	W = maxW - W;
	W = where( W < 748, 748, W );
	W = where( W > 751, 751, W );
	W = W - min(W);
	W = W / max(W) + 1e-5;
	
	converter.vtkToBlitz( reader->GetOutput(), V );
	W = where( V == numeric_limits< short > ::max(), numeric_limits< PIXEL > :: max(), W );
	cout << min(W) << ", " << max(W) << endl;

	// convert back to blitz for writing result
	vtkImageData * out = vtkImageData::New();
	converter.blitzToVtk( W, out );

	// write result into file
	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName( argc[2] );
	writer->SetInput( out );
	writer->Write();

#elif 0

	// load masked background image
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName( argc[1] );
	reader->Update();

	// convert to blitz format
	Array< PIXEL, 3 > image;
	nbfVTKInterface converter;
	converter.vtkToBlitz( reader->GetOutput(), image );
	cout << min(image) << ", " << max(image) << endl;
	
	// reduce computational domain
	// image( Range(290,387), Range(54,154), Range(250,511) ) = numeric_limits< PIXEL >::max();

	// adjust contrast
	//image = where( image < 748, 748, image );
	//image = where( ( image > 751 ) && ( image < numeric_limits<short>::max() ), 751, image );
	//image = where( image == numeric_limits< PIXEL > ::max(), numeric_limits< PIXEL > :: max(), image );

	//image( Range(0,289), Range::all(), Range::all() ) = numeric_limits< PIXEL >::max();
	//image( Range(388,toEnd), Range::all(), Range::all() ) = numeric_limits< PIXEL >::max();
	//image( Range::all(), Range(0,53), Range::all() ) = numeric_limits< PIXEL >::max();
	//image( Range::all(), Range(155,toEnd), Range::all() ) = numeric_limits< PIXEL >::max();
	//image( Range::all(), Range::all(), Range(0,249) ) = numeric_limits< PIXEL >::max();
	//image( Range::all(), Range::all(), Range(0,460) ) = numeric_limits< PIXEL >::max();

	//Array< short, 3 > vImage( image( Range(290,387), Range(54,154), Range(250,511) ) );
	//Array< short, 3 > sImage( vImage.shape() );
	//sImage = vImage;

	//sImage = sImage - 660 + 1e-1;
	//sImage = sImage / max(sImage) + 1e-1;

	//cout << min(image) << ", " << max(image) << endl;

	//nbfFastMarching3Dccc< PIXEL > fastMarching( image );

	//int tSize = 11;
	//Array< PIXEL, 3 > T( tSize, tSize, tSize );
	//TinyVector< PIXEL, 3 > center = ( T.shape() - 1 ) / 2.0;
	//firstIndex i;
	//secondIndex j;
	//thirdIndex k;
	//T = sqrt( 1.0 * pow2( i - center[firstDim] ) + pow2( j - center[secondDim] ) + pow2( k - center[thirdDim] ) );
	//T = where( T < 5, 0 , 1 );

	//fastMarching.setTemplate( T );
	//Array< PIXEL, 3 > ccc;
	//cout << "applying ccc..." << endl;
	//fastMarching.applyTemplate(ccc);

	image( Range(0,55), Range(0,20), 60 ) = numeric_limits< PIXEL > :: max();
	image( Range(0,55), Range(30,toEnd), 60 ) = numeric_limits< PIXEL > :: max();
	image( Range(65,toEnd), Range(0,20), 60 ) = numeric_limits< PIXEL > :: max();
	image( Range(65,toEnd), Range(30,toEnd), 60 ) = numeric_limits< PIXEL > :: max();

	image( Range(0,13), 67, Range(0,15) ) = numeric_limits< PIXEL > :: max();
	image( Range(0,13), 67, Range(25,toEnd) ) = numeric_limits< PIXEL > :: max();
	image( Range(26,toEnd), 67, Range(0,15) ) = numeric_limits< PIXEL > :: max();
	image( Range(26,toEnd), 67, Range(25,toEnd) ) = numeric_limits< PIXEL > :: max();

	//nbfMatlabWriter mw;
	//mw.setFileName("blitz");
	//mw.write(ccc);

	//return;

	image = image + 1;

	nbfFastMarching3D< PIXEL > fm3( image );

	//TinyVector< int, DIM > min1(353,113, image.ubound(thirdDim) );
	TinyVector< int, DIM > min1(353-290,113-54, image.ubound(thirdDim));
	PIXEL ini = 0;

	
	fm3.setAliveSet( min1, ini );

	//TinyVector<int,3> stop(295,129,293);
	TinyVector<int,3> stop(295-290,129+3-54,293-250);

	//fastMarching.setStopPoint( stop );

	Array< PIXEL, DIM > distance( image.shape() );

	nbfTimer tEvolution, tReinit;

	cout << "Running fast marching..." << endl;
	tEvolution.start();
	fm3.execute( distance );
	tEvolution.stop();
	cout << "Done." << endl;

	cout << "Tiempo total del algoritmo (secs)= " << tEvolution.elapsedSeconds() << endl;

	nbfGeodesicPath3D< PIXEL > geodesic( distance );
	vector< TinyVector< int, 3 > > path;
	geodesic.getPath( distance, stop, path );
	//geodesic.getPath( stop, path );

	cout << path.size() << endl;

	vector< PIXEL > dist;
	vector< TinyVector< int, 3 > > points;
	for ( int i = 0; i < path.size(); i++ ){
		dist.push_back(0);
		points.push_back( path[i] );
	}
	image = 1;
	nbfFastMarching3D< PIXEL > fm3d( image );
	fm3d.setAliveSet( points, dist );
	Array< PIXEL, 3 > ipath( image.shape() );
	fm3d.setStopDistance(20);
	TinyVector<int,3> last = fm3d.execute( ipath );

	ipath = where( ipath == numeric_limits< PIXEL > :: max(), 20, ipath );

	//TinyVector< int, 2 > fin( input.ubound(firstDim),input.ubound(secondDim));
	//vector< TinyVector< PIXEL, DIM > > points;
	//fastMarching.getPath( fin, points, input );

	//cout << points.size() << endl;

	 //convert back to blitz for writing result
	
	//Array< PIXEL, 3 > subDistance( image( Range(290,387), Range(54,154), Range(250,511) ) );
	//Array< PIXEL, 3 > subDistance( image( Range(290,387), Range(54,154), Range(461,511) ) );

	vtkImageData * flagela = vtkImageData::New();
	//PIXEL dMax = distance(last);
	//image = where( image < numeric_limits<PIXEL>::max(), image, 0 );
	//cout << min(image) << ", " << max(image) << endl;

	converter.blitzToVtk( ipath, flagela );

	// write result into file
	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName( argc[2] );
	writer->SetInput( flagela );
	writer->Write();

	return;

#elif 0

	// load masked background image
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName( argc[1] );
	reader->Update();

	// convert to blitz format
	Array< PIXEL, 3 > image;
	nbfVTKInterface converter;
	converter.vtkToBlitz( reader->GetOutput(), image );

	TinyVector<int,3> stop(295-290,129+8-54,293-2-250);

	nbfGeodesicPath3D< PIXEL > geodesic( image );
	vector< TinyVector< int, 3 > > path;
	geodesic.getPath( image, stop, path );
	//geodesic.getPath( stop, path );

	cout << path.size() << endl;

	vector< PIXEL > dist;
	vector< TinyVector< int, 3 > > points;
	for ( int i = 0; i < path.size(); i++ ){
		dist.push_back(0);
		points.push_back( path[i] );
	}
	image = 1;
	nbfFastMarching3D< PIXEL > fm3d( image );
	fm3d.setAliveSet( points, dist );
	Array< PIXEL, 3 > ipath( image.shape() );
	fm3d.setStopDistance(20);
	TinyVector<int,3> last = fm3d.execute( ipath );

	ipath = where( ipath == numeric_limits< PIXEL > :: max(), 20, ipath );

	vtkImageData * flagela = vtkImageData::New();
	converter.blitzToVtk( ipath, flagela );

	// write result into file
	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName( argc[2] );
	writer->SetInput( flagela );
	writer->Write();

#else
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName( argc[1] );
	reader->Update();

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetInput( reader->GetOutput() );
	cast->SetOutputScalarTypeToShort();
	cast->Update();

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName( argc[2] );
	writer->SetInput( cast->GetOutput() );
	writer->Write();

#endif
}
