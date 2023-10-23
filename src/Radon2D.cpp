#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkMath.h>

#include <nbfRadon.h>
#include <io/nbfImageReader.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <nbfTimer.h>
#include <bs/nbfBordStrategyMirror.h>

#define PIXEL float

void main( int argv, char ** argc )
{
	Array< PIXEL, 2 > image;
	nbfMatlabReader reader;
	reader.setFileName( argc[1] );

	Array< PIXEL, 2 > targetRadon;

	reader.read( targetRadon );
	//image = image / 256;
	image.resize( targetRadon.rows(), targetRadon.rows() );

	//nbfMatlabReader mreader;
	//mreader.setFileName( argc[1] );
	//mreader.read( image );

	nbfRadon< PIXEL > radoner;
	radoner.setImage( image );

	Timer timer;

	image.resize(11000,11000);
	Range I(0,image.rows()-2);
	Array< PIXEL, 2 > indirection( image.shape() );
	timer.start();
	indirection(I,Range::all()) = image(I+1,Range::all()) - image(I,Range::all());
	timer.stop();
	cout << timer.elapsedSeconds() << endl;

	BordStrategyMirrorSimple< float, 2 > bsForS( image, 1 );
	bsForS.refresh();

	timer.start();
	indirection = forward11n(image,firstDim);
	timer.stop();
	cout << timer.elapsedSeconds() << endl;

	return;

	Array< PIXEL, 1 > angles( 181 );
	firstIndex i;
	angles = ( i - 90.0 ) * vtkMath::DegreesToRadians();
	radoner.setAngles( angles );
	
	timer.start();

	radoner.forwardProjection();

	timer.stop();

	cout << timer.elapsedSeconds() << endl;

	Array< PIXEL, 2 > tmpRadon;
	radoner.getRadon( tmpRadon );

	tmpRadon = targetRadon;

	cout << "[" << min(targetRadon) << ", " << max(targetRadon) << "]\n";
	Array< PIXEL, 2 > Phi( image.shape() );
	Phi = .1;
	timer.start();
	radoner.art( targetRadon, Phi, atoi(argc[3]) );
	timer.stop();

	cout << timer.elapsedSeconds() << endl;

	nbfMatlabWriter writer;
	writer.setFileName( argc[2] );
	writer.write( Phi );
}