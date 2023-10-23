#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkMath.h>

#include <bs/nbfBordStrategyMirror.h>
#include <nbfDifferentials.h>

#include <em/nbfRadonStructure.h>
#include <io/nbfImageReader.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <nbfTimer.h>

#define PIXEL float

#include <random/discrete-uniform.h>

void main( int argv, char ** argc )
{
	Array< PIXEL, 2 > projections;
	nbfMatlabReader reader;
	reader.setFileName( argc[1] );
	reader.read( projections );

	cout << projections.shape() << endl;
	Array< PIXEL, 2 > image( projections.cols(), projections.cols() );
	image = 0;

	Array< PIXEL, 1 > angles( projections.rows() );
	firstIndex i;
	PIXEL size = ( angles.numElements() - 1 ) / 2.0;
	angles = ( i - size ) / size * 60 * vtkMath::DegreesToRadians();
	//angles = 15 * vtkMath::DegreesToRadians();

	Array< PIXEL, 2 > tmp;
	reader.setFileName( argc[2] );
	reader.read( tmp );
	cout << tmp.shape() << endl;
	angles.resize( tmp.cols() );
	angles = tmp( 0, Range::all() ) * vtkMath::DegreesToRadians();

	cout << angles << endl;

	nbfRadonStructure< PIXEL > radoner;
	radoner.setImage( image );
	radoner.setProjections( projections );
	radoner.setAngles( angles );

	nbfTimer t;
	t.start();
	radoner.sirt(50);
	t.stop();
	cout << "t = " << t.elapsedSeconds() << endl;

	nbfMatlabWriter writer;
	writer.setFileName("tmp.array");
	writer.write(image);

	return;
}