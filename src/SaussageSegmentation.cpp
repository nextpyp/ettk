#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfImageWriter.h>
#include <vtkMath.h>
#include <nbfVeselnessFilter.h>
#include <nbfGaussianFilter.h>
#include <nbfPolarDomain.h>
#include <nbfMinimalSurfaceVideo.h>

void main( int argv, char ** argc )
{
	nbfMrcReader reader;

	Timer t;

	// 1. V - input image
	Array< float, 3 > V;

	Array< float, 3 > implicit;

	nbfMinimalSurfaceSaussage< float > ms;

	ms.addPointToAxis( TinyVector<float,3>(29,32,0) );
	ms.addPointToAxis( TinyVector<float,3>(44,55,249) );
	ms.addPointToAxis( TinyVector<float,3>(51,81,279) );

	reader.setFileName( fileName );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	implicit.resize( V.shape() );

	nbfMatlabWriter writer;
	writer.setFileName( resultFileName );

	t.start();
	ms.search(V,implicit);
	t.stop();
	cout << "t = " << t.elapsedSeconds() << endl;
	cout << "Saving result to " << resultFileName << " ..." << endl;
	writer.write(V);
}