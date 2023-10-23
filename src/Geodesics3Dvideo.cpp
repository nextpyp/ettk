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
	nbfMatlabReader reader;

	Timer t;

	// 1. V - input image
	Array< float, 3 > V;

	Array< float, 3 > implicit;

	nbfMinimalSurfaceVideo< float > ms;

	int movie = atoi( argc[1] );

	char * fileName, * fileNameStart, * fileNameEnd;
	char * resultFileName;
	switch ( movie ){
		case 1:
			fileName = "C:/cygwin/home/abarte/workspace/electron/data/video/fedex.array";
			fileNameStart = "C:/cygwin/home/abarte/workspace/electron/data/video/fedexStart.array";
			fileNameEnd = "C:/cygwin/home/abarte/workspace/electron/data/video/fedexEnd.array";
			resultFileName = "C:/cygwin/home/abarte/workspace/electron/data/video/fedex";
			ms.addPointToAxis( TinyVector<float,3>(29,32,0) );
			ms.addPointToAxis( TinyVector<float,3>(44,55,249) );
			ms.addPointToAxis( TinyVector<float,3>(51,81,279) );
			break;
		case 2:
			fileName = "C:/cygwin/home/abarte/workspace/electron/data/video/mini.array";
			fileNameStart = "C:/cygwin/home/abarte/workspace/electron/data/video/miniStart.array";
			fileNameEnd = "C:/cygwin/home/abarte/workspace/electron/data/video/miniEnd.array";
			resultFileName = "C:/cygwin/home/abarte/workspace/electron/data/video/mini";
			ms.addPointToAxis( TinyVector<float,3>(56,41,0) );
			ms.addPointToAxis( TinyVector<float,3>(56,85,29) );
			ms.addPointToAxis( TinyVector<float,3>(59,134,49) );
			ms.addPointToAxis( TinyVector<float,3>(65,209,69) );
			ms.addPointToAxis( TinyVector<float,3>(68,251,77) );
			break;
		case 3:
			fileName = "C:/home/project/HIV/data/video/liron.array";
			fileNameStart = "C:/home/project/HIV/data/video/lironStart.array";
			fileNameEnd = "C:/home/project/HIV/data/video/lironEnd.array";
			resultFileName = "C:/home/project/HIV/data/video/liron";
			ms.addPointToAxis( TinyVector<float,3>(253,551,0) );
			ms.addPointToAxis( TinyVector<float,3>(254,528,5) );
			ms.addPointToAxis( TinyVector<float,3>(254,507,11) );
			ms.addPointToAxis( TinyVector<float,3>(258,486,23) );
			ms.addPointToAxis( TinyVector<float,3>(258,481,27) );
			ms.addPointToAxis( TinyVector<float,3>(258,458,35) );
			ms.addPointToAxis( TinyVector<float,3>(264,449,39) );
			ms.addPointToAxis( TinyVector<float,3>(261,399,57) );
			ms.addPointToAxis( TinyVector<float,3>(266,366,66) );
			break;
	}

	cout << "Processing " << fileName << " ..." << endl;

	reader.setFileName( fileName );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	//nbfImageWriter imwriter;
	//imwriter.setFileName("liron.vtk");
	//imwriter.write(V);

	implicit.resize( V.shape() );

	reader.setFileName( fileNameStart );
	if ( reader.read( implicit( Range::all(), Range::all(), V.lbound(thirdDim) ) ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	reader.setFileName( fileNameEnd );
	if ( reader.read( implicit( Range::all(), Range::all(), V.ubound(thirdDim) ) ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	nbfMatlabWriter writer;
	writer.setFileName( resultFileName );

	t.start();
	ms.search(V,implicit);
	t.stop();
	cout << "t = " << t.elapsedSeconds() << endl;
	cout << "Saving result to " << resultFileName << " ..." << endl;
	writer.write(V);
}