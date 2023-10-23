#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkImageData.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkTransform.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfBlitzReader.h>

#include <nbfTimer.h>
#include <em/nbfDensityLocalMaximaClustering.h>
#include <em/nbfAngularSearchMetricEM3D.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfSphericalImageMetric3D.h>
#include <em/nbfIterativeAveraging.h>

#include <em/nbfSubVolume.h>

#define PIXEL double

/// Compute distance matrix from volume list and associated geometry file.

void main( int argv, char ** argc )
{
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName(argc[1]);
	reader->Update();

	vector< nbfSubVolume< PIXEL > > subvolumes;

	string outFile( argc[2] );
	nbfSubVolume< PIXEL > :: read( outFile, subvolumes );

	cout << subvolumes.size() << endl;

	TinyVector< int, 3 > size = 64;
	for ( unsigned int i = 0; i < subvolumes.size(); i++ ){
		subvolumes[i].setVolume( reader->GetOutput() );
		subvolumes[i].setDimensions( size );
		subvolumes[i].setCutOffset( 32 );
	}

	nbfWedgeFilter3D< PIXEL > w;

	nbfProjectionRotationMetric3D< PIXEL > m(&w);
	
	nbfAngularSearchMetricEM3D< PIXEL > metrika(&m);
	metrika.setAngleSearchX(-2,2);
	metrika.setAngleSearchY(-2,2);

	Array< PIXEL, 2 > A(2,2);
	A = 1, 2,
		1, 1;
	cout << A << endl;

	metrika.setStrategy(A);

	subvolumes[0].setNormal( TinyVector< PIXEL,3>(0,0,1) );
	subvolumes[1] = subvolumes[0];
	PIXEL phi = 45 * vtkMath::DegreesToRadians();
	subvolumes[1].setRotation( TinyVector< PIXEL, 3 >( 3, -4, 27 ) );

	metrika.setInput1( &subvolumes[0] );
	metrika.setInput2( &subvolumes[1] );
	cout << "d = " << metrika.getDistance() << endl;
	cout << "d = " << metrika.getCorrelationPeak() << endl;

	return;

#if 0
	// retrieve geometry from file
	nbfBlitzReader mreader;
	mreader.setFileName( argc[1] );
	Array< float, 2 > G;
	mreader.read(G);

	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();

	vector< vtkImageData * > input( G.rows() );
	vector< vtkTransform * > transforms;
	vector< double > wedgeLs;
	vector< double > wedgeUs;

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToDouble();

	vtkImageChangeInformation * center;
	vtkImageReslice * reslice;

	for ( int i = 0; i < G.rows(); i++ ){
		stringstream fileName;
		if ( i < 9 ){
			fileName << argc[2] << "_0" << i + 1 << ".vtk";
		} else{
			fileName << argc[2] << "_" << i + 1 << ".vtk";
		}
		cout << fileName.str() << endl;
		reader->SetFileName( fileName.str().c_str() );
		reader->Update();

		input[i] = vtkImageData::New();
		cast->SetInput( reader->GetOutput() );
		cast->Update();
		input[i]->DeepCopy( cast->GetOutput() );

		wedgeLs.push_back( G(i,0) );
		wedgeUs.push_back( G(i,1) );

		vtkTransform * t = vtkTransform::New();
		t->RotateX( G(i,2) );
		t->RotateY( G(i,3) );
		t->RotateZ( G(i,4) );
		t->Translate( G(i,5), G(i,6), G(i,7) );
		transforms.push_back( t );
	}

	vtkImageData * image = vtkImageData::New();
	image->DeepCopy( reader->GetOutput() );

	nbfWedge< PIXEL > w;
	w.setGeometry( image );
	w.bandPassOn(.2,.5,.01);

	/// metric to be used
	nbfProjectionRotationMetric3D< PIXEL > metric(&w);

	nbfAngularSearchMetric3D< PIXEL > metricSearch( &metric );
	metricSearch.setAngleSearchX(-2,2);
	metricSearch.setAngleSearchY(-2,2);
	metricSearch.setAngleSearchZ(0,0);

	Array< PIXEL, 2 > A(2,2);
	A = 1, 2,
		1, 1;
	cout << A << endl;

	metricSearch.setStrategy(A);
	
	Array< PIXEL, 2 > D;

	nbfTimer t;
	t.start();
	metricSearch.getMatrix( input, wedgeLs, wedgeUs, transforms, D );
	t.stop();
	cout << "running time = " << t.elapsedSeconds() << endl;

	cout << D << endl;

	nbfMatlabWriter writer;
	writer.setFileName( argc[3] );
	writer.write(D);

	//nbfDensityLocalMaximaClustering< PIXEL > cluster;
	//cluster.setInput( input );
	//cluster.setMetric( &metric );
	//cluster.setDistanceRadius(.12);

	//vector< vector< int > > classes = cluster.execute();

	//cout << "clusters = " << classes.size() << endl;

	//for ( int i = 0; i < classes.size(); i++ ){
	//	cout << "class " << i << ", elements = ";
	//	for ( int j = 0; j < classes[i].size(); j++ ){
	//		cout << classes[i][j] << ", ";
	//	}
	//	cout << endl;
	//}
#endif
}