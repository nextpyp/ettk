#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>

#include <io/nbfMrcReader.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfImageWriter.h>
#include <io/nbfImageReader.h>
#include <io/nbfVTKInterface.h>

void main( int argv, char ** argc )
{
	//if ( argv != 3 ){
	//	cout << "Parameters:" << endl;
	//	cout << "\t 1 : input file" << endl;
	//	cout << "\t 2 : output file (.vtk)" << endl;
	//	exit(0);
	//}


	//nbfMrcReader reader;

	//// 1. V - input image
	//Array< float, 3 > V;

#if 0
	nbfMatlabReader mreader;
	mreader.setFileName(argc[1]);
	mreader.read(V);
#elif 0
	reader.setFileName( argc[1] );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}
#else
	//vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	//reader->SetFileName( argc[1] );
	//reader->Update();
	//vtkImageData * image = vtkImageData::New();
	//int dim[3];
	//reader->GetOutput()->GetDimensions(dim);
	//image->SetDimensions( dim );
	//image->GetPointData()->SetScalars( reader->GetOutput()->GetPointData()->GetScalars() );

	Timer t;

	Array< float, 3 > V;
	//t.start();
	//nbfVTKInterface::vtkToBlitz( reader->GetOutput(), V );
	//t.stop();
	//cout << t.elapsedSeconds() << endl;

	//cout << V.shape() << ", " << V(100,100,100) << endl;

	//vtkImageData * back = vtkImageData::New();

	//Range I(0,V.rows()-1,2);
	//Range J(0,V.cols()-1,2);
	//Range K(0,V.depth()-1,2);
	//Array< short, 3 > P( V );
	////P = V;


	//t.start();
	//nbfVTKInterface::blitzToVtk( P, back );
	//t.stop();
	//cout << t.elapsedSeconds() << endl;


#endif

	//vtkImageReader * r = vtkBMPReader::New();
	//r->SetFileName( argc[1] );
	//r->Update();
	//nbfVTKInterface::vtkToBlitz( r->GetOutput(), V );
	//cout << V.shape() << endl;
	//cout << min(V) << ", " << max(V) << endl;
	
	//vtkStructuredPointsWriter * w = vtkStructuredPointsWriter::New();
	//w->SetInput( r->GetOutput() );
	//w->SetFileName( argc[2] );
	//w->Write();
	//w->Delete();

	//r->Delete();

	nbfImageReader imr;
	imr.setFileName( argc[1] );
	imr.read( V );

	nbfImageWriter imw;
	imw.setFileName( argc[2] );
	//Array< float, 2 > A( V.rows(), V.depth() );
	//A = cast<float>( V( Range::all(), 83 , Range::all() ) );

	//A = A - min(A);
	//A = A / max(A) * 255;

	//Array< float, 2 > B( V.shape() );
	//B = cast<float>(V);

	imw.write( V );

	//nbfMatlabWriter writer;
	//writer.setFileName(argc[2]);
	//writer.write(V);
}