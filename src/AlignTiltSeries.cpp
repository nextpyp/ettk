#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkImageData.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMrcWriter.h>

#include <fm/nbfFastMarching3D.h>
#include <bp/nbfFastGeodesicPath3D.h>
#include <bp/nbfGeodesicPath3D.h>

#include <vtkBMPWriter.h>
#include <vtkImageReslice.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>

#include <em/nbfCorrelationMetricFFTW2D.h>

#define PIXEL double

void main( int argv, char ** argc )
{
	vtkStructuredPointsReader * read = vtkStructuredPointsReader::New();
	read->SetFileName( argc[1] );
	read->Update();
	vtkImageData * data = vtkImageData::New();
	data->DeepCopy( read->GetOutput() );

#if 1
	Array< PIXEL, 3 > tiltSeries;
	nbfVTKInterface::vtkToBlitz( data, tiltSeries );
	
	int size = tiltSeries.rows();
	Array< PIXEL, 2 > A( size, size );
	Array< PIXEL, 2 > B( size, size );

	vtkImageData * input1 = vtkImageData::New();
	vtkImageData * input2 = vtkImageData::New();

	A = cast< PIXEL >( tiltSeries( Range(fromStart,size-1), Range(fromStart,size-1), 0 ) );
	nbfVTKInterface::blitzToVtk( A, input1 );

	nbfFourierFilter2D< PIXEL > filter;
	filter.setGeometry( input1 );
	filter.bandPassOn(.02,.5,.01);
	nbfCorrelationMetricFFTW2D< PIXEL > ccc( &filter );

	Array< PIXEL, 3 > ccc3D( A.rows(), A.cols(), tiltSeries.depth() - 1 );

	PIXEL tx = 0;
	PIXEL ty = 0;

	Array< PIXEL, 3 > alignedSeries( tiltSeries.shape() );
	alignedSeries( Range::all(), Range::all(), ( tiltSeries.depth() - 1 ) / 2 ) = tiltSeries( Range::all(), Range::all(), ( tiltSeries.depth() - 1 ) / 2 );

	for ( int i = ( tiltSeries.depth() - 1 ) / 2; i < tiltSeries.depth() - 1; i++ ){
		A = cast< PIXEL >( tiltSeries( Range::all(), Range::all(), i ) );
		nbfVTKInterface::blitzToVtk( A, input1 );
		B = cast< PIXEL >( tiltSeries( Range::all(), Range::all(), i + 1 ) );
		nbfVTKInterface::blitzToVtk( B, input2 );

		ccc.setInput1( input1 );
		ccc.setInput2( input2 );
		ccc.execute();
		double t[3];
		ccc.getTransform()->GetPosition(t);
		cout << t[0] << ", " << t[1] << ", " << t[2] << ", " << ccc.getCorrelationPeak() << endl;

		//tx += t[0];
		ty += t[1];

		// transform image if ccc succesful
		vtkImageReslice * reslice = vtkImageReslice::New();
		reslice->SetInterpolationModeToLinear();
		reslice->SetBackgroundLevel(mean(B));
		vtkImageData * aligned = vtkImageData::New();
		reslice->SetInput( input2 );

		vtkTransform * trans = vtkTransform::New();
		trans->Translate( -tx, -ty, 0 );
		reslice->SetResliceTransform( trans );
		reslice->Update();
		aligned->DeepCopy( reslice->GetOutput() );
		
		Array< PIXEL, 2 > Bi;
		nbfVTKInterface::vtkToBlitz( aligned, Bi );
		alignedSeries( Range::all(), Range::all(), i + 1 ) = Bi;

		//ccc.getCCCImage( B );
		//ccc3D( Range::all(), Range::all(), i ) = B;
	}

	ty = tx = 0;
	for ( int i = ( tiltSeries.depth() - 1 ) / 2; i > 0; i-- ){
		A = cast< PIXEL >( tiltSeries( Range::all(), Range::all(), i ) );
		nbfVTKInterface::blitzToVtk( A, input1 );
		B = cast< PIXEL >( tiltSeries( Range::all(), Range::all(), i - 1 ) );
		nbfVTKInterface::blitzToVtk( B, input2 );

		ccc.setInput1( input1 );
		ccc.setInput2( input2 );
		ccc.execute();
		double t[3];
		ccc.getTransform()->GetPosition(t);
		cout << t[0] << ", " << t[1] << ", " << t[2] << ", " << ccc.getCorrelationPeak() << endl;

		//tx += t[0];
		ty += t[1];

		// transform image if ccc succesful
		vtkImageReslice * reslice = vtkImageReslice::New();
		reslice->SetInterpolationModeToLinear();
		reslice->SetBackgroundLevel(mean(B));
		vtkImageData * aligned = vtkImageData::New();
		reslice->SetInput( input2 );

		vtkTransform * trans = vtkTransform::New();
		trans->Translate( -tx, -ty, 0 );
		reslice->SetResliceTransform( trans );
		reslice->Update();
		aligned->DeepCopy( reslice->GetOutput() );
		
		Array< PIXEL, 2 > Bi;
		nbfVTKInterface::vtkToBlitz( aligned, Bi );
		alignedSeries( Range::all(), Range::all(), i - 1 ) = Bi;

		//ccc.getCCCImage( B );
		//ccc3D( Range::all(), Range::all(), i ) = B;
	}


	//ccc3D = max(ccc3D) - ccc3D;
	//ccc3D = 10.0 * ccc3D / max(ccc3D) + 1e-3;
	//cout << "ccc = [" << min(ccc3D) << "," << max(ccc3D) << "]" << endl;

	TinyVector< int, 2 > init = minIndex( ccc3D( Range::all(), Range::all(), 0 ) );
	cout << init << endl;
	TinyVector< PIXEL, 3 > init3d( init[0], init[1], 0 );

	Array< PIXEL, 3 > distance( ccc3D.shape() );

	nbfFastMarching3D< PIXEL > fm( ccc3D );
	fm.setAliveSet( init3d, 0 );
	// fm.execute( distance );

	Array< short, 3 > sal( alignedSeries.shape() );
	sal = cast< short >( alignedSeries );

	vtkImageData * out = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( sal, out );
	nbfMrcWriter mrc;
	mrc.setFileName( argc[2] );
	mrc.write( out );

	return;

#else

	Array< PIXEL, 3 > distance;
	nbfVTKInterface::vtkToBlitz( data, distance );
	cout << "distance = " << distance.shape() << endl;
	cout << "[" << min(distance) << "," << max(distance) << "]" << endl;

	cout << distance( 157,156,95 ) << endl;

	Array< PIXEL, 3 > ccc3D;
	read->SetFileName( argc[2] );
	read->Update();
	vtkImageData * ccc = vtkImageData::New();
	ccc->DeepCopy( read->GetOutput() );
	nbfVTKInterface::vtkToBlitz( data, ccc3D );

#endif

	nbfGeodesicPath3D< PIXEL > path3D( distance );
	vector< TinyVector< int, 3 > > path;
	
	TinyVector< int, 2 > start2D = minIndex( ccc3D( Range::all(), Range::all(), ccc3D.ubound(thirdDim) ) );
	TinyVector< int, 3 > start( start2D[0], start2D[1], distance.ubound(thirdDim) );

	cout << distance( start ) << endl;

	path3D.getPath( distance, start, path );

	//vtkImageData * aligned = vtkImageData::New();
	//nbfVTKInterface::blitzToVtk( ccc3D, aligned );

	//vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	//writer->SetFileName( argc[2] );
	//writer->SetOutputTypeToBinary();
	//writer->SetInput( aligned );
	//writer->Write();
	//writer->Delete();

	//nbfMatlabWriter writer;
	//writer.setFileName( argc[3] );
	//writer.write( A );

}