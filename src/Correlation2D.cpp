#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkImageData.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <vtkBMPWriter.h>
#include <vtkImageReslice.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkImageMedian3D.h>

#include <em/nbfCorrelationImageMetric2D.h>
#include <em/nbfCorrelationImageMetric3D.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfAngularSearchMetric3D.h>

#define PIXEL double

void main( int argv, char ** argc )
{
	//Array< PIXEL, 2 > P;
	//nbfMatlabReader mr;
	//mr.setFileName( argc[1] );
	//mr.read(P);

	//Array< PIXEL, 2 > A( P.shape() );
	//A = P;
	//Array< PIXEL, 2 > B( P.shape() );
	//B = P;

	vtkStructuredPointsReader * read = vtkStructuredPointsReader::New();
	read->SetFileName( argc[1] );
	read->Update();
	vtkImageData * data = vtkImageData::New();
	data->DeepCopy( read->GetOutput() );

	Array< PIXEL, 3 > tiltSeries;
	nbfVTKInterface::vtkToBlitz( data, tiltSeries );

	int size = min( tiltSeries.rows(), tiltSeries.cols() );

	Array< PIXEL, 2 > A( size, size );
	Array< PIXEL, 2 > B( size, size );

	A = cast< PIXEL >( tiltSeries( Range(fromStart,size-1), Range(fromStart,size-1), 10 ) );
	B = cast< PIXEL >( tiltSeries( Range(fromStart,size-1), Range(fromStart,size-1), 10 ) );

	vtkImageData * input1 = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( A, input1 );

	vtkImageData * input2 = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( B, input2 );

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetInput( data );
	cast->SetOutputScalarTypeToDouble();
	cast->Update();

	input1->DeepCopy( cast->GetOutput() );

	vtkImageMedian3D * median = vtkImageMedian3D::New();
	median->SetKernelSize(5,5,5);
	median->SetInput( input1 );
	median->Update();

	input1->DeepCopy( median->GetOutput() );


////////////////////////////////

	vtkImageChangeInformation * center = vtkImageChangeInformation::New();
	center->CenterImageOn();
	center->SetInput( input1 );
	center->Update();

	vtkImageReslice * reslice = vtkImageReslice::New();
	reslice->WrapOff();
	reslice->MirrorOff();
	reslice->SetBackgroundLevel( mean(tiltSeries) );
	reslice->SetInput( center->GetOutput() );
	reslice->SetInterpolationModeToCubic();
	
	vtkTransform * transform = vtkTransform::New();
	transform->RotateX( 0 );
	transform->RotateY( 0 );
	transform->RotateZ( 50 );
	transform->Translate(0,0,0);
	reslice->SetResliceTransform( transform );
	reslice->Update();

	input2->DeepCopy( reslice->GetOutput() );

	read->SetFileName( argc[2] );
	read->Update();
	input2->DeepCopy( read->GetOutput() );

	median->SetInput( input2 );
	median->Update();

	input2->DeepCopy( median->GetOutput() );

#if 1
//////////////////////////////

	nbfWedgeFilter3D< PIXEL > filter;
	filter.setDimensions( input1 );

	//filter.bandPassOn(.02,.5,.01);
	nbfProjectionRotationMetric3D< PIXEL > ccc( &filter );
	ccc.setInput1( input1 );
	ccc.setInput2( input2 );
	//ccc.edgeFilterOn();
	//ccc.windowOff();
	ccc.execute();
	double t[3];
	ccc.getTransform()->GetPosition(t);
	cout << t[0] << ", " << 2*t[1] << ", " << t[2] << ", " << ccc.getCorrelationPeak() << endl;
#else

	nbfWedgeFilter3D< PIXEL > filter;
	filter.setDimensions( input1 );
	filter.bandPassOn(.02,.5,.01);
	vtkTransform * t1 = vtkTransform::New();
	t1->RotateY(71.9768);
	t1->RotateZ(-5.79041);
	//filter.wedge1On( TinyVector<PIXEL,2>(-60,60), t1 );
	 //transform->Concatenate( t1 );
	//filter.wedge2On( TinyVector<PIXEL,2>(-60,60), transform );

	nbfCorrelationImageMetric3D< PIXEL > metric( &filter );

	nbfAngularSearchMetric3D< PIXEL > search( &metric );
	search.setInput1( input1 );
	search.setInput2( input2 );
	search.setAngleSearchZ(-180,180);

	Array< PIXEL, 2 > S(1,2);
	S = 1, 1;

	search.setStrategy(S);
	search.execute();
	double t[3];
	search.getTransform()->GetPosition(t);
	cout << t[0] << ", " << 2*t[1] << ", " << t[2] << ", " << search.getCorrelationPeak() << endl;

	//	// define wedge filter
	//nbfWedge< double > wedge;
	//wedge.setGeometry( image1 );
	//wedge.bandPassOn(.2,.5,.01);

	//nbfAngularSearchCorrelationMetricFFTW< double, 3 > metric( &wedge );
	//metric.setStrategy( strategy );
	//metric.setAngleSearchX(0,0);
	//metric.setAngleSearchY(0,0);
	//metric.setAngleSearchZ(0,0);
	//metric.setInput1( image1 );
	//metric.setInput2( image2 );
	//
	//metric.execute();
	//cout << metric.getCorrelationPeak() << endl;


#endif

	//vtkImageData * aligned = ccc.getAligned();
	//nbfVTKInterface::vtkToBlitz( aligned, A );

	//nbfMatlabWriter writer;
	//writer.setFileName( argc[3] );
	//writer.write( A );

/*
	nbfMatlabReader reader;

	// multi purpose writer
	nbfMatlabWriter writer;
	writer.setFileName( "ipath" );

	Timer t;

	// 1. V - input image
	Array< float, 3 > V;

	reader.setFileName( argc[1] );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	// 1. V - segmentation image
	Array< float, 3 > S;

	reader.setFileName( argc[1] );
	if ( reader.read( S ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	int xTemplate = 30;
	int yTemplate = 30;
	int zTemplate = 60;
	Array< float, 3 > gTemplate( xTemplate, yTemplate, zTemplate );
	firstIndex i; secondIndex j;  thirdIndex k;

	// theta runs in [0,2\pi]
	Array< float, 3 > theta( gTemplate.shape() );
	theta = 2.0 * j * vtkMath::Pi() / ( theta.ubound(secondDim) + 0.0 );
	
	// extent in the radial direction
	Array< float, 3 > rho( gTemplate.shape() );
	rho = i;

	// extent in the vertical direction
	Array< float, 3 > height( gTemplate.shape() );
	height = k;

	vector< TinyVector< int, 3 > > rootPoints;
	vector< TinyVector< float, 3 > > rootNormals;

	nbfLinearInterpolator<> interp( V );

	Array< float, 3 > X( gTemplate.shape() );
	Array< float, 3 > Y( gTemplate.shape() );
	Array< float, 3 > Z( gTemplate.shape() );

	Array< float, 3 > crop( gTemplate.shape() );
	Array< bool, 3 > B( gTemplate.shape() );

	Array< float, 3 > ccc( V.shape() );
	ccc = - numeric_limits< float > :: max();

	Array< float, 3 > :: iterator iter = S.begin();
	while ( iter != S.end() ){
		float slbound = 1;
		float subound = o;
		TinyVector< int, 3 > p = iter.position();
		// if close to the surface
		if ( ( (*iter) < subound ) && ( (*iter) > slbound ) ){
			rootPoints.push_back( iter.position() );
			float x = p(firstDim);
			float y = p(secondDim);
			float z = p(thirdDim);
			float dx = ( S(x+1,y,z) - S(x-1,y,z) ) / 2.0;
			float dy = ( S(x,y+1,z) - S(x,y-1,z) ) / 2.0;
			float dz = ( S(x,y,z+1) - S(x,y,z-1) ) / 2.0;
			rootNormals.push_back( TinyVector< float, 3 >( dx, dy, dz );
			X = x + ;
			Y = y + ;
			Z = z + ;
			interp.interpolate(X,Y,Z,crop,B);
			float mean = mean( crop );
			ccc( iter.position() ) = crop * gTemplate / pow2( crop - mean );
		}
		++iter;
	}

	// look for local maxima of correlation CCC
	iter = ccc.begin();
	vector< TinyVector< int, 3 > > :: iterator iterPoints = rootPoints.begin();
	while ( iterPoints != rootPoints.end() ){
		float x = (*iterPoints)(firstDim);
		float y = (*iterPoints)(secondDim);
		float z = (*iterPoints)(thirdDim);		
		// neighborhood size
		int m = 5;
		Range I( max( x - m, V.lbound(0) ), min( x + m, V.ubound(0) ) );
		Range I( max( y - m, V.lbound(1) ), min( y + m, V.ubound(1) ) );
		Range I( max( z - m, V.lbound(2) ), min( z + m, V.ubound(2) ) );
		float minCCCth = .5;
		if ( ( max( ccc(I,J,K) ) != (*iter) ) && ( ccc(x,y,z) > minCCCth ) ){
			rootPoints.erase(iter);
		}
		else{
			++iter;
		}
	}
	*/
}