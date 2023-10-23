#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMrcReader.h>
#include <io/nbfMrcWriter.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfVTKInterface.h>
#include <nbfMaximalFlowMW.h>
// #include <nbfDifferentials.h>

#include <vtkImageCast.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageMathematics.h>
#include <vtkMath.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>

#include <vtkPoints.h>
#include <vtkImageData.h>

#define PIXEL float

int main( int argc, char ** argv )
{

	if ( argc < 8 ){
		cout << "USAGE: virus_segment_membrane input.mrc iradius oradius weight iterations variances output" << endl;
		return 1;
	}

	float innerRadius = atof( argv[2] );
	float outerRadius = atof( argv[3] );
	float weight = atof(argv[4]);
	int iterations = atoi( argv[5] );
	float varFactor = atof( argv[6] );
	float smooth = 0*1.0;

	// read image volume
	nbfMrcReader mrc;
	mrc.setFileName( argv[1] );

	vtkImageData * data = vtkImageData :: New();
	mrc.read(data);

	Array< PIXEL, 3 > V;
	nbfVTKInterface :: vtkToBlitz( data, V );

	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");
	//w.write(V);

	// store result in blitz array
	Array< float, 3 > R;
	R.resize( V.shape() );

	TinyVector< float, 3 > center = R.shape() / 2.0;

	firstIndex i;
	secondIndex j;
	thirdIndex k;
	R = sqrt( pow2(i-center[0]+0.0) + pow2(j-center[1]+0.0) + pow2(k-center[2]+0.0) );

	if ( argc > 8 ){
		Array< PIXEL, 3 > saved( V.shape() );
		saved = where( ( R > innerRadius ) && ( R < outerRadius ), V, 0 );
		nbfMrcWriter writermw;
		writermw.setFileName( argv[8] );
		writermw.write( saved );
	}

	//// normalize
	//float m = mean(V);
	//float v = sum( pow2( V - m ) ) / V.size() / varFactor;
	//float lower = m - v;
	//float upper = m + v;
	//V = where( V < lower, lower, V );
	//V = where( V > upper, upper, V );
	//V = V - min(V);
	//V = V / max(V);
	//V = V + weight;

	// normalize within inner and outer ring
	float m = sum( where( ( R > innerRadius ) && ( R < outerRadius ), V, 0 ) ) / sum( where( ( R > innerRadius ) && ( R < outerRadius ), 1, 0 ) );
	vector< float > pixeles;
	for ( int i = 0; i < V.rows(); i++ ){
		for ( int j = 0; j < V.cols(); j++ ){
			for ( int k = 0; k < V.depth(); k++ ){
				if ( ( R(i,j,k) > innerRadius ) && ( R(i,j,k) < outerRadius ) ){
					pixeles.push_back( V(i,j,k) );
				}
			}
		}
	}
	sort( pixeles.begin(), pixeles.end() );
	m = pixeles[ floor( pixeles.size() / 2.0 ) ];

	float v = sum( where( ( R > innerRadius ) && ( R < outerRadius ), pow2( V - m ), 0 ) ) / sum( where( ( R > innerRadius ) && ( R < outerRadius ), 1, 0 ) ) / varFactor;
	float lower = m - v;
	float upper = m + v;
	V = where( V < lower, lower, V );
	V = where( V > upper, upper, V );
	V = V - min(V);
	V = V / max(V);
	V = V + weight;

	// set pressure volume
	Array< float, 3 > P( V.shape() );
	P = 0;

	// set outer sphere as park of SINK
	P = where( R > outerRadius, -1, P );

	// set volume walls as SINK
	P( 0, Range::all(), Range::all() ) = -1;
	P( P.ubound(0), Range::all(), Range::all() ) = -1;
	P( Range::all(), 0, Range::all() ) = -1;	
	P( Range::all(), P.ubound(1), Range::all() ) = -1;
	P( Range::all(), Range::all(), 0 ) = -1;
	P( Range::all(), Range::all(), P.ubound(2) ) = -1;
	
	// set source as interior sphere
	P = where( R < innerRadius, 1, P );

	nbfMaximalFlow< float > flow;

	nbfTimer t;
	t.start();
	flow.execute(P,V,iterations);

	t.stop();

	//Array< float, 3 > G( V.shape() );
	//BordStrategyMirrorDouble< float, 3 > bsForV( V, 1 );

	//// parameters
	//float timeStep = .25;	// evolution time step
	//int eiterations = 3;	// evolution iterations
	//for ( int i = 0; i < eiterations; i++ ){
	//	V = P;
	//	bsForV.refresh();
	//	P = P + timeStep * curvature3D(V);
	//}

	//P = P - mean(P);

	//// fix geometry for MRC writer
	//Array< PIXEL, 3 > B;
	//P.transposeSelf(thirdDim,secondDim,firstDim);
	//B.resize( P.shape() );
	//B = P.reverse(secondDim);
	//nbfVTKInterface::blitzToVtk( B, data );

	vtkImageGaussianSmooth * smoothf = vtkImageGaussianSmooth :: New();
	smoothf->SetDimensionality(3);
	smoothf->SetRadiusFactors( smooth, smooth, smooth );
	nbfVTKInterface :: blitzToVtk( P, data );
	smoothf->SetInput( data );
	smoothf->Update();
	nbfVTKInterface :: vtkToBlitzReference( smoothf->GetOutput(), P );

	nbfMrcWriter writerm;
	writerm.setFileName( argv[7] );
	writerm.write( P );
	cout << "File " << argv[7] << " written.\nElapsed time = " << t.elapsedSeconds() << " seconds." << endl;
	
	data->Delete();
	smoothf->Delete();
	return 0;
}