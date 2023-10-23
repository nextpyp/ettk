#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

//#include <vtkMath.h>
//#include <vtkImageData.h>
//#include <vtkPointData.h>
//#include <vtkImageReader.h>
//#include <vtkStructuredPoints.h>
//#include <vtkStructuredPointsReader.h>
//#include <vtkStructuredPointsWriter.h>
//#include <vtkImageReslice.h>
//#include <vtkImageGaussianSmooth.h>
//#include <vtkImageChangeInformation.h>
//#include <vtkTransform.h>
//#include <vtkTransformFilter.h>
//#include <vtkImageNonMaximumSuppression.h>
//#include <vtkImageMathematics.h>
//#include <vtkImageContinuousDilate3D.h>

//#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <nbfVeselnessFilter.h>
#include <bs/nbfBordStrategyMirror.h>

#include <fm/nbfFastFastMarching2D.h>

void main( int argv, char ** argc )
{
	Array< float, 2 > I;
	nbfMatlabReader reader;
	reader.setFileName("slice.array");
	reader.read(I);

	I = I + 1e-2;

	cout << I.shape() << endl;

	vector< TinyVector< int, 2 > > points;
	points.push_back( TinyVector<int,2>(159/2,421/2) );
	points.push_back( TinyVector<int,2>(224/2,68/2) );
	points.push_back( TinyVector<int,2>(415/2,209/2) );
	points.push_back( TinyVector<int,2>(345/2,443/2) );
	points.push_back( TinyVector<int,2>(150,25) );
	points.push_back( TinyVector<int,2>(187,201) );
	points.push_back( TinyVector<int,2>(145,238) );

	vector< float > distances;
	distances.push_back(0);
	distances.push_back(0);
	distances.push_back(0);
	distances.push_back(0);
	distances.push_back(0);
	distances.push_back(0);
	distances.push_back(0);

	// compute centroid
	TinyVector< float, 2 > center(0,0);
	for ( int i = 0; i < points.size(); i++ ){
		center = center + 1.0 * points[i];
	}
	center = center / ( points.size() + 0.0 );

	//center[0] = 133;
	//center[1] = 133;

	cout << center << endl;

	// compute minimum circle radius
	float meanRho = numeric_limits<float > :: max();
	for ( int i = 0; i < points.size(); i++ ){
		float rho = sqrt( pow2( points[i][0] - center[0] ) + pow2( points[i][1] - center[1] ) );
		if ( rho < meanRho ){
			meanRho = rho;
		}
	}
	// meanRho = meanRho / ( 0.0 + points.size() );

	firstIndex i; 
	secondIndex j;

	Array< float, 2 > rho( I.shape() );
	rho = sqrt( pow2( i - center[firstDim] ) + pow2( j - center[secondDim] ) );

	//BordStrategyMirrorDouble< float, 2 > BSforI( I, 1 );
	//BSforI.refresh();
	//Array< float, 2 > Ix( I.shape() ), Iy( I.shape() );
	//Ix = - central12n(I,firstDim);
	//Iy = - central12n(I,secondDim);
	//Ix = sqrt( Ix * Ix + Iy * Iy );
	//Ix = Ix / max(Ix);

	//I = ( 1 - Ix ) + 1e-1;

	I = where( rho > 10, I / rho, numeric_limits< float > :: max() );

	nbfMatlabWriter mwriter;
	mwriter.setFileName("ipath");
	mwriter.write(I);

	nbfFastFastMarching2D< float > fm(I);
	fm.setAliveSet( points, distances );

	Array< float, 2 > D( I.shape() );
	fm.execute(D);

	BordStrategyMirrorDouble< float, 2 > BSforD( D, 1 );
	BSforD.refresh();

	Array< float, 2 > Dx( D.shape() ), Dy( D.shape() );
	Dx = - central12n(D,firstDim);
	Dy = - central12n(D,secondDim);
	D = sqrt( Dx * Dx + Dy * Dy );
	Dx = Dx / D;
	Dy = Dy / D;

	Array< float, 2 > phi( D.shape() );

	// initialize to mean circle
	phi = sqrt( pow2( i - center[0] ) + pow2( j - center[1] ) ) - meanRho * 2 / 3;

	BordStrategyMirrorDouble< float, 2 > BSforPhi( phi, 1 );
	BSforPhi.refresh();

	Array< float, 2 > R( I.shape() );

	int iterations = 2000;
	float timeStep = 0.25;
	for ( int i = 0; i < iterations; i++ ){

		R = curvature2D(phi);
		// R = where( R < 0, R, 0 );
		D = 1.0 * R - 1.0 * advection2D(phi,Dx,Dy);

		//cout << min(D) << ", " << max(D) << endl;

		phi = phi + timeStep * D;

		BSforPhi.refresh();

		//// reinitialization
		//if ( i % 10 == 0 ){
		//	int iterMorel = 10;

		//	// compute distances
		//	S = sussmanSign( phi );
		//	double stop = 1;
		//	do {
		//		R = morelSussman2D(phi,S);
		//		stop = abs( max(R) );
		//		phi = phi + 0.1 * R;
		//		BSforPhi.refresh();
		//		iterMorel--;
		//	} while( ( iterMorel > 0 ) & ( stop > 1 ) );

		//}
		cout << "iter = " << i << endl;
	}

	mwriter.write(phi);
}