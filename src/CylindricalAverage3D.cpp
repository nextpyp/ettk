#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkImageReader.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageReslice.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkImageNonMaximumSuppression.h>
#include <vtkImageMathematics.h>
#include <vtkCylindricalTransform.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <nbfCylindricalDomain3.h>
#include <bs/nbfBordStrategyMirror.h>

#include <nbfCutSubVolumes.h>

void main( int argv, char ** argc )
{
	nbfMatlabWriter mwriter;
	mwriter.setFileName("ipath");

	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName(argc[1]);
	reader->Update();

	Array< short, 3 > A;
	nbfVTKInterface converter;
	converter.vtkToBlitz( reader->GetOutput(), A );

	//A.resize(100,100,60);
	Array< float, 3 > Af( cast<float>(A) );

	TinyVector< int, 3 > center = A.shape() / 2;
	firstIndex i;
	secondIndex j;
	//for ( int k = 0; k < Af.depth(); k++ ){
	//	Af( Range::all(), Range::all(), k ) = sqrt( pow2(i-center[0]+0.0)+pow2(j-center[1]+0.0) +0.0 );
	//}

	//nbfPolarDomain< float, 2 > polar;
	//Array< bool, 2 > B2;
	//polar.setCenter( TinyVector<int,2>( center[0], center[1] ) );
	//for ( int k = 0; k < Af.depth(); k++ ){
	//	polar.cartesian2polar( Af( Range::all(), Range::all(), k ), C( Range::all(), Range::all(), k ), B2 );
	//}

	nbfCylindricalDomain3< float > cyl;
	cyl.setCenter( center );
	//cyl.setMinRho();
	cyl.setMaxRho( A.ubound(firstDim) / 2 );
	cyl.setResRho( A.ubound(firstDim) / 2 );
	cyl.setResTheta(180);

	Array< float, 3 > C;
	Array< bool, 3 > B;
	cyl.cartesian2cylindrical(Af,C,B);

	//C = where( B == true, C, 0 );

	Array< float, 2 > S( C.rows(), C.depth() );
	thirdIndex k;
	S = sum( C(i,k,j), k );

	//C = where( B == true, C, 0 );
	//cout << C.shape() << endl;

	mwriter.write(S);

	reader->Delete();
}