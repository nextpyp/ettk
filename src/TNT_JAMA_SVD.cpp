#define NBF_PARALLEL_IMPLEMENTATION_MPI 1

#define NBF_AVERAGE_IN_RECIPROCAL_SPACE

#ifdef WIN32
	 #define NBF_VERBOSE
#endif

#include "mpi.h"

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <string.h>

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
#include <vtkImageContinuousDilate3D.h>
#include <vtkImageCast.h>
#include <vtkImageExtractComponents.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMrcWriter.h>
#include <bs/nbfBordStrategyMirror.h>

#include <em/nbfImageMetric.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfWedgedAverageImage3D.h>

#include <em/nbfUnsupervisedLoopClustering.h>

#include <mxml.h>

#include <../tnt/tnt.h>
#include <../tnt/jama_svd.h>

#include <nbfTimer.h>

#define PIXEL float

int main( int argc, char ** argv )
{
	nbfMatlabReader r;
	r.setFileName("pca_iteration_1_lowerDimensionalRepresentation.matlab");
	Array< PIXEL, 3 > Ag;
	r.read(Ag);
	Ag.transposeSelf(secondDim,firstDim,thirdDim);
	Array< PIXEL, 2 > A( 500, Ag.cols() );
	A = Ag( Range(0,499), Range :: all(), 0 );
	
	firstIndex i;
	secondIndex j;


	// compute A tilde
	Array< PIXEL, 1 > M( A.rows() );
	M = 1.0 / sqrt( sum( pow2( A(i,j) ), j ) );
	Array< PIXEL, 1 > N( A.cols() );
	N = 1.0 / sqrt( sum( pow2( A(j,i) ), j ) );
	for ( int i = 0; i < A.cols(); i++ ){
		A( Range :: all(), i ) *= sqrt(M);
	}
	for ( int i = 0; i < A.rows(); i++ ){
		A( i, Range :: all() ) *= sqrt(N);
	}
	nbfMatlabWriter w;
	w.setFileName("p.matlab");
	w.write(A);

	TNT :: Array2D< PIXEL > Atnt(A.rows(),A.cols(),A.data());
	nbfTimer t;
	t.start();
	JAMA :: SVD< PIXEL > csvd(Atnt);

	TNT :: Array2D< PIXEL > Stnt;
	csvd.getS(Stnt);
	Array< PIXEL, 2 > S( A.cols(), A.cols() );
	for ( int i = 0; i < A.cols(); i++ ){
		S(i,i) = Stnt[i][i];
	}
	w.write(S);

	TNT :: Array2D< PIXEL > Vtnt;
	csvd.getV(Vtnt);

	Array< PIXEL, 2 > V( A.cols(), A.cols() );
	for ( int i = 0; i < V.rows(); i++ ){
		for ( int j = 0; j < V.cols(); j++ ){
			V(i,j) = Vtnt[i][j];
		}
	}
	//for ( int i = 0; i < V.rows(); i++ ){
	//	V( Range :: all(), i ) /= sqrt(N);
	//}

	w.write(V);

	TNT :: Array2D< PIXEL > Utnt;
	csvd.getU(Utnt);

	Array< PIXEL, 2 > U( A.rows(), A.cols() );
	for ( int i = 0; i < U.rows(); i++ ){
		for ( int j = 0; j < U.cols(); j++ ){
			U(i,j) = Utnt[i][j];
		}
	}
	w.write(U);


	t.stop();
	cout << "t = " << t.elapsedSeconds() << " seconds." << endl;
	return 1;
}