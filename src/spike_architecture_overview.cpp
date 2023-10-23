#define BZ_GENERATE_GLOBAL_INSTANCES

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
#include <vtkPolyData.h>
#include <vtkImageCast.h>
#include <vtkImageExtractComponents.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMODReader.h>
#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMrcReader.h>

#include <em/nbfImageMetric.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfWedgedAverageImage3D.h>

#include <em/nbfExtractPointsAndNormals3D.h>

#define PIXEL float

template< class Pixel >
void NBFnormalize( Array< Pixel, 2 > & A, float factor = 1.5 )
{
	PIXEL m = mean(A);
	PIXEL v = sqrt( sum( pow2( A - m ) ) / ( 1.0 * A.size() ) );
	float lower = m - factor * v;
	float upper = m + factor * v;
	A = where( A < lower, lower, A );
	A = where( A > upper, upper, A );
}


int main( int argc, char ** argv )
{
	cout << argv[0] << "\n";
	cout << " Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";

	if ( argc != 4 ){
		cout << "USAGE: " << argv[0] << " volumesFile gaussianDenoise output \n" << endl;
		return 0;
	}


	// read volume list
	vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
	nbfWedgedSubImage3D< PIXEL > :: read( argv[1], volumeList );

	if ( volumeList.size() == 0 ){
		cout << "Volume list is empty" << endl;
		return 1;
	}

	// initialize output png file (number of volumes x 3 othogonal views)
	int dimensions = volumeList[0].getDimensions()[0];
	Array< float, 2 > montage( dimensions * 3, dimensions * ( volumeList.size() + 0 ) );


	nbfImageFilter< PIXEL, 3 > imfilter;
	imfilter.gaussianFilterOn( atof(argv[2]) );

	vtkImageData * data = vtkImageData :: New();
	for ( int i = 0; i < volumeList.size(); i++ ){
		//TinyVector< PIXEL, 3 > norm(0,100,0);
		//volumeList[i].setNormal( norm );
		volumeList[i].setCutOffset( 0.0 );
		// read volume
		volumeList[i].getImage(data);

		// apply filter
		imfilter.execute( data );

		Array< PIXEL, 3 > A;
		nbfVTKInterface :: vtkToBlitz( data, A );

		PIXEL m = mean(A);
		PIXEL v = sqrt( sum( pow2( A - m ) ) / A.size() );
		float lower = m - 1.5 * v;
		float upper = m + 1.5 * v;

		// extract x,y,z slices
		//Ax.reverseSelf(secondDim);
		//Ay.transposeSelf(1,0);
		//Ay.reverseSelf(secondDim);
		Array< PIXEL, 2 > Ax( A( Range :: all(), A.cols() / 2, Range :: all() ) );
		Ax.reverseSelf(1);
		montage( Range( fromStart, A.rows() - 1 ), Range( i * A.rows(), ( i + 1 ) * A.rows() - 1 ) ) = Ax;
		Array< PIXEL, 2 > Ay( A( A.rows() / 2, Range :: all(), Range :: all() ) );
		Ay.reverseSelf(1);
		montage( Range( A.rows(), 2 * A.rows() - 1 ), Range( i * A.rows(), ( i + 1 ) * A.rows() - 1 ) ) = Ay;
		montage( Range( 2 * A.rows(), toEnd ), Range( i * A.rows(), ( i + 1 ) * A.rows() - 1 ) ) = A( Range :: all(), Range :: all(), A.depth() / 2 );
		Array< PIXEL, 2 > view( montage( Range::all(), Range( i * A.rows(), ( i + 1 ) * A.rows() - 1 ) ) );
		NBFnormalize(view,2.5);

		//if ( i == volumeList.size() - 1 ){
		//	montage( Range::all(), Range( ( i + 1 ) * A.rows(), ( i + 2 ) * A.rows() - 1 ) ) = 0;
		//}
	}
	montage.reverseSelf(1);
	data->Delete();

	montage = montage - min(montage);
	montage = montage / max(montage);
	montage = montage * numeric_limits< unsigned char > :: max();

	vtkImageData * colorMontage = vtkImageData :: New();
	colorMontage->SetScalarTypeToUnsignedChar();
	colorMontage->SetNumberOfScalarComponents(3);
	colorMontage->SetDimensions( montage.rows(), montage.cols(), 1 );
	colorMontage->AllocateScalars();
	
	for ( int j = 0; j <  montage.cols(); j++ )
	{
		int jOffset = j * montage.rows();
		for ( int i = 0; i < montage.rows() ; i++ )
		{
			int offset = jOffset + i;
			colorMontage->GetPointData()->GetScalars()->InsertTuple3( offset, montage( i, j ), montage( i, j ), montage( i, j ) );
		}
	}

	vtkPNGWriter * writer = vtkPNGWriter :: New();
	writer->SetFileName( argv[3] );
	writer->SetInput( colorMontage );
	writer->Write();
	colorMontage->Delete();
	writer->Delete();

	return 0;
}
