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

	if ( argc != 2 ){
		cout << "USAGE: " << argv[0] << " mrcfile (without mrc extension) \n" << endl;
		return 0;
	}

	vtkImageData * data = vtkImageData :: New();

	stringstream imagefile;
	imagefile << argv[1] << ".mrc";

	nbfMrcReader reader;
	reader.setFileName( imagefile.str().c_str() );
	reader.read( data );

	Array< PIXEL, 3 > A;
	nbfVTKInterface :: vtkToBlitz( data, A );
	PIXEL m = mean(A);
	PIXEL v = sqrt( sum( pow2( A - m ) ) / A.size() );
	//PIXEL minimo = min(A);
	//PIXEL maximo = max(A);
	//PIXEL minimohalf = ( m + minimo ) / 2;
	//PIXEL maximohalf = ( m + maximo ) / 2;

	//A = ( A - m ) / v;
	//PIXEL cercano;

	////A = where( A < cercano, cercano, A );
	////A = where( A > cercano, cercano, A );
	float lower = m - 1.5 * v;
	float upper = m + 1.5 * v;
	//A = where( A < lower, lower, A );
	//A = where( A > upper, upper, A );

	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");
	//w.write(A);

	//cout << "A = " << min(A) << ", " << max(A) << endl;

	stringstream segfile;
	segfile << argv[1] << "_seg.mrc";

	reader.setFileName( segfile.str().c_str() );
	reader.read( data );
	
	Array< PIXEL, 3 > A_seg;
	nbfVTKInterface :: vtkToBlitz( data, A_seg );

	// check if all dimensions are equal
	if ( ( A.rows() != A.cols() ) || ( A.rows() != A.depth() ) ){
		int maxdim = A.rows();
		for ( int i = 1; i < 2; i++ ){
			if ( A.extent(i) > maxdim ){
				maxdim = A.extent(i);
			}
		}
		Array< PIXEL, 3 > oldA( A.shape() );
		oldA = A;
		A.resize( maxdim, maxdim, maxdim );
		A = m;
		TinyVector< int, 3 > diff = ( maxdim - oldA.shape() ) / 2;
		A( Range( diff(0), diff(0) + oldA.rows() - 1 ), Range( diff(1), diff(1) + oldA.cols() - 1 ), Range( diff(2), diff(2) + oldA.depth() - 1 ) ) = oldA;

		oldA = A_seg;
		A_seg.resize( maxdim, maxdim, maxdim );
		A_seg = min( oldA );
		A_seg( Range( diff(0), diff(0) + oldA.rows() - 1 ), Range( diff(1), diff(1) + oldA.cols() - 1 ), Range( diff(2), diff(2) + oldA.depth() - 1 ) ) = oldA;
	}

	Array< PIXEL, 1 > values(9);
	values = .1, .01, .0050, .0025, .0010, .00050, .00025, .0001, min(A_seg);

	cout << "Computing iso-surfaces at following thresholds: " << values << endl;

	Array< unsigned char, 2 > XYZall( A.rows() * values.size(), 3 * A.rows() );
	Array< bool, 2 > XYZsegmentedall( A.rows() * values.size(), 3 * A.rows() );

	vtkImageData * mySurface = vtkImageData::New();
	nbfVTKInterface :: blitzToVtk( A_seg, mySurface );

	PIXEL globalMinima = max(A);

	for ( int i = 0; i < values.size(); i++ ){

		// generate list of points + normals from input surface for constrained template search
		vtkPolyData * points = vtkPolyData::New();
		nbfExtractPointsAndNormals3D< PIXEL > extract;
		extract.setSurface( mySurface, values(i) );
		extract.setMagnification( 1 );

		Array< bool, 3 > Asegmented( A.shape() );
		Asegmented = false;

		if ( values(i) != min(A_seg) ){
			extract.execute( points );

			for ( int j = 0; j < points->GetNumberOfPoints(); j++ ){
				double point[3];
				points->GetPoints()->GetPoint( j, point );
				Asegmented( (int)point[0], (int)point[1], (int)point[2] ) = true;
			}
		}

		Array< bool, 2 > XYZsegmented( Asegmented.rows(), 3 * Asegmented.rows() );
		Array< bool, 2 > viewXsegmented( XYZsegmented( Range :: all(), Range( fromStart, A.rows() - 1 ) ) );
		Array< bool, 2 > viewYsegmented( XYZsegmented( Range :: all(), Range( A.rows(), 2 * A.rows() - 1 ) ) );
		Array< bool, 2 > viewZsegmented( XYZsegmented( Range :: all(), Range( 2 * A.rows(), toEnd ) ) );

		viewXsegmented = Asegmented( Asegmented.rows() / 2, Range :: all(), Range :: all() );
		viewYsegmented = Asegmented( Range :: all(), Asegmented.cols() / 2, Range :: all() );
		viewZsegmented = Asegmented( Range :: all(), Range :: all(), Asegmented.depth() / 2 );

		Array< PIXEL, 2 > XYZ( A.rows(), 3 * A.rows() );
		Array< PIXEL, 2 > viewX( XYZ( Range :: all(), Range( fromStart, A.rows() - 1 ) ) );
		Array< PIXEL, 2 > viewY( XYZ( Range :: all(), Range( A.rows(), 2 * A.rows() - 1 ) ) );
		Array< PIXEL, 2 > viewZ( XYZ( Range :: all(), Range( 2 * A.rows(), toEnd ) ) );

		//if ( i == values.ubound(firstDim) ){
			viewX = A( A.rows() / 2, Range :: all(), Range :: all() );
			// PIXEL ml = max( A( A.rows() / 2, Range :: all(), Range :: all() ) );
			// viewX = where( viewX == globalMinima, ml, viewX );
			// viewX = viewX - min(viewX);
			// viewX = viewX / max(viewX);
			NBFnormalize(viewX,2.5);

			viewY = A( Range :: all(), A.cols() / 2, Range :: all() );
			//ml = max( A( Range :: all(), A.cols() / 2, Range :: all() ) );
			//viewY = where( viewY == globalMinima, ml, viewY );
			//viewY = viewY - min(viewY);
			//viewY = viewY / max(viewY);
			NBFnormalize(viewY,2.5);

			viewZ = A( Range :: all(), Range :: all(), A.depth() / 2 );
			//ml = max( A( Range :: all(), Range :: all(), A.depth() / 2 ) );
			//viewZ = where( viewZ == globalMinima, ml, viewZ );
			//viewZ = viewZ - min(viewZ);
			//viewZ = viewZ / max(viewZ);
			NBFnormalize(viewZ,2.5);
			
		//} else {
		//	PIXEL epsilon = 1e-1;
		//	viewX = where( ( A_seg( A.rows() / 2, Range :: all(), Range :: all() ) < values(i) + epsilon ) && ( A_seg( A.rows() / 2, Range :: all(), Range :: all() ) > values(i) - epsilon ), m, A( A.rows() / 2, Range :: all(), Range :: all() ) );
		//	viewY = where( ( A_seg( Range :: all(), A.cols() / 2, Range :: all() ) < values(i) + epsilon ) && ( A_seg( Range :: all(), A.cols() / 2, Range :: all() ) > values(i) - epsilon ), m, A( Range :: all(), A.cols() / 2, Range :: all() ) );
		//	viewZ = where( ( A_seg( Range :: all(), Range :: all(), A.depth() / 2 ) < values(i) + epsilon ) && ( A_seg( Range :: all(), Range :: all(), A.depth() / 2 ) > values(i) - epsilon ), m, A( Range :: all(), Range :: all(), A.depth() / 2 ) );
		//}

		XYZ = XYZ - min(XYZ);
		XYZ = XYZ / max(XYZ);
		XYZ = XYZ * numeric_limits< unsigned char > :: max();

		Array< unsigned char, 2 > XYZuchar( XYZ.shape() );
		XYZuchar = cast< unsigned char >(XYZ);
		// XYZuchar.transposeSelf(1,0);

		XYZall( Range( i * A.rows(), (i+1)*A.rows()-1), Range :: all() ) = XYZuchar;
		XYZsegmentedall( Range( i * A.rows(), (i+1)*A.rows()-1), Range :: all() ) = XYZsegmented;

		points->Delete();
	}

	mySurface->Delete();

	vtkImageData * colorMontage = vtkImageData :: New();
	//nbfVTKInterface::blitzToVtk( XYZall, colorMontage );
	colorMontage->SetScalarTypeToUnsignedChar();
	colorMontage->SetNumberOfScalarComponents(3);
	colorMontage->SetDimensions( XYZall.rows(), XYZall.cols(), 1 );
	colorMontage->AllocateScalars();
	
	for ( int j = 0; j <  XYZall.cols(); j++ )
	{
		int jOffset = j * XYZall.rows();
		for ( int i = 0; i < XYZall.rows() ; i++ )
		{
			int offset = jOffset + i;
			if ( XYZsegmentedall( i, j ) == false ){
				colorMontage->GetPointData()->GetScalars()->InsertTuple3( offset, XYZall( i, j ), XYZall( i, j ), XYZall( i, j ) );
			} else {
				colorMontage->GetPointData()->GetScalars()->InsertTuple3( offset, numeric_limits< unsigned char > :: max(), numeric_limits< unsigned char > :: max(), 0 );
			}
		}
	}

	vtkPNGWriter * writer = vtkPNGWriter :: New();
	stringstream filename;
	filename << argv[1] << ".png";
	writer->SetFileName( filename.str().c_str() );
	writer->SetInput( colorMontage );
	writer->Write();
	colorMontage->Delete();
	writer->Delete();

	data->Delete();

	return 0;
}
