#include "mpi.h"
#define NBF_VERBOSE

#define BZ_GENERATE_GLOBAL_INSTANCES

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
#include <vtkImageContinuousDilate3D.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkXMLPImageDataWriter.h>
#include <vtkPoints.h>
#include <vtkPolyVertex.h>
#include <vtkProbeFilter.h>
#include <vtkContourFilter.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMrcWriter.h>

#include <em/nbfFourierImageMetricCore.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfExtractPointsAndNormals3D.h>
#include <em/nbfTemplateSearchEM.h>
#include <em/nbfCorrelationImageMetric.h>

#include <em/nbfCutSubVolumes.h>


#define PIXEL double

int main( int argc, char ** argv )
{
#if 1
	nbfMrcReader reader;	
	vtkImageData * data = vtkImageData::New();
	
	reader.setFileName( argv[1] );
	reader.read(data);

	Array< float, 3 > series1, series2;
	nbfVTKInterface::vtkToBlitz(data,series1);

	printf("%s=[%f,%f]\n",argv[1],min(series1), max(series1) );
	cout << series1.shape() << endl;
	
	reader.setFileName( argv[2] );
	reader.read(data);
	nbfVTKInterface::vtkToBlitz(data,series2);
	printf("%s=[%f,%f]\n",argv[2],min(series2), max(series2) );
	cout << series2.shape() << endl;

	series1 -= series2;
	printf("%s=[%f,%f]\n",argv[1],min(series1), max(series1) );
	
	nbfMrcWriter writer;
	writer.setFileName( argv[3] );
	writer.write(series1);	
#else

	nbfMrcReader reader;
	
	vtkImageData * data = vtkImageData::New();
	
	Array< float, 3 > series, virion;
	
	reader.setFileName( argv[1] );
	Array< float, 1 > tilts, means;
	
	reader.read(data,tilts,means);

	ifstream input( argv[2], ios::in );

	stringstream tilt_file;
	input >> tilt_file.rdbuf();
	
	cout << tilts << endl;
	
	cout << "read done" << endl;
	nbfVTKInterface::vtkToBlitzReference(data,series);
	
	cout << "copy done" << endl;
	
	reader.setFileName( argv[3] );
	
	vtkImageData * data1 = vtkImageData::New();
	reader.read(data1);
	nbfVTKInterface::vtkToBlitzReference(data1,virion);

	//virion=virion*mean(series)/mean(virion);
	//cout << "Virion = " << virion.shape() << ", [" << min(virion) << "," << max(virion) << "]" << endl;
	
	Range J(2048-1671,2048-1192);

	float virion_x=1024+128;
	float virion_z=-20;

	float tilt_x, angle;

	Array< float, 3 > dif( 480, 480, series.depth() );
	
	for ( int i = 0; i < tilts.size(); i++ ){
	//for ( int i = 0; i < 1; i++ ){
		
		tilt_file >> angle;
		cout << "angle= " << angle << endl;
	
		angle = angle*vtkMath::Pi()/180;
		tilt_x=(virion_x-1024) - tan( angle ) * virion_z;
		tilt_x= round( tilt_x/sqrt(1+pow2( tan( angle ) ) ) );
	
		tilt_x+=1024;
		Range I( tilt_x-240, tilt_x+239);
		
		Array< float, 3 > H( series(I,J,Range::all() ) );
		cout << "Sub-series = " << H.shape() << ", [" << min(H) << "," << max(H) << "]" << endl;

		Array< float, 2 > slice( series( I, J, i ) );
		Array< float, 2 > vslice( virion( Range::all(),Range::all(), i ) );

		vslice = ( vslice - min(vslice) ) * ( max(slice) - min(slice) ) / ( max(vslice) - min(vslice) ) + min(slice) ;
		
		//cout << "CCC1(" << i << ") = " << mean( pow2( vslice - slice ) ) << endl;
		//slice = ( slice - min(slice) ) * ( max(vslice) - min(vslice) ) / ( max(slice) - min(slice) ) + min(vslice) ;				
		//cout << "CCC1(" << i << ") = " << mean( pow2( vslice - slice ) ) << endl;		

		dif( Range :: all(), Range :: all(), i ) = slice - vslice;
		Array< float, 2 > A( slice.shape() );
		A = slice;
		slice = slice - vslice;

		Range Iprev( tilt_x-241-479, tilt_x-241);
		Range Ipost( tilt_x+240, tilt_x+240+479);
		series( Iprev, J, i ) = A - vslice;
		series( Ipost, J, i ) = A - vslice;
		

		}
	
	nbfMrcWriter writer;
	writer.setFileName( argv[4] );
	writer.write(series);	

	writer.setFileName( argv[5] );
	writer.write(dif);	
#endif	
}
