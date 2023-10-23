#define BZ_GENERATE_GLOBAL_INSTANCES

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#ifndef WIN32
	#include "mpi.h"
#endif

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
#include <vtkXMLPPolyDataWriter.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMODReader.h>
#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMrcWriter.h>

#include <em/nbfImageMetric.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfWedgedAverageImage3D.h>

#include <em/nbfExtractPointsAndNormals3D.h>

#define PIXEL double

#define USING_PHANTOM_LIST 0
#define LOOP_CREATE_VOLUME_LIST_UNDEF -1
#define LOOP_CREATE_VOLUME_LIST_POS 0
#define LOOP_CREATE_VOLUME_LIST_EUL 1
#define LOOP_CREATE_VOLUME_LIST_MRC 2
#define LOOP_CREATE_VOLUME_LIST_MOD 3

int main( int argc, char ** argv )
{
	cout << argv[0] << "\n";
	cout << " Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";

	if ( argc < 2 ){
		cout << "USAGE:\neul2cmm eul_file (no extension)\n" << endl;
		return 0;
	}

	// PARSE PARAMETERS

	stringstream positionsFile;
	positionsFile << argv[ 1 ] << ".eul";

	nbfMODReader mr;
	mr.setFileName( positionsFile.str().c_str() );
	Array< float, 2 > points;
	mr.readRawEul( points );
	if ( points.cols() != 6 ){
		cout << "ERROR: eul file format not recognized." << endl;
		return 0;
	}

	if ( points.size() > 0 ){
		cout << points.rows() << " locations converted.\n";
	}
	else{
		cout << "ERROR reading " << positionsFile.str() << ". No locations will be extracted." << endl;
		return 0;
	}

	stringstream chimeraFile;
	chimeraFile << argv[1] << ".cmm";
	ofstream modelFile( chimeraFile.str().c_str(), ios::out );
	stringstream outputStream;
	outputStream << "<marker_sets>\n<marker_set name=\"spikes detected automatically\">\n";

	vtkPoints * vpoints = vtkPoints :: New();

	for ( int i = 0; i < points.rows(); i++ ){
		outputStream << "<marker id=\"" << i + 1 << "\" x=\"" << points(i,0) << "\" y=\"" << points(i,1) << "\" z=\"" << points(i,2) << "\" r=\"1\" g=\"0.2\" b=\"1\" radius=\"2\"/>\n";
		vpoints->InsertNextPoint( points(i,0),  points(i,1),  points(i,2) );
	}
	outputStream << "</marker_set>\n</marker_sets>";
	modelFile << outputStream.str();
	modelFile.close();

	vtkPolyData * data = vtkPolyData :: New();
	data->SetPoints( vpoints );
	vpoints->Delete();

	vtkXMLPPolyDataWriter * writer = vtkXMLPPolyDataWriter::New();
	writer->SetInput( data );
	stringstream vtkFile;
	vtkFile << argv[1] << ".xml";
	writer->SetFileName( vtkFile.str().c_str() );
	writer->Write();
	data->Delete();
	
	return 0;
}