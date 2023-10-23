#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

#include <nbfTimer.h>

using namespace blitz;

#define PIXEL float

#include <vtkImageData.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>

#include <vector>

#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfImageFilter.h>
#include <em/nbfFourierFilter.h>
#include <em/nbfFourierImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>

#include <vnl/algo/vnl_amoeba.h>

int main( int argc, char *argv[] )
{
	vector< nbfWedgedSubImage3D< PIXEL > > volumeList;
	
	stringstream volumeFile;
	volumeFile << argv[1] << ".txt";
	nbfWedgedSubImage3D< PIXEL > :: read( volumeFile.str().c_str(), volumeList );

	nbfImageFilter< PIXEL, 3 > imageFilter;
	nbfFourierFilter< PIXEL, 3 > fourierFilter;
	fourierFilter.bandPassOn(.05,.2,.01);

	nbfProjectionRotationMetric3D< PIXEL > metric( &imageFilter, &fourierFilter );
	metric.setNumberOfCandidatePeaksToSearch(15);

	//nbfFourierImageMetric< PIXEL, 3 > metric( &imageFilter, &fourierFilter );
	metric.setInput1( &volumeList[0] );
	for ( int i = 1; i < volumeList.size(); i++ ){
	//for ( int i = 1; i < 2; i++ ){
		volumeList[i].setTransform( (vtkTransform*)NULL);
		//volumeList[i].getWedge()->set( -90, 90 );
		metric.setInput2( &volumeList[i] );
		nbfTimer timer;
		timer.start();
		metric.getDistance();
		timer.stop();
		cout << "D = " << metric.getDistance() << ", (elapsed time " << timer.elapsedSeconds() << ")" << endl;
		volumeList[i].setTransform( metric.getTransform() );
		//return 0;
	}

	// save aligned volume series
	stringstream alignedFile;
	alignedFile << argv[1] << ".aligned.txt";
	nbfWedgedSubImage3D< PIXEL > :: write( alignedFile.str().c_str(), volumeList );

	// retrieve ground_truth
	vector< nbfWedgedSubImage3D< PIXEL > > groundVolumeList;
	nbfWedgedSubImage3D< PIXEL > :: read( volumeFile.str().c_str(), groundVolumeList );

	// compute error measure
	Array< PIXEL, 1 > error( volumeList.size() - 1 );
	error = 0;
	for ( int vols = 0; vols < volumeList.size(); vols++ ){
		vtkTransform * ground_truth = vtkTransform::New();
		groundVolumeList[vols].getTransform( ground_truth );
		vtkTransform * result = vtkTransform::New();
		volumeList[vols].getTransform( result );
		for ( int i = 0; i < 3; i++ ){
			for ( int j = 0; j < 3; j++ ){
				error(vols-1) += pow2( ground_truth->GetMatrix()->GetElement(i,j) - result->GetMatrix()->GetElement(i,j) );
			}
		}
		ground_truth->Delete();
		result->Delete();
	}
	error = sqrt( error );

	nbfMatlabWriter mw;
	stringstream errorFile;
	errorFile << argv[1] << ".error.matlab";
	mw.setFileName( errorFile.str().c_str() );
	mw.write( error );

	cout << "Error = " << error << endl;

	return 0;

	vnl_vector< double > seed(3);
	//seed[0] = 4.60347;
	//seed[1] = -1.90428;
	//seed[2] = 151.81;
	//seed[0] = 80;
	//seed[1] = 0;
	//seed[2] = 120;
	//seed[0] = 3.51027;
	//seed[1] = -0.418487;
	//seed[2] = 149.77;
	seed[0] = seed[1] = seed[2] = 0;

	cout << "seed = " << seed[0] << ", " << seed[1] << ", " << seed[2] << endl;
#if 0
	PIXEL size = 4.0;
	PIXEL delta = .2;

	Array< PIXEL, 3 > C( 2*size+1, 2*size+1, 2*size+1 );

	int indexI, indexJ, indexK;
	indexI = indexJ = indexK = 0;
	for ( PIXEL i = - size; i <= size; i++ ){
		for ( PIXEL j = - size; j <= size; j++ ){
			for ( PIXEL k = - size; k <= size; k++ ){
				vnl_vector< PIXEL > current(3);
				current[0] = seed[0] + delta * i;
				current[1] = seed[1] + delta * j;
				current[2] = seed[2] + delta * k;
				C(indexI,indexJ,indexK) = metric.f(current);
				cout << "[" << i << "," << j << "," << k << "] - " << C(indexI,indexJ,indexK) << endl; 
				indexK++;
			}
			indexK = 0;
			indexJ++;
		}
		indexJ = 0;
		indexI++;
	}

	nbfMatlabWriter w;
	w.setFileName("p.matlab");
	w.write(C);
//#else
	//vnl_amoeba optimizer( metric );
	vnl_powell optimizer( &metric );


	optimizer.set_max_function_evals(5);
	optimizer.set_verbose(true);
	optimizer.set_trace(true);
	optimizer.set_x_tolerance(1e-4);
	optimizer.minimize(seed);

	cout << "minimzer = " << seed[0] << ", " << seed[1] << ", " << seed[2] << ", F = " << optimizer.get_end_error() << endl;
#endif
  return 0;
}

