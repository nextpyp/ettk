#define BZ_GENERATE_GLOBAL_INSTANCES

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <io/nbfMrcReader.h>
#include <io/nbfMrcWriter.h>
#include <io/nbfImageWriter.h>
#include <fm/nbfFastMarchingFool3D.h>
#include <fm/nbfFastMarching.h>

#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkImageGaussianSmooth.h>

int main( int argc, char ** argv )
{
	if ( argc < 4 ){
		cout << "Usage: virus_dilate_membrane virion seg_th band_width" << endl;
		cout << "\tvirion (no extension)\n";
		cout << "\tseg_th (segmentation threshold)\n";
		cout << "\tvirion (distance to dilate membrane)\n";
		return 0;
	}

	nbfMrcReader reader;

	vtkImageData * data = vtkImageData :: New();
	Array< float, 3 > I, S;

	// read image to be masked
	stringstream input;
	input << argv[1] << ".mrc";

	reader.setFileName( input.str().c_str() );
	reader.read( data );

	nbfVTKInterface :: vtkToBlitz( data, I );

	// read implicit surface
	stringstream surface;
	surface << argv[1] << "_seg.mrc";
	reader.setFileName( surface.str().c_str() );
	reader.read( data );

	// if surface has different binning than volume
	if ( argc > 4 ){
		float bin_factor = atof( argv[4] );
		vtkImageResample * resample = vtkImageResample :: New();
		resample->SetInput( data );
		resample->SetAxisMagnificationFactor( 0, bin_factor );
		resample->SetAxisMagnificationFactor( 1, bin_factor );
		resample->SetAxisMagnificationFactor( 2, bin_factor );
		resample->SetInterpolationModeToCubic();
		resample->Update();
		data->DeepCopy( resample->GetOutput() );
		resample->Delete();
	}

	nbfVTKInterface :: vtkToBlitz( data, S );

	float threshold = atof( argv[2] );

	vector< TinyVector< int, 3 > > aliveP;
	vector< float > aliveD;

	// find seed points which are inside the threshold
	Array< float, 3 > :: iterator iterS = S.begin();
	while ( iterS != S.end() ){
		if ( (*iterS) > threshold ){
				aliveP.push_back( iterS.position() );
				aliveD.push_back( -1 );
		}
		++iterS;
	}

	// fast marching weights
	Array< float, 3 > W( S.shape() );
	W = 1;

	nbfFastMarching3D< float > grow( W );
	grow.setAliveSet( aliveP, aliveD );

	float dmax = atof( argv[3] );
	grow.setStopDistance( dmax );
	grow.execute(S);

	S = where( S < dmax, 1, 0 );

	vtkImageGaussianSmooth * filter = vtkImageGaussianSmooth::New();
	nbfVTKInterface :: blitzToVtk( S, data );
	filter->SetInput( data );
	filter->SetRadiusFactors( 2, 2, 2 );
	filter->Update();
	nbfVTKInterface :: vtkToBlitz( filter->GetOutput(), S );
	filter->Delete();

	I = 2 * ( I - min(I) ) / ( max(I) - min(I) );
	I = - ( I - mean(I) );
	I = I * S;
	// I = - where( S < dmax, I, max(I) );

	// multi purpose writer
	nbfMrcWriter writer;
	stringstream output;
	output << argv[1] << "_T" << argv[2] << ".mrc";

	writer.setFileName( output.str().c_str() );
	writer.write(I);
	data->Delete();

	return 1;
}