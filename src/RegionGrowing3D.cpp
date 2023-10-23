#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <io/nbfImageReader.h>
#include <io/nbfImageWriter.h>
#include <fm/nbfFastMarchingFool3D.h>
#include <fm/nbfFastMarchingCurvature3D.h>

void main( int argv, char ** argc )
{
	if ( argv != 4 ){
		exit(0);
	}

	nbfImageReader reader;

	Timer t;

	float bandWidth = 10;
	float grayThres = .6588;

	// V - surface membrane
	Array< float, 3 > V;

	reader.setFileName( argc[1] );
	if ( reader.read( V ) ){
		cout << "Error reading file.\n";
		exit(0);
	}

	// input image
	reader.setFileName( argc[2] );
	Array< float, 3 > I;
	if ( reader.read( I ) ){
		cout << "Error reading file.\n";
		exit(0);
	}
	I = I - min(I);
	I = I / max(I);


	vector< TinyVector< int, 3 > > aliveP;
	vector< float > aliveD;

	// find seed points which are immediately outside the membrane
	Array< float, 3 > :: iterator iterV = V.begin();
	Array< float, 3 > :: iterator iterI = I.begin();
	while ( iterV != V.end() ){
		//if ( ( (*iterV) < 0 ) && ( (*iterV) > -.75 ) ){
		if ( ( (*iterV) < - 1 ) && ( (*iterV) > - 2 ) && ( (*iterI) < .5 ) ){
			if ( iterV.position()[secondDim] == 71 ){
				aliveP.push_back( iterV.position() );
				aliveD.push_back( -1 );
			}
			//cout << iterV.position() << endl;
		}
		++iterV; ++iterI;
	}


	Array< float, 3 > W( V.shape() );
	W = where( ( V < -1 ) && ( V > - bandWidth ), I, numeric_limits< float > :: max() );
	//W( Range::all(), Range(fromStart,52), Range::all() ) = numeric_limits< float > :: max();
	//W( Range::all(), Range(86,toEnd), Range::all() ) = numeric_limits< float > :: max();
	W[ aliveP ] = I;

	//nbfFastMarchingFool3D< float > grow( W, grayThres );
	//grow.setAliveSet( aliveP, aliveD );
	//
	Array< float, 3 > S( W.shape() );
	//grow.execute(S);

	nbfFastMarchingCurvature3D< float, float > fms( W );
	fms.setAliveSet( aliveP, aliveD );
	fms.setAlpha(0.05);
	fms.setNeighborhoodSize(1);
	fms.execute( S );

	// multi purpose writer
	nbfImageWriter writer;
	writer.setFileName( argc[3] );
	writer.writeFast(I);

}