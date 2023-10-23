#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <nbf3DReader.h>
#include <nbfMatlabWriter.h>

void main( int argv, char ** argc )
{
	if ( argv != 10 ){
		cout << "Parameters:" << endl;
		cout << "\t 1 : input raw file (.raw)" << endl;
		cout << "\t 2 : rows      (=1024)" << endl;
		cout << "\t 3 : cols      (= 130)" << endl;
		cout << "\t 4 : depth     (=1024)" << endl;
		cout << "\t 5 : center[x] (= 550)" << endl;
		cout << "\t 6 : center[y] (= 472)" << endl;
		cout << "\t 7 : center[z] (=  65)" << endl;
		cout << "\t 8 : rho       (=  92)" << endl;
		cout << "\t 9 : output file (.array)" << endl;
		exit(0);
	}

	nbfMatlabWriter writer;

	Timer t;
	double total;
	total = 0;

	// 1. V - input image

	Array< float, 3 > V;

	nbf3DReader rawReader;

	// set mode to big endian if needed
	//rawReader.setBigEndian(true);

	// specify volume size (normally from file name)
	rawReader.setDimensions( atoi(argc[2]), atoi(argc[3]), atoi(argc[4]));

	// set type to short (signed) 16 bits
	rawReader.setDataType(1);
	rawReader.setFileName( argc[1] );

	cout << "\tVraw<short> = ";

	Array< short, 3 > Vs;
	rawReader.read( Vs );
	Vs.transposeSelf(firstDim,thirdDim,secondDim);

	cout << "\tVraw<short> = [" << min(Vs) << "," << max(Vs) << "]" << endl;
	
	// B2
	int centerX = atoi(argc[5]); // 2*275;
	int centerY = atoi(argc[6]); // 2*236;
	int centerZ = atoi(argc[7]); // 65;

	int rho = atof( argc[8] ); // 92
	int lboundX = blitz::minmax::max( centerX - rho, Vs.lbound(firstDim) );
	int lboundY = blitz::minmax::max( centerY - rho, Vs.lbound(secondDim) );
	int uboundX = blitz::minmax::min( centerX + rho, Vs.ubound(firstDim) );
	int uboundY = blitz::minmax::min( centerY + rho, Vs.ubound(secondDim) );

	if ( lboundX == Vs.lbound(firstDim) ){
		uboundX = 2 * centerX;
	}
	if ( uboundX == Vs.ubound(firstDim) ){
		lboundX = centerX - ( Vs.ubound(firstDim) - centerX );
	}

	if ( lboundY == Vs.lbound(secondDim) ){
		uboundY = 2 * centerY;
	}
	if ( uboundY == Vs.ubound(secondDim) ){
		lboundY = centerY - ( Vs.ubound(secondDim) - centerY );
	}

//	cout << "center = [" << centerX << "," << centerY << "," << centerZ << "]" << endl;
//	cout << "V<float> = V<short>[" << lboundX << ":" << uboundX << "," << lboundY << ":" << uboundY << ",:]" << endl;

	// crop 3D data
	V.resize( uboundX - lboundX + 1, uboundY - lboundY + 1, Vs.depth() );
	V = cast<float>( Vs( Range(lboundX,uboundX), Range(lboundY,uboundY), Range::all() ) );
	Vs.free();

	cout << "\tVcrop<float> = " << V.shape() << endl;
	cout << "\tVcrop<float> = [" << min(V) << "," << max(V) << "]" << endl;

	writer.setFileName( argc[9] );
	writer.write( V );
}