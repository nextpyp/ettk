#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfDifferentials.h>
#include <io/nbfMatlabWriter.h>
#include <bs/nbfBordStrategyMirror.h>

void main(){

	Array< float, 2 > input( 100, 100 );
	input = 5;
	Range I( input.rows() / 3, input.rows() * 2 / 3 );
	input( I, I ) = 2;

	Range J( input.rows() / 6, input.rows() * 2 / 6 );
	input( J, J ) = 2;

	Array< float, 2 > phi( input.shape() );
	firstIndex i;
	secondIndex j;
	phi = sqrt( pow2( i - input.rows() / 2.5 ) + pow2( j - input.cols() / 2.0 ) );
	phi = 10 - phi;

	nbfMatlabWriter writer;

	writer.setFileName("image");
	writer.write(input);

	writer.setFileName("phi");
	writer.write(phi);

	BordStrategyMirrorDouble< float, 2 > BSforPhi( phi, 2 );
	BSforPhi.refresh();

	Array< int, 2 > inside( input.shape() );
    inside = where( phi > 0, 1, 0 );

	float beta = 1.0;
	float nu = 0;
	float inSideMean = sum( input * inside ) / sum( inside );
	float outSideMean = sum( ( 1 - inside ) * input ) / sum( ( 1 - inside ) );

	int iters = 300;
	float lambda = 1e-2;

	Array< float, 2 > update( input.shape() );

	for ( int i = 0; i < iters; i++ ){
		update = beta * curvature2D(phi) + ( 
			-nu +
			-pow2( input - inSideMean ) + 
			pow2( input - outSideMean ) 
			) * sqrt( pow2( central12n(phi,firstDim) ) + pow2( central12n(phi,secondDim) ) );
		// update inside
        inside = where( phi > 0, 1, 0 );
		// update means
		inSideMean = sum( input * inside ) / sum( inside );
		outSideMean = sum( ( 1 - inside ) * input ) / sum( ( 1 - inside ) );
		// update phi
		phi = phi + lambda * update;
		BSforPhi.refresh();
	}

	cout << min(phi) << ", " << max(phi) << endl;

	writer.setFileName("phi_after");
	writer.write(phi);
}