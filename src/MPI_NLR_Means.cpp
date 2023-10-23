#define NBF_PARALLEL_IMPLEMENTATION_MPI 1

#include "mpi.h"

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <string.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <nbfPolarDomain.h>
#include <nbfMaximalFlow.h>
#include <bs/nbfBordStrategyMirror.h>

#include <fftw3.h>

#define PIXEL double

int main( int argc, char ** argv )
{
	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	//  Initialize MPI.
	MPI_Init ( &argc, &argv );

	int my_id, num_procs;

	//  Get this processes's rank.
	MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

	//  Find out how many processes are available.
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	cout << "Process " << my_id << " is active.\n";

	// define control tags for all processes
	int tag_invariant_representation = 0;
	int tag_nlmeans = 1;
	int tag_nlmeans_initialization = 2;
	int tag_done = 3;

	// define size of rotation invariant representation
	int fftSize = 32;
	
	// define size of control data packets
	int controlSize = 1;

	// define number of pixels handled by each process at a time
	int pixelsPerProcess = 100;

	// START

	if ( argc != 5 ){
		cout << "Usage: input_image patch_radius exponential_decay output_image" << endl;
		return 0;
	}

	// input is shared by all processes
	Array< PIXEL, 2 > input;

	char * inputFile = argv[1];
	nbfMatlabReader mr;
	mr.setFileName( inputFile );
	Array< float, 2 > tmp;
	mr.read( tmp );
	input.resize( tmp.shape() );
	input = cast<double>(tmp);

	cout << "Input image geometry = " << input.shape() << endl;

	// radius of local patches
	int rsize = atoi( argv[2] );

	// exponential decay for computing pixel weights
	PIXEL decay = atof( argv[3] );

	// filename to save rotation invariant image representation
	stringstream invariantFileName;
	invariantFileName << inputFile << ".RI.matlab";

	if ( my_id == 0 ){ // Master process

		cout << "MPI - Master process:\n";
		cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
		cout << "  The number of processes is " << num_procs << "\n";

		double * controlData = new double[ controlSize ];

		// initialize output image
		Array< PIXEL, 2 > output( input.shape() );

		// compute rotation invariant image representation
		Array< PIXEL, 3 > invariant;
		invariant.resize( input.rows(), input.cols(), 2 + fftSize );
		cout << "invariant = " << invariant.shape() << endl;

		int index = 0;

		cout << "computing RI representation...";

		// start all available processes
		for ( int process = 1; process < num_procs; process++ ){
			if ( index < input.size() ){
				// send control information: index
				controlData[0] = index;
				MPI_Send ( controlData, controlSize, MPI_DOUBLE, process, tag_invariant_representation, MPI_COMM_WORLD );
				index += pixelsPerProcess; 
			}
		}

		int jobsDone = 0;
		
		double * results = new double[ 1 + ( 2 + fftSize ) * pixelsPerProcess ];

		// loop while not all jobs submitted
		while ( jobsDone < input.size() ){

			MPI_Status status;

			// receive result from slaves
			MPI_Recv( results, 1 + ( 2 + fftSize ) * pixelsPerProcess, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			int source = status.MPI_SOURCE;
			int tag = status.MPI_TAG;

			if ( tag != tag_invariant_representation ){
				cerr << "ERROR in message sequence. Expecting the result of RI operation. Received = " << tag << endl;
				return 1;
			}

			jobsDone += pixelsPerProcess;

			// retrieve current index
			int currentIndex = results[0];

			// assign RI representations
			int count = 1;
			for ( int pixel = currentIndex; pixel < currentIndex + pixelsPerProcess; pixel++ ){
				// if still inside image
				if ( pixel < input.size() ){
					int row = floor( 1.0 * pixel / input.rows() );
					int col = fmod( 1.0 * pixel, input.rows() );
					for ( int i = 0; i < 2 + fftSize; i++ ){
						invariant( row, col, i ) = results[count];
						count++;
					}
				}
			}

			// send remaining jobs if not already submitted
			if ( index < input.size() ){
				controlData[0] = index; // pixel index
				MPI_Send ( controlData, controlSize, MPI_DOUBLE, source, tag_invariant_representation, MPI_COMM_WORLD ); 
				index += pixelsPerProcess;
			}
		}

		delete [] results;

		w.setFileName( invariantFileName.str().c_str() );
		w.write(invariant);

		// DONE computing invariant representation, NOW paralelize denoising algorithm

		cout << "done.\ndenoising...";

		// signal all processes that invariant representation is ready to be retrieved from file
		for ( int process = 1; process < num_procs; process++ ){	
			// send control info
			MPI_Send ( controlData, controlSize, MPI_DOUBLE, process, tag_nlmeans_initialization, MPI_COMM_WORLD );
		}

		index = 0; 
	
		// start all available processes
		for ( int process = 1; process < num_procs; process++ ){
			if ( index < input.size() ){
				controlData[0] = index;
				MPI_Send ( controlData, controlSize, MPI_DOUBLE, process, tag_nlmeans, MPI_COMM_WORLD );
				index += pixelsPerProcess;
			}
		}

		int indexesDone = 0;

		results = new double[ 1 + pixelsPerProcess ];

		// loop while not all jobs submitted
		while ( indexesDone < input.size() ){

			MPI_Status status;

			// receive results from slaves
			MPI_Recv( results, 1 + pixelsPerProcess, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			int source = status.MPI_SOURCE;
			int tag = status.MPI_TAG;

			if ( tag != tag_nlmeans ){
				cerr << "ERROR in message sequence. Expecting the result of NLRmeans operation (here). Received = " << tag << endl;
				return 1;
			}
			
			indexesDone += pixelsPerProcess;

			for ( int i = 0; i < pixelsPerProcess; i++ ){
				int next = results[0] + i;
				if ( next < output.size() ){
					output.data()[ next ] = results[ 1 + i ];
				}
			}

			// send remaining jobs if not already submitted
			if ( index < input.size() ){
				controlData[0] = index;
				MPI_Send ( controlData, controlSize, MPI_DOUBLE, source, tag_nlmeans, MPI_COMM_WORLD );
				index += pixelsPerProcess;
			}
		}

		delete [] results;
		
		cout << "done. terminating..." << endl;

		// end MPI

		// send termination message to all nodes
		for ( int process = 1; process < num_procs; process++ ){
			cout << "finishing " << process << "...";
			MPI_Send ( controlData, controlSize, MPI_DOUBLE, process, tag_done, MPI_COMM_WORLD );
			cout << "done." << endl;
		}

		// save output to file
		cout << "done. Writing output...";
		w.setFileName( argv[4] );
		w.write(output);
		cout << "done." << endl;
		
		delete [] controlData;

		// end Master process

	} else {

		// Slave process

		int source = 0;
		int tag = -1;

		MPI_Status status;

		double * controlData = new double[ controlSize ];

		fftw_plan  fftplan;
		fftw_plan ifftplan;

		Array< complex< PIXEL >, 1 > data1, data2;
		fftw_complex * inFFT;
		fftw_complex * outFFT;

		Array< PIXEL, 3 > invariant;

		while ( tag != tag_done ){

			// retrieve pixel index and tag information
			MPI_Recv( controlData, controlSize, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			tag = status.MPI_TAG;

			if ( tag == tag_invariant_representation ){

				// index + fftSize invariant representations
				double * result = new double[ 1 + ( 2 + fftSize ) * pixelsPerProcess ];
				result[0] = controlData[0];

				// store local image patch
				Array< PIXEL, 2 > local( rsize * 2 + 1, rsize * 2 + 1 );
				TinyVector< PIXEL, 2 > center( ( local.shape() - 1 ) / 2 );
						
				int count = 1;
				
				// for all pixels this process is assigned
				for ( int pixelOffset = 0; pixelOffset < pixelsPerProcess; pixelOffset++ ){

					int current = controlData[0] + pixelOffset;
					int row = floor( 1.0 * current / input.rows() );
					int col = current % input.rows();

					if ( current < input.size() ){

						// get local patch
						for ( int p = -rsize; p <= rsize; p++ ){
							for ( int q = -rsize; q <= rsize; q++ ){
								if ( input.isInRange( row + p, col + q ) ){
									local( p + rsize, q + rsize ) = input( row + p, col + q );
								} else {
									local( p + rsize, q + rsize ) = - numeric_limits< PIXEL > :: max();
								}
							}
						}

						// replace outside pixels with mean
						PIXEL background = sum( where( local == - numeric_limits< PIXEL > :: max(), 0, local ) ) / sum( where( local == - numeric_limits< PIXEL > :: max(), 0, 1 ) );
						local = where( local == - numeric_limits< PIXEL > :: max(), background, local );

						// set center to 0
						// local( center ) = 0;

						nbfMaximalFlow< PIXEL > flow;

						Array< PIXEL, 2 > press( local.shape() );

						// set source to image center
						press( center ) = 1;
						// set sink to image boundary
						press( press.lbound(0), Range::all() ) = -1;
						press( press.ubound(0), Range::all() ) = -1;
						press( Range::all(), press.lbound(1) ) = -1;
						press( Range::all(), press.ubound(1) ) = -1;

						// convert to polar and integrate
						Array< PIXEL, 2 > P;
						Array< bool, 2 > B;
						nbfPolarDomain< PIXEL, 2 > polar;

						PIXEL rho = rsize;
						polar.setCenter( center );
						polar.setMaxRho( rho );
						polar.setResRho( 32 );
						polar.setResTheta( 32 );
						polar.cartesian2polar( local, P, B );

						// edge detector
						BordStrategyMirrorDouble< PIXEL, 2 > bsForLocal( local, 1 );
						bsForLocal.refresh();
						// local = max(local) - local;
						Array< PIXEL, 2 > edges( local.shape() );
						edges = sqrt( pow2( forward11n(local,firstDim) ) + pow2( forward11n(local,secondDim) ) );
						edges = max(edges) - edges;

						// compute scaling profile (max-flow/min-cut)
						flow.execute( press, edges, 100 );

						// make interior to 1 and exterior to 0
						// press = where( press >= 0, 1, 0 );

						// convert to polar
						Array< PIXEL, 2 > Pscale;
						polar.cartesian2polar( press, Pscale, B );
						Pscale = where( Pscale >= 0, 1, 0 );

						Array< PIXEL, 1 > N( P.rows() );
						firstIndex i; secondIndex j;
						N = sum( Pscale(j,i), j );
						int nmax = max(N);
						int n = Pscale.cols() - nmax;

						Array< PIXEL, 2 > Mi( Pscale.shape() );
						Mi( Range(n,toEnd), Range::all() ) = Pscale( Range(fromStart,nmax-1), Range::all() );
						Mi( Range(0,n-1), Range::all() ) = 1;

						// compute integral up to scale factor
						//P = P * Mi;
						Array< PIXEL, 1 > L( P.cols() );
						L = sum( P(j,i), j );

						//N = sum( Mi(j,i), j );

						// normalize line integrals by scale factor
						//L = L / N;

						// store result
						result[count] = mean( local ); 
						result[count+1] = sqrt( sum( ( local - result[count] ) * ( local - result[count] ) ) );
						count+=2;
						for ( int k = 0; k < fftSize; k++ ){
							result[ count ] = L(k);
							count++;
						}
					}
				} // end for

				// send result and notify master node we have finished processing
				MPI_Send ( result, 1 + ( 2 + fftSize ) * pixelsPerProcess, MPI_DOUBLE, source, tag_invariant_representation, MPI_COMM_WORLD );

				delete [] result;

			} else if ( tag == tag_nlmeans_initialization ){
				
				// read rotation invariant representation from file
				nbfMatlabReader reader;
				reader.setFileName( invariantFileName.str().c_str() );
				reader.read( invariant );

				// initialize fftw pipeline
				data1.resize( fftSize );
				data2.resize( fftSize );
				inFFT  = reinterpret_cast<fftw_complex*>( data1.data() );
				outFFT  = reinterpret_cast<fftw_complex*>( data2.data() );

				fftplan = fftw_plan_dft_1d( fftSize, inFFT, outFFT, FFTW_FORWARD, FFTW_MEASURE );
				ifftplan = fftw_plan_dft_1d( fftSize, outFFT, inFFT, FFTW_BACKWARD, FFTW_MEASURE );

			} else if ( tag == tag_nlmeans ){
				
				int currentIndex = controlData[0];

				double * result = new double[ 1 + pixelsPerProcess ];
				result[0] = currentIndex;
				int count = 1;
				
				for ( int pixelOffset = 0; pixelOffset < pixelsPerProcess; pixelOffset++ ){

					int current = currentIndex + pixelOffset;
					int row = floor( 1.0 * current / input.rows() );
					int col = current % input.rows();

					if ( current < input.size() ){
						PIXEL average = 0;
						PIXEL wtotal = 0;

						//PIXEL mean_1 = mean( invariant( row, col, Range(2,toEnd) ) );
						//PIXEL sigma_1 = sum( pow2( invariant( row, col, Range(2,toEnd) ) - mean_1 ) );

						real( data1 ) = invariant( row, col, Range(2,toEnd) );
						//real( data1 ) = ( invariant( row, col, Range(2,toEnd) ) - mean_1 ) / sqrt( sigma_1 );
						//real( data1 ) = invariant( row, col, Range(2,toEnd) ) - mean_1;
						imag( data1 ) = 0;
						fftw_execute(fftplan);

						Array< complex< PIXEL >, 1 > fft1( data1.shape() );
						fft1 = data2;

						Array< PIXEL, 2 > weights( input.shape() );

						for ( int i = 0; i < invariant.rows(); i++ ){
							for ( int j = 0; j < invariant.cols(); j++ ){

								//PIXEL mean_2 = mean( invariant( i, j, Range(2,toEnd) ) );
								//PIXEL sigma_2 = sum( pow2( invariant( i, j, Range(2,toEnd) ) - mean_2 ) );
								
								real( data1 ) = invariant( i, j, Range(2,toEnd) );
								//real( data1 ) = ( invariant( i, j, Range(2,toEnd) ) - mean_2 ) / sqrt( sigma_2 );
								//real( data1 ) = invariant( i, j, Range(2,toEnd) ) - mean_2;
								fftw_execute(fftplan);
								
								data2 = fft1 * conj( data2 );
								fftw_execute(ifftplan);

								// Normalized cross-correlation
								//PIXEL d = max(  real( data1 ) / ( 1.0 * data1.size() ) );
								//PIXEL d = .5 * ( 1 - max(  real( data1 ) / ( 1.0 * data1.size() ) ) );
								

								// L^2 norm
								PIXEL d = min( sum( pow2( invariant( row, col, Range(2,toEnd) ) ) ) + sum( pow2( invariant( i, j, Range(2,toEnd) ) ) ) - 2.0 * real( data1 ) / ( 1.0 * data1.size() ) );
								//PIXEL d = min( sum( pow2( invariant( row, col, Range(2,toEnd) ) - mean_1 ) ) + sum( pow2( invariant( i, j, Range(2,toEnd) ) - mean_2 ) ) - 2.0 * real( data1 ) / ( 1.0 * data1.size() ) );

								PIXEL w = exp( - d / decay / decay );

								weights(i,j) = w;
								// average += ( w * input(i,j) );
								
								////average += ( w * ( ( input(i,j) - invariant(i,j,0) ) / invariant(i,j,1) * invariant(row,col,1) + invariant(row,col,0) ) );
								//if ( fabs( invariant(i,j,0) - invariant(row,col,0) ) < 1.5e-1 ){
								//	average += w * input(i,j);
								//} else {
									//average += ( w * ( input(i,j) - invariant(i,j,0) + invariant(row,col,0) ) );
								//}
								//average += ( w * ( input(i,j) / mean_2 * mean_1 ) );
								average += w * input(i,j);

								wtotal += w;
							}
						}

						result[ count ] = average / wtotal;
						result[ count ] = input(row,col);
						count++;
					}
				}

				MPI_Send ( result, 1 + pixelsPerProcess, MPI_DOUBLE, source, tag_nlmeans, MPI_COMM_WORLD );
				delete [] result;
			}
		} // end while

		delete [] controlData;

		fftw_destroy_plan( fftplan );
		fftw_destroy_plan( ifftplan );

	}

	MPI_Finalize();

	return 0;
}