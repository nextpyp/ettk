#define NBF_PARALLEL_IMPLEMENTATION_MPI 1
#define NBF_VERBOSE 1

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
#include <vtkImageContinuousDilate3D.h>
#include <vtkImageCast.h>
#include <vtkImageExtractComponents.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMrcWriter.h>
#include <bs/nbfBordStrategyMirror.h>

#include <em/nbfImageMetric.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfWedgedAverageImage3D.h>

#include <em/nbfLoopClustering.h>

#define PIXEL float

int main(int argc, char **argv)
{
	FILE *fp ;
	int i ;
	int bwIn, bwOut, degOut ;
	REAL alpha, beta, gamma ;
	REAL *sigInR, *sigInI, *sigOutR, *sigOutI ;
	REAL *scratch ;
	double tstart, tstop;
	double *seminaive_naive_tablespace, *trans_seminaive_naive_tablespace2;
	double *seminaive_naive_tablespace2 ;
	double **seminaive_naive_table2,**seminaive_naive_table ;
	double **trans_seminaive_naive_table2;

	if (argc < 9)
	{
		fprintf(stdout, "Usage: test_s2_rotate bwIn bwOut degOut ");
		fprintf(stdout, "alpha beta gamma ");
		fprintf(stdout, "input_filename output_filename\n");
		exit(0);
	}


	bwIn = atoi( argv[ 1 ] );
	bwOut = atoi( argv[ 2 ] );
	degOut = atoi( argv[ 3 ] );
	alpha = (REAL) atof( argv[ 4 ] );
	beta = (REAL) atof( argv[ 5 ] );
	gamma = (REAL) atof( argv[ 6 ] );

	for ( int rep = 0; rep < 1000; rep++ ){
		sigInR = (REAL *) malloc(sizeof(REAL)*(4*bwIn*bwIn));
		sigInI = (REAL *) malloc(sizeof(REAL)*(4*bwIn*bwIn));
		sigOutR = (REAL *) malloc(sizeof(REAL)*(4*bwOut*bwOut));
		sigOutI = (REAL *) malloc(sizeof(REAL)*(4*bwOut*bwOut));

		if ( bwOut > bwIn )
			scratch = (REAL *) malloc(sizeof(REAL)*((16*bwOut*bwOut) + (48 * bwOut)));
		else
			scratch = (REAL *) malloc(sizeof(REAL)*((16*bwIn*bwIn) + (48 * bwIn)));

		seminaive_naive_tablespace =
			(double *) malloc(sizeof(double) *
			(Reduced_Naive_TableSize(bwIn,bwIn) +
			Reduced_SpharmonicTableSize(bwIn,bwIn)));

		trans_seminaive_naive_tablespace2 =
			(double *) malloc(sizeof(double) *
			(Reduced_Naive_TableSize(bwOut,bwOut) +
			Reduced_SpharmonicTableSize(bwOut,bwOut)));

		seminaive_naive_tablespace2 =
			(double *) malloc(sizeof(double) *
			(Reduced_Naive_TableSize(bwOut,bwOut) +
			Reduced_SpharmonicTableSize(bwOut,bwOut)));



		/****
		At this point, check to see if all the memory has been
		allocated. If it has not, there's no point in going further.
		****/

		if ( (scratch == NULL) || 
			(sigInR == NULL ) || (sigInI == NULL ) ||
			(sigOutR == NULL ) || (sigOutI == NULL ) ||
			(seminaive_naive_tablespace == NULL) ||
			(trans_seminaive_naive_tablespace2 == NULL) )
		{
			perror("Error in allocating memory");
			exit( 1 ) ;
		}


		fprintf(stdout,"Generating seminaive_naive tables...\n");
		seminaive_naive_table = SemiNaive_Naive_Pml_Table(bwIn, bwIn,
			seminaive_naive_tablespace,
			scratch);


		fprintf(stdout,"Generating seminaive_naive tables...\n");
		seminaive_naive_table2 = SemiNaive_Naive_Pml_Table(bwOut, bwOut,
			seminaive_naive_tablespace2,
			scratch);


		fprintf(stdout,"Generating trans_seminaive_naive tables...\n");
		trans_seminaive_naive_table2 =
			Transpose_SemiNaive_Naive_Pml_Table(seminaive_naive_table2,
			bwOut, bwOut,
			trans_seminaive_naive_tablespace2,
			scratch);

		fprintf(stdout,"reading in signal ...\n");

		/* read in signal */
		//fp = fopen(argv[7], "r");
		//for ( i = 0 ; i < (4*bwIn*bwIn) ; i ++ )
		//  {
		//    fscanf(fp,"%lf",sigInR+i); cout << "R=\t" << *(sigInR+i) << endl;
		//    fscanf(fp,"%lf",sigInI+i); cout << "I=\t" << *(sigInI+i) << endl;
		//  }
		//fclose( fp ) ;

		cout << sizeof(REAL) << endl;
		Array< float, 2 > A;
		nbfMatlabReader reader;
		reader.setFileName("a.matlab");
		reader.read(A);

		Array< double, 2 > B( A.shape() );
		B = cast< double >( A );

		for ( int i = 0; i < B.size(); i++ ){
			*(sigInR+i) = *( B.data() + i );
			*(sigInI+i) = 0;
		}

		fprintf(stdout,"about to rotate ...\n");
		tstart = csecond();

		rotateFct( bwIn, bwOut, degOut,
			sigInR, sigInI,
			sigOutR, sigOutI,
			alpha, beta, gamma,
			scratch,
			seminaive_naive_table,
			trans_seminaive_naive_table2 ) ;

		tstop = csecond();
		fprintf(stdout,"finished rotating ...\n");
		fprintf(stdout,"rotation time \t = %.4e\n", tstop - tstart);

		/* write out rotated signal */
		//fp = fopen(argv[8], "w");
		//for ( i = 0 ; i < (4*bwOut*bwOut) ; i ++ )
		//  {
		//    fprintf(fp,"%.15f\n%.15f\n",sigOutR[i],sigOutI[i]);
		//  }
		//fclose( fp ) ;

		fprintf(stdout,"finished writing ...\n");

		for ( int i = 0; i < B.size(); i++ ){
			*( B.data() + i ) = sigOutR[i];
		}

		nbfMatlabWriter writer;
		writer.setFileName(argv[8]);
		writer.write(B);

		free(trans_seminaive_naive_table2);
		free(seminaive_naive_table2);
		free(seminaive_naive_table);

		free(seminaive_naive_tablespace2);
		free(trans_seminaive_naive_tablespace2);
		free(seminaive_naive_tablespace);

		free(scratch);
		free(sigOutI);
		free(sigOutR);
		free(sigInI);
		free(sigInR);
	}
	return 0 ;
}