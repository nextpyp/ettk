/***************************************************************************
  **************************************************************************
  
                SOFT: SO(3) Fourier transform code

                Version 1.0

  
   Peter Kostelec, Dan Rockmore
   {geelong,rockmore}@cs.dartmouth.edu
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
  
   Copyright 2003 Peter Kostelec, Dan Rockmore
  
  
     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
  
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
  
     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
  
  
   Commercial use is absolutely prohibited.
  
   See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/

/*

 a somewhat memory-friendly test routine to rotate a spherical function
 by massaging its S^2 Fourier coefficients with Wigner-D functions

 bw = bandwidth of input signal
 degOut = max degree of spherical harmonic you want to use ( < bw )
 alpha, beta, gamma -> the three Euler angles

             0 <= alpha, gamma < 2*pi
             0 <= beta <= pi

 inputSamples -> filename of input samples in INTERLEAVED format
 outputSamples -> filename of output (rotated) samples in INTERLEAVED format

 Here are order of rotation events:
  1) rotate by gamma about the z-axis
  2) rotate by beta about the y-axis
  3) rotate by alpha about the z-axis.

 example: test_s2_rotate_mem bw degOut alpha beta gamma inputSamples outputSamples

 example: test_s2_rotate_mem 32 31 0.37 2.32 4.37 fctIn.dat fctOut.dat


 NOTE: Sometimes there is a segmentation fault *after* all the rotating and
 writing out of the output file is complete. I haven't tracked this down yet,
 but I believe it has to do with freeing up the memory associated with doing
 the S^2 transforms ... my array of double pointers are not pointing in the
 right places when I try to free memory. However, the rotation itself is
 correct.

*/

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "FST_semi_memo.h"
#include "csecond.h"
#include "cospmls.h"
#include "primitive_FST.h"
#include "seminaive.h"
#include "rotate_so3_mem.h"


/* #define max(A, B) ((A) > (B) ? (A) : (B)) */

/**************************************************************/
/**************************************************************/


int main(int argc, char **argv)
{
  FILE *fp ;
  int i ;
  int bw, degOut ;
  REAL alpha, beta, gamma ;
  REAL *sigR, *sigI ;
  REAL *scratch ;
  double tstart, tstop ;
  double *seminaive_naive_tablespace ;
  double *trans_seminaive_naive_tablespace;
  double **seminaive_naive_table ;
  double **trans_seminaive_naive_table;

  if (argc < 3)
    {
      fprintf(stdout, "Usage: test_s2_rotate_mem bw degOut ");
      fprintf(stdout, "alpha beta gamma  ");
      fprintf(stdout, "input_filename output_filename\n");
      exit(0);
    }


  bw = atoi( argv[ 1 ] );
  degOut = atoi( argv[ 2 ] );
  alpha = (REAL) atof( argv[ 3 ] );
  beta = (REAL) atof( argv[ 4 ] );
  gamma = (REAL) atof( argv[ 5 ] );

  sigR = (REAL *) malloc(sizeof(REAL)*(4*bw*bw));
  sigI = (REAL *) malloc(sizeof(REAL)*(4*bw*bw));
  scratch = (REAL *) malloc(sizeof(REAL)*((10*bw*bw) + (48 * bw)));


  seminaive_naive_tablespace =
    (double *) malloc(sizeof(double) *
		      (Reduced_Naive_TableSize(bw,bw) +
		       Reduced_SpharmonicTableSize(bw,bw)));

  trans_seminaive_naive_tablespace =
    (double *) malloc(sizeof(double) *
		      (Reduced_Naive_TableSize(bw,bw) +
		       Reduced_SpharmonicTableSize(bw,bw)));

  /****
       At this point, check to see if all the memory has been
       allocated. If it has not, there's no point in going further.
  ****/

  if ( (scratch == NULL) || 
       (sigR == NULL ) || (sigI == NULL ) ||
       (seminaive_naive_tablespace == NULL) ||
       (trans_seminaive_naive_tablespace == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }
  

  fprintf(stdout,"Generating seminaive_naive tables...\n");
  seminaive_naive_table = SemiNaive_Naive_Pml_Table(bw, bw,
						    seminaive_naive_tablespace,
						    scratch);


  fprintf(stdout,"Generating trans_seminaive_naive tables...\n");
  trans_seminaive_naive_table =
    Transpose_SemiNaive_Naive_Pml_Table(seminaive_naive_table,
					bw, bw,
					trans_seminaive_naive_tablespace,
					scratch);


  fprintf(stdout,"reading in signal ...\n");

  /* read in signal */
  fp = fopen(argv[6], "r");
  for ( i = 0 ; i < (4*bw*bw) ; i ++ )
    {
      fscanf(fp,"%lf",sigR+i);
      fscanf(fp,"%lf",sigI+i);
    }
  fclose( fp ) ;

  fprintf(stdout,"about to rotate ...\n");
  tstart = csecond() ;

  rotateFct_mem( bw, degOut,
		 sigR, sigI,
		 alpha, beta, gamma,
		 scratch,
		 seminaive_naive_table,
		 trans_seminaive_naive_table ) ;

  tstop = csecond();
  fprintf(stdout,"finished rotating ...\n");
  fprintf(stdout,"rotation time \t = %.4e\n", tstop - tstart);
 
  /* write out rotated signal */
  fp = fopen(argv[7], "w");
  for ( i = 0 ; i < (4*bw*bw) ; i ++ )
    {
      fprintf(fp,"%.15f\n%.15f\n",sigR[i],sigI[i]);
    }
  fclose( fp ) ;
 
  fprintf(stdout,"finished writing ...\n");
 
  free(trans_seminaive_naive_table);
  free(seminaive_naive_table);
  free(seminaive_naive_tablespace);
 
  free(scratch);
  free(sigI);
  free(sigR);
 
  return 0 ;
 
}


