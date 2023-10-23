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
  to test the correlation routines

  - uses fftw for the fft in the inverse SO(3) transform, but
  - it's still the same ol' plain fft in the S^2 routines

  - ASSUMES bwIn >= bwOut



  - uses the Wigner-d symmetries and FFTW
  - uses part of SpharmonicKit
  - INTERLEAVED (i.e. real/imaginary) SAMPLES of signal and pattern files
  - [result] -> optional -> filename of all the correlation values
                (if you want all of them)
  - bwIn -> bw of input spherical signals
  - bwOut -> bw of so(3) transform you want to do
  - degLim -> max degree of Wigner-D functions you'll be using


  ASSUMES bwIn >= bwOut

  example: test_soft_fftw_correlate2 signalFile patternFile bwIn bwOut degLim [result]



*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fftw3.h"
#include "complex.h" 
#include "csecond.h"
#include "so3_correlate_fftw.h"
#include "soft_fftw.h"

#include "FST_semi_memo.h"
#include "cospmls.h"
#include "primitive_FST.h"
#include "seminaive.h"
#include "utils_so3.h"


int main ( int argc,
	   char **argv )
{
  FILE *fp ;
  int i ;
  int n, bwIn, bwOut, degLim ;
  double tstart, tstop ;
  fftw_complex *workspace1, *workspace2  ;
  REAL *workspace3 ;
  REAL *sigR, *sigI ;
  REAL *sigCoefR, *sigCoefI ;
  REAL *patCoefR, *patCoefI ;
  fftw_complex *so3Sig, *so3Coef ;
  fftw_plan p1 ;
  int na[2], inembed[2], onembed[2] ;
  int rank, howmany, istride, idist, ostride, odist ;
  int tmp, maxloc, ii, jj, kk ;
  REAL maxval ;
  double *seminaive_naive_tablespace  ;
  double **seminaive_naive_table ;


  if (argc < 6 )
    {
      printf("test_soft_fftw_correlate2 signalFile patternFile ");
      printf("bwIn bwOut degLim [result]\n");
      exit(0) ;
    }

  bwIn = atoi( argv[3] );
  bwOut = atoi( argv[4] );
  degLim = atoi( argv[5] );

  n = 2 * bwIn ;

  sigR = (REAL *) calloc( n * n, sizeof(REAL) );
  sigI = (REAL *) calloc( n * n, sizeof(REAL) );
  so3Sig = fftw_malloc( sizeof(fftw_complex) * (8*bwOut*bwOut*bwOut) );
  workspace1 = fftw_malloc( sizeof(fftw_complex) * (8*bwOut*bwOut*bwOut) );
  workspace2 = fftw_malloc( sizeof(fftw_complex) * ((14*bwIn*bwIn) + (48 * bwIn)));
  workspace3 = (REAL *) malloc( sizeof(REAL) * (12*n + n*bwIn));
  sigCoefR = (REAL *) malloc( sizeof(REAL) * bwIn * bwIn ) ;
  sigCoefI = (REAL *) malloc( sizeof(REAL) * bwIn * bwIn ) ;
  patCoefR = (REAL *) malloc( sizeof(REAL) * bwIn * bwIn ) ;
  patCoefI = (REAL *) malloc( sizeof(REAL) * bwIn * bwIn ) ;
  so3Coef = fftw_malloc( sizeof(fftw_complex) * ((4*bwOut*bwOut*bwOut-bwOut)/3) ) ;


  seminaive_naive_tablespace =
    (double *) malloc(sizeof(double) *
		      (Reduced_Naive_TableSize(bwIn,bwIn) +
		       Reduced_SpharmonicTableSize(bwIn,bwIn)));


  /****
       At this point, check to see if all the memory has been
       allocated. If it has not, there's no point in going further.
  ****/

  if ( (seminaive_naive_tablespace == NULL) ||
       (sigR == NULL) || (sigI == NULL) ||
       (so3Coef == NULL) ||
       (workspace1 == NULL) || (workspace2 == NULL) ||
       (workspace3 == NULL) ||
       (sigCoefR == NULL) || (sigCoefI == NULL) ||
       (patCoefR == NULL) || (patCoefI == NULL) ||
       (so3Sig == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }


  /* create plan for inverse SO(3) transform */
  n = 2 * bwOut ;
  howmany = n*n ;
  idist = n ;
  odist = n ;
  rank = 2 ;
  inembed[0] = n ;
  inembed[1] = n*n ;
  onembed[0] = n ;
  onembed[1] = n*n ;
  istride = 1 ;
  ostride = 1 ;
  na[0] = 1 ;
  na[1] = n ;

  p1 = fftw_plan_many_dft( rank, na, howmany,
			   workspace1, inembed,
			   istride, idist,
			   so3Sig, onembed,
			   ostride, odist,
			   FFTW_FORWARD, FFTW_ESTIMATE );


  fprintf(stdout,"Generating seminaive_naive tables...\n");
  seminaive_naive_table = SemiNaive_Naive_Pml_Table(bwIn, bwIn,
						    seminaive_naive_tablespace,
						    (double *) workspace2);

  n = 2 * bwIn ;

  printf("Reading in signal file\n");
  /* read in SIGNAL samples */
  /* first the signal */
  fp = fopen(argv[1],"r");
  for ( i = 0 ; i < 2*n * n ; i ++ )
    {
      fscanf(fp,"%lf", sigR + i);
      fscanf(fp,"%lf", sigI + i);
    }
  fclose( fp );

  printf("now taking spherical transform of signal\n");
  FST_semi_memo( sigR, sigI,
		 sigCoefR, sigCoefI,
		 n, seminaive_naive_table,
		 (double *) workspace2, 1, bwIn ) ;

  printf("Reading in pattern file\n");
  /* read in SIGNAL samples */
  /* first the signal */
  fp = fopen(argv[2],"r");
  for ( i = 0 ; i < 2 * n * n ; i ++ )
    {
      fscanf(fp,"%lf", sigR + i);
      fscanf(fp,"%lf", sigI + i);
    }
  fclose( fp );

  printf("now taking spherical transform of pattern\n");
  FST_semi_memo( sigR, sigI,
		 patCoefR, patCoefI,
		 n, seminaive_naive_table,
		 (double *) workspace2, 1, bwIn ) ;

  printf("freeing seminaive_naive_table and seminaive_naive_tablespace\n");
  
  free( seminaive_naive_table ) ;
  free( seminaive_naive_tablespace ) ;


  printf("about to combine coefficients\n");

  /* combine coefficients */
  tstart = csecond() ;
  so3CombineCoef_fftw( bwIn, bwOut, degLim,
		       sigCoefR, sigCoefI,
		       patCoefR, patCoefI,
		       so3Coef ) ;
  tstop = csecond();
  fprintf(stderr,"combine time \t = %.4e\n", tstop - tstart);
  
  printf("about to inverse so(3) transform\n");

  tstart = csecond();
  /* now inverse so(3) */
  Inverse_SO3_Naive_fftw( bwOut,
			  so3Coef,
			  so3Sig,
			  workspace1,
			  workspace2,
			  workspace3,
			  &p1,
			  1 ) ;
  tstop = csecond();
  printf("finished inverse so(3) transform\n");
  fprintf(stderr,"inverse so(3) time \t = %.4e\n", tstop - tstart);


  /* now find max value */
  maxval = 0.0 ;
  maxloc = 0 ;
  for ( i = 0 ; i < 8*bwOut*bwOut*bwOut ; i ++ )
    {
      if (so3Sig[i][0] >= maxval)
	{
	  maxval = so3Sig[i][0];
	  maxloc = i ;
	}
    }

  ii = floor( maxloc / (4.*bwOut*bwOut) );
  tmp = maxloc - (ii*4.*bwOut*bwOut);
  jj = floor( tmp / (2.*bwOut) );
  tmp = maxloc - (ii *4*bwOut*bwOut) - jj*(2*bwOut);
  kk = tmp ;

  printf("ii = %d\tjj = %d\tkk = %d\n", ii, jj, kk);

  printf("alpha = %f\nbeta = %f\ngamma = %f\n",
	 M_PI*jj/((REAL) bwOut),
	 M_PI*(2*ii+1)/(4.*bwOut),
	 M_PI*kk/((REAL) bwOut) );


  /* now save data -> just the real part because the
   imaginary parts should all be 0 ... right?!? */
  if ( argc == 7 )
    {
      printf("about to save data\n");
      fp = fopen( argv[6], "w" );
      for( i = 0 ; i < 8*bwOut*bwOut*bwOut ; i ++ )
	fprintf(fp,"%.16f\n", so3Sig[i][0]);
      fclose( fp );
    }

  fftw_free( so3Coef ) ;
  free( patCoefI );
  free( patCoefR );
  free( sigCoefI );
  free( sigCoefR );
  free( workspace3 );
  fftw_free( workspace2 );
  fftw_free( workspace1 );
  fftw_free( so3Sig ) ;
  free( sigI );
  free( sigR );

  return 0 ;

}
