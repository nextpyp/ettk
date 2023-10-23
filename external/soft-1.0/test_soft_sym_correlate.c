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


  - uses the Wigner-d symmetries
  - uses part of SpharmonicKit
  - INTERLEAVED (i.e. real/imaginary) COEFFICIENTS of signal and pattern files

    Note: the coefficients have to be in the order that SpharmonicKit
          spatial -> spectral S^2 Fourier transform routines produce.

  - [result] -> optional -> filename of all the correlation values
                (if you want all of them)
  - bw -> bw of input spherical signals
  - degLim -> max degree of Wigner-D functions you'll be using


  example: test_soft_sym_correlate signalCoefsFile patternCoefsFile bw degLim [result]


*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "complex.h" 
#include "csecond.h"
#include "so3_correlate_sym.h"
#include "soft_sym.h"


int main ( int argc,
	   char **argv )
{
  FILE *fp ;
  int i ;
  int bw, degLim ;
  double tstart, tstop ;
  REAL *rdata, *idata ;
  REAL *workspace1, *workspace2 ;
  REAL *sigCoefR, *sigCoefI ;
  REAL *patCoefR, *patCoefI ;
  REAL *so3CoefR, *so3CoefI ;
  int tmp, maxloc, ii, jj, kk ;
  REAL maxval ;


  if (argc < 5 )
    {
      printf("test_soft_sym_correlate sigCoef patCoef bw degLim [result]\n");
      exit(0) ;
    }

  bw = atoi( argv[3] );
  degLim = atoi( argv[4] );

  rdata = (REAL *) malloc( sizeof(REAL) * (8*bw*bw*bw) );
  idata = (REAL *) malloc( sizeof(REAL) * (8*bw*bw*bw) );
  workspace1 = (REAL *) malloc( sizeof(REAL) * (16*bw*bw*bw) );
  workspace2 = (REAL *) malloc( sizeof(REAL) * (24*bw + 2*bw*bw) );
  sigCoefR = (REAL *) malloc( sizeof(REAL) * bw * bw ) ;
  sigCoefI = (REAL *) malloc( sizeof(REAL) * bw * bw ) ;
  patCoefR = (REAL *) malloc( sizeof(REAL) * bw * bw ) ;
  patCoefI = (REAL *) malloc( sizeof(REAL) * bw * bw ) ;
  so3CoefR = (REAL *) malloc( sizeof(REAL) * ((4*bw*bw*bw-bw)/3) ) ;
  so3CoefI = (REAL *) malloc( sizeof(REAL) * ((4*bw*bw*bw-bw)/3) ) ;


  /****
       At this point, check to see if all the memory has been
       allocated. If it has not, there's no point in going further.
  ****/

  if ( (rdata == NULL) || (idata == NULL) ||
       (workspace1 == NULL) || (workspace2 == NULL) ||
       (sigCoefR == NULL) || (sigCoefI == NULL) ||
       (patCoefR == NULL) || (patCoefI == NULL) ||
       (so3CoefR == NULL) || (so3CoefI == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }


  /* read in coefficients */
  /* first the signal */
  fp = fopen(argv[1],"r");
  for ( i = 0 ; i < bw*bw ; i ++ )
    {
      fscanf(fp,"%lf", sigCoefR + i);
      fscanf(fp,"%lf", sigCoefI + i);
    }
  fclose( fp );

  /* now the pattern */
  fp = fopen(argv[2],"r");
  for ( i = 0 ; i < bw*bw ; i ++ )
    {
      fscanf(fp,"%lf", patCoefR + i );
      fscanf(fp,"%lf", patCoefI + i );
    }
  fclose( fp );

  printf("about to combine coefficients\n");

  /* combine coefficients */
  tstart = csecond() ;
  so3CombineCoef( bw, bw, degLim,
		  sigCoefR, sigCoefI,
		  patCoefR, patCoefI,
		  so3CoefR, so3CoefI ) ;
  tstop = csecond();
  fprintf(stderr,"combine time \t = %.4e\n", tstop - tstart);

  printf("about to inverse so(3) transform\n");

  tstart = csecond();
  /* now inverse so(3) */
  Inverse_SO3_Naive_sym( bw,
			 so3CoefR, so3CoefI,
			 rdata, idata,
			 workspace1, workspace2,
			 1 );
  tstop = csecond();
  printf("finished inverse so(3) transform\n");
  fprintf(stderr,"inverse so(3) time \t = %.4e\n", tstop - tstart);


  /* now find max value */
  maxval = 0.0 ;
  maxloc = 0 ;
  for ( i = 0 ; i < 8*bw*bw*bw ; i ++ )
    {
      if (rdata[i] >= maxval)
	{
	  maxval = rdata[i];
	  maxloc = i ;
	}
    }

  ii = floor( maxloc / (4.*bw*bw) );
  tmp = maxloc - (ii*4.*bw*bw);
  jj = floor( tmp / (2.*bw) );
  tmp = maxloc - (ii *4*bw*bw) - jj*(2*bw);
  kk = tmp ;

  printf("ii = %d\tjj = %d\tkk = %d\n", ii, jj, kk);

  printf("alpha = %f\nbeta = %f\ngamma = %f\n",
	 M_PI*jj/((REAL) bw),
	 M_PI*(2*ii+1)/(4.*bw),
	 M_PI*kk/((REAL) bw) );

  /* now save data -> just the real part because the
     imaginary parts should all be 0 */
  if ( argc == 6 )
    {
      printf("about to save data\n");
      fp = fopen( argv[5], "w" );
      for( i = 0 ; i < 8*bw*bw*bw ; i ++ )
	fprintf(fp,"%.16f\n", rdata[i]);
      fclose( fp );
    }


  free( so3CoefI );
  free( so3CoefR );
  free( patCoefI );
  free( patCoefR );
  free( sigCoefI );
  free( sigCoefR );
  free( workspace2 );
  free( workspace1 );
  free( idata );
  free( rdata );

  return 0 ;

}
