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

  test program to see if the INVERSE full SO3 transform
  actually works.


  -uses Wigner-d symmetries

  I'll generate in Mma my sample values, and will see what
  coefficients I get at the end

  input: - bandwidth bw 
         - input file names of real and imaginary parts of
	   coefficients ... should be 1/3 bw (4 bw^2-1) values
	   each
	 - output files for real and imaginary parts of samples,
	   each of size (2 bw)^3

  example: test_soft_sym_inv.c 16 rCoeff.dat iCoeff.dat rSam.dat iSam.dat

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils_so3.h"
#include "soft_sym.h"
#include "csecond.h"

int main( int argc,
	  char **argv )

{
  int i, bw, n, n3 ;
  double *rsignal, *isignal ;
  double *rcoeffs, *icoeffs ;
  double *workspace1, *workspace2 ;
  double tstart, tstop, runtime ;
  FILE *fp ;
  
  if (argc < 2)
    {
      fprintf(stdout,"Usage: test_soft_sym_inv bw realCoeff imagCoeff ");
      fprintf(stdout,"realSam imagSam\n");
      exit(0);
    }

  bw = atoi( argv[1] );
  n = 2 * bw ;
  n3 = n * n * n ;

  /* real and imaginary parts of signal each need n^3 space */
  rsignal = ( double * ) calloc( n3, sizeof( double ) ) ;
  isignal = ( double * ) calloc( n3, sizeof( double ) ) ;

  /* real and imaginary parts of coeffs each need
     totalCoeffs_so3( bw) amount of space */
  rcoeffs = ( double * ) malloc(sizeof( double ) * totalCoeffs_so3( bw ) ) ;
  icoeffs = ( double * ) malloc(sizeof( double ) * totalCoeffs_so3( bw ) ) ;

  /* now for LOTS OF workspace */
  workspace1 = ( double * ) malloc(sizeof( double ) * 4 * n3 ) ;
  workspace2 = ( double * ) malloc(sizeof( double ) * ( 24 * bw + 2 * bw * bw) );


  /* check if allocated all memory without problems */
  if ( ( rsignal == NULL ) || ( isignal == NULL ) ||
       ( rcoeffs == NULL ) || ( icoeffs == NULL ) ||
       ( workspace1 == NULL ) || ( workspace2 == NULL ) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }

  /* read in coefficient files */
  fp = fopen( argv[2] , "r" );
  for ( i = 0 ; i < totalCoeffs_so3( bw ) ; i ++ )
    fscanf(fp, "%lf", rcoeffs+i) ;
  fclose ( fp ) ;

  fp = fopen( argv[3] , "r" );
  for ( i = 0 ; i < totalCoeffs_so3( bw ) ; i ++ )
    fscanf(fp, "%lf", icoeffs+i) ;
  fclose ( fp ) ;


  /* turn on stopwatch */
  tstart = csecond( ) ;

  /* now do the forward transform */
  Inverse_SO3_Naive_sym( bw,
			 rcoeffs, icoeffs,
			 rsignal, isignal,
			 workspace1, workspace2,
			 0 ) ;

  /* turn off stopwatch */
  tstop = csecond( ) ;
  runtime = tstop - tstart ;

  fprintf(stdout,"runtime: %.5f seconds\n", runtime);

  /* write out samples to disk */
  fp = fopen( argv[ 4 ], "w" );
  for ( i = 0 ; i < n3 ; i ++ )
    fprintf(fp, "%.16f\n", rsignal[ i ] ) ;
  fclose( fp ) ;

  fp = fopen( argv[ 5 ], "w" );
  for ( i = 0 ; i < n3 ; i ++ )
    fprintf(fp, "%.16f\n", isignal[ i ] ) ;
  fclose( fp ) ;


  /* free up memory (and there's lots of it) */
  free( workspace2 );
  free( workspace1 );
  free( icoeffs );
  free( rcoeffs );
  free( isignal );
  free( rsignal );

  return 0 ;
}
