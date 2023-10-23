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

  test program to see if the FORWARD full SO3 transform
  actually works.

  I'll generate in Mma my sample values, and will see what
  coefficients I get at the end

  input: - bandwidth bw 
         - input file names of real and imaginary parts of
	   sample values ... should be (2*bw)^3 values each
	 - output files for real and imaginary parts of coefficients;
	   there are bw*(4*bw^2-1)/3 many coefficients
         - output format:
              0 -> ordered as the INVERSE TRANSFORM would expect the
                   coefficients
	      1 -> ordered (and labeled) by degree, i.e.

                      for l = 0 : bw - 1
		       for m1 = -l : l
                        for m2 = -l : l
                         coefficient of degree l, orders m1, m2

  example: test_soft_for 16 rSam.dat iSam.dat rCoeffs.dat iCoeffs.dat 1

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils_so3.h"
#include "soft.h"
#include "csecond.h"

int main( int argc,
	  char **argv )

{
  int i, bw, n, n3 ;
  int l, m1, m2, dummy, format ;
  double *rsignal, *isignal ;
  double *rcoeffs, *icoeffs ;
  double *workspace1, *workspace2 ;
  double tstart, tstop, runtime ;
  FILE *fp ;
  
  if (argc < 2)
    {
      fprintf(stdout,"Usage: test_soft_for bw realSam_file imagSam_file ");
      fprintf(stdout,"realCoef_file imagCoef_file order_flag\n");
      exit(0);
    }

  bw = atoi( argv[1] );
  n = 2 * bw ;
  n3 = n * n * n ;
  format = atoi( argv[6] );

  /* real and imaginary parts of signal each need n^3 space */
  rsignal = ( double * ) malloc( sizeof( double ) * n3 ) ;
  isignal = ( double * ) malloc( sizeof( double ) * n3 ) ;

  /* real and imaginary parts of coeffs each need
     totalCoeffs_so3( bw) amount of space */
  rcoeffs = ( double * ) malloc(sizeof( double ) * totalCoeffs_so3( bw ) ) ;
  icoeffs = ( double * ) malloc(sizeof( double ) * totalCoeffs_so3( bw ) ) ;

  /* now for LOTS OF workspace */
  workspace1 = ( double * ) malloc(sizeof( double ) * 4 * n3 ) ;
  workspace2 = ( double * ) malloc(sizeof( double ) * ( 24 * bw + 2 * bw * bw) );


  if ( ( rsignal == NULL ) || ( isignal == NULL ) ||
       ( rcoeffs == NULL ) || ( icoeffs == NULL ) ||
       ( workspace1 == NULL ) || ( workspace2 == NULL ) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }



  /* read in signal files */
  fp = fopen( argv[2] , "r" );
  for ( i = 0 ; i < n3 ; i ++ )
    fscanf(fp, "%lf", rsignal+i) ;
  fclose ( fp ) ;

  fp = fopen( argv[3] , "r" );
  for ( i = 0 ; i < n3 ; i ++ )
    fscanf(fp, "%lf", isignal+i) ;
  fclose ( fp ) ;


  /* turn on stopwatch */
  tstart = csecond( ) ;

  /* now do the forward transform */
  Forward_SO3_Naive( bw,
		     rsignal, isignal,
		     rcoeffs, icoeffs,
		     workspace1, workspace2 ) ;

  /* turn off stopwatch */
  tstop = csecond( ) ;
  runtime = tstop - tstart ;

  fprintf(stdout,"runtime: %.5f seconds\n", runtime);

  /* write out coefficients to disk */
  /* ordered as the inverse transform expects the coefficients */
  if ( format == 0 )
    {
      fp = fopen( argv[ 4 ], "w" );
      for ( i = 0 ; i < totalCoeffs_so3( bw ) ; i ++ )
	fprintf(fp, "%.16f\n", rcoeffs[ i ] ) ;
      fclose( fp ) ;
      
      fp = fopen( argv[ 5 ], "w" );
      for ( i = 0 ; i < totalCoeffs_so3( bw ) ; i ++ )
	fprintf(fp, "%.16f\n", icoeffs[ i ] ) ;
      fclose( fp ) ;
    }
  else /* ordered in a more human friendly way */
    {
      fp = fopen( argv[ 4 ], "w" );
      for ( l = 0 ; l < bw ; l ++ )
	for ( m1 = -l ; m1 < l + 1 ; m1 ++ )
	  for ( m2 = -l ; m2 < l + 1 ; m2 ++ )
	    {
	      dummy = so3CoefLoc( m1, m2, l, bw ) ;
	      fprintf(fp, "l = %d m1 = %d m2 = %d\t%.15f\n",
		      l, m1, m2, rcoeffs[ dummy ] ) ;
	    }
      fclose( fp ) ;

      fp = fopen( argv[ 5 ], "w" );
      for ( l = 0 ; l < bw ; l ++ )
	for ( m1 = -l ; m1 < l + 1 ; m1 ++ )
	  for ( m2 = -l ; m2 < l + 1 ; m2 ++ )
	    {
	      dummy = so3CoefLoc( m1, m2, l, bw ) ;
	      fprintf(fp, "l = %d m1 = %d m2 = %d\t%.15f\n",
		      l, m1, m2, icoeffs[ dummy ] ) ;
	    }
      fclose( fp );
    }

  /* free up memory (and there's lots of it) */
  free( workspace2 );
  free( workspace1 );
  free( icoeffs );
  free( rcoeffs );
  free( isignal );
  free( rsignal );

  return 0 ;
}
