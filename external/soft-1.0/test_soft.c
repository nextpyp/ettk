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

  test program to loop through inverse-forward SO3 transforms
  lots of times to run error checks.

  - this is a plain vanilla so(3) fourier routine - not using
  FFTW and not using Wigner-d symmetries

  spectral - spatial - spectral

  input: - bandwidth bw
         - loops
  	 - output files for real and imaginary parts of errors


  example: test_soft bw loops realError imagError

  example: test_soft 16 10 rError.dat iError.dat

*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "utils_so3.h"
#include "soft.h"
#include "csecond.h"

int main( int argc,
	  char **argv )

{
  int l, k, j, bw, n, n3 ;
  int loops ;
  double *rsignal, *isignal ;
  double *rcoeffsIn, *icoeffsIn ;
  double *rcoeffsOut, *icoeffsOut ;
  double *workspace1, *workspace2 ;
  double tstartF, tstopF, runtimeF ;
  double tstartI, tstopI, runtimeI ;
  double total_time ;
  double ave_error;
  double ave_relerror;
  double stddev_error, stddev_relerror;
  double *relerror;
  double *curmax ;
  double granderror, grandrelerror;

  double realtmp, imagtmp ;
  double origmag, tmpmag ;

  long int seed ;
  FILE *fp ;
  
  if (argc < 2)
    {
      fprintf(stdout, "Usage: test_soft bw loops ");
      fprintf(stdout, "[realError_file imagError_file]\n");
      exit(0);
    }

  bw = atoi( argv[1] );
  loops = atoi( argv[2] );
  n = 2 * bw ;
  n3 = n * n * n ;

  /* real and imaginary parts of signal each need n^3 space */
  rsignal = ( double * ) malloc( sizeof( double ) * n3 ) ;
  isignal = ( double * ) malloc( sizeof( double ) * n3 ) ;

  /* real and imaginary parts of coeffs each need
     totalCoeffs_so3( bw) amount of space */
  rcoeffsIn = ( double * ) malloc(sizeof( double ) * totalCoeffs_so3( bw ) ) ;
  icoeffsIn = ( double * ) malloc(sizeof( double ) * totalCoeffs_so3( bw ) ) ;

  rcoeffsOut = ( double * ) malloc(sizeof( double ) * totalCoeffs_so3( bw ) ) ;
  icoeffsOut = ( double * ) malloc(sizeof( double ) * totalCoeffs_so3( bw ) ) ;

  /* now for LOTS OF workspace */
  workspace1 = ( double * ) malloc(sizeof( double ) * 4 * n3 ) ;
  workspace2 = ( double * ) malloc(sizeof( double ) * ( 24 * bw + 2 * bw * bw) );


  /** space for errors **/
  relerror = (double *) malloc(sizeof(double) * loops);
  curmax = (double *) malloc(sizeof(double) * loops);

  /* check if any problems allocating memory */
  if ( ( rsignal == NULL) || ( isignal == NULL ) ||
       ( rcoeffsIn == NULL ) || ( icoeffsIn == NULL ) ||
       ( rcoeffsOut == NULL ) || ( icoeffsOut == NULL ) ||
       ( workspace1 == NULL ) || ( workspace2 == NULL ) ||
       ( relerror == NULL ) || ( curmax == NULL ) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }


  /* generate seed for random number generator */
  time ( &seed ) ;
  srand48( seed ) ;

  /* initialize error */
  granderror = 0.0 ;
  grandrelerror = 0.0 ;

  /* initialize time */
  runtimeF = 0.0 ;
  runtimeI = 0.0 ;

  fprintf(stderr,"About to enter for-loop\n");

  for( k = 0 ; k < loops ; k ++ )
    {
      /* generate random coefficients */
      for( l = 0 ; l < totalCoeffs_so3( bw ) ; l++ )
	{
	  rcoeffsIn[ l ] = 2.0 * ( drand48() - 0.5 ) ;
	  icoeffsIn[ l ] = 2.0 * ( drand48() - 0.5 ) ;
	}

      /* turn on stopwatch */
      tstartI = csecond( ) ;

      /* now do inverse transform */
      Inverse_SO3_Naive( bw,
			 rcoeffsIn, icoeffsIn,
			 rsignal, isignal,
			 workspace1, workspace2 ) ;

      /* turn off stopwatch */
      tstopI = csecond( ) ;
      runtimeI += tstopI - tstartI ; 
      fprintf(stderr,"inv time \t = %.4e\n", tstopI - tstartI);
 
      /* turn on stopwatch */
      tstartF = csecond( ) ;
      
      /* now do the forward transform */
      Forward_SO3_Naive( bw,
			 rsignal, isignal,
			 rcoeffsOut, icoeffsOut,
			 workspace1, workspace2 ) ;

      /* turn off stopwatch */
      tstopF = csecond( ) ;
      runtimeF += tstopF - tstartF ;
      fprintf(stderr,"for time \t = %.4e\n", tstopF - tstartF);

      relerror[ k ] = 0.0 ;
      curmax[ k ] = 0.0 ;
      /* now figure out errors */
      for( j = 0 ; j < totalCoeffs_so3( bw ) ; j ++ )
	{
	  realtmp = rcoeffsIn[ j ] - rcoeffsOut[ j ] ;
	  imagtmp = icoeffsIn[ j ] - icoeffsOut[ j ] ;
	  origmag = sqrt((rcoeffsIn[ j ]*rcoeffsIn[ j ]) +
			 (icoeffsIn[ j ]*icoeffsIn[ j ]));
	  tmpmag = sqrt((realtmp*realtmp) +
			(imagtmp*imagtmp));
	  relerror[ k ] = MAX(relerror[ k ] , tmpmag/(origmag + pow(10.0, -50.0)));
	  curmax[ k ] = MAX( curmax[ k ], tmpmag );
	}

    
      fprintf(stderr,"r-o error\t = %.12f\n", curmax[ k ]);
      fprintf(stderr,"(r-o)/o error\t = %.12f\n\n", relerror[ k ]);
      
      granderror += curmax[ k ];
      grandrelerror += relerror[ k ];
    }


  total_time = runtimeF + runtimeI ;

  ave_error = granderror / ( (double) loops );
  ave_relerror = grandrelerror / ( (double) loops );
  stddev_error = 0.0 ; stddev_relerror = 0.0;
  for( k = 0 ; k < loops ; k ++ )
    {
      stddev_error += pow( ave_error - curmax[ k ] , 2.0 );
      stddev_relerror += pow( ave_relerror - relerror[ k ] , 2.0 );
    }
  /*** this won't work if loops == 1 ***/
  if( loops != 1 )
    {
      stddev_error = sqrt(stddev_error / ( (double) (loops - 1) ) );
      stddev_relerror = sqrt(stddev_relerror / ( (double) (loops - 1) ) );
    }


  fprintf(stderr,"Program: test_soft\n");
  fprintf(stderr,"Bandwidth = %d\n", bw);

#ifndef WALLCLOCK
  fprintf(stderr,"Total elapsed cpu time :\t\t %.4e seconds.\n",
	  total_time);
  fprintf(stderr,"Average cpu forward per iteration:\t %.4e seconds.\n",
	  runtimeF/((double) loops));  
  fprintf(stderr,"Average cpu inverse per iteration:\t %.4e seconds.\n",
	  runtimeI/((double) loops));
#else
  fprintf(stderr,"Total elapsed wall time :\t\t %.4e seconds.\n",
	  total_time);
  fprintf(stderr,"Average wall forward per iteration:\t %.4e seconds.\n",
	  runtimeF/((double) loops));  
  fprintf(stderr,"Average wall inverse per iteration:\t %.4e seconds.\n",
	  runtimeI/((double) loops));
#endif

  fprintf(stderr,"Average r-o error:\t\t %.4e\t",
	  granderror/((double) loops));
  fprintf(stderr,"std dev: %.4e\n",stddev_error);
  fprintf(stderr,"Average (r-o)/o error:\t\t %.4e\t",
	  grandrelerror/((double) loops));
  fprintf(stderr,"std dev: %.4e\n\n",stddev_relerror);


  /* if an error file is asked for (probably just after one loop) */
  if( argc > 3 )
    {
      fp = fopen( argv[ 3 ] , "w" );
      for( k = 0 ; k < totalCoeffs_so3( bw ) ; k++ )
	fprintf( fp, "%.15f\n", rcoeffsIn[ k ] - rcoeffsOut[ k ]);
      fclose( fp ) ;
      
      fp = fopen( argv[ 4 ] , "w" );
      for( k = 0 ; k < totalCoeffs_so3( bw ) ; k++ )
	fprintf( fp, "%.15f\n", icoeffsIn[ k ] - icoeffsOut[ k ]);
      fclose( fp ) ;      
    }

  /* free up memory (and there's lots of it) */
  free( curmax );
  free( relerror );
  free( workspace2 );
  free( workspace1 );
  free( icoeffsOut );
  free( rcoeffsOut );
  free( icoeffsIn );
  free( rcoeffsIn );
  free( isignal );
  free( rsignal );

  return 0 ;
}
