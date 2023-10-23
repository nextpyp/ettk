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

  a routine to see if the C function wigNaiveAnalysis actually works.

  for orders m1, m2, let l = max(|m1|, |m2|)
  let bw = bandwidth

  this routine should project a signal of length 2*bw onto
  all the Wigner functions

  d_{m1,m2}^l, d_{m1,m2}^{l+1}, ..., d_{m1,m2}^{bw-1}

  In other words, it will determine the coefficients
  < signal, d_{m1,m2}^{j} > for
  j = l ... bw - 1


  input: orders m1, m2
         bandwidth bw,
	 name of file containing sample values (length 2*bw)
	 name of output file (to write out the coefficients)


  example: test_Wigner_Analysis 3 4 16 signal.dat splat.dat


*/

#include <stdio.h>
#include <stdlib.h>
#include "utils_so3.h"
#include "makeWigner.h"
#include "wignerTransforms.h"


int main ( int argc ,
	   char **argv )
{
  int i, m1, m2, bw, n ;
  int m ;
  double *signal, *result, *wigners ;
  double *workspace, *scratch ;
  double *sinPts, *cosPts ;
  double *sinPts2, *cosPts2 ;
  FILE *fp ;

  
  if (argc < 6)
    {
      fprintf(stdout,"Usage: test_Wigner_Analysis m1 m2 bw input_file output_file\n");
      exit(0);
    }

  m1 = atoi( argv[1] );
  m2 = atoi( argv[2] );
  bw = atoi( argv[3] );

  m = MAX( ABS( m1 ) , ABS( m2 ) ) ;

  n = 2 * bw ;


  signal = ( double * ) malloc(sizeof( double ) * n ) ;
  result = ( double * ) malloc(sizeof( double ) * ( bw - m ) ) ;
  wigners = ( double * ) malloc( sizeof( double ) * ( bw - m ) * n ) ;
  workspace = (double *) malloc(sizeof( double ) * (4 + 9) * n ) ;
  sinPts = workspace ;
  cosPts = sinPts + n ;
  sinPts2 = cosPts + n ;
  cosPts2 = sinPts2 + n ;
  scratch = cosPts2 + n ; /* scratch needs to be of size n */

  /* note that the definition of wigSpec requires that instead of
     evaluating at beta, I need to evaluate at beta/2; ergo I call
     SinEvalPts2 instead of SinEvalPts, etc etc
  */


  /* read in sample values */
  fp = fopen( argv[4] , "r") ;
  for ( i = 0 ; i < n ; i ++ )
    fscanf( fp, "%lf", signal + i );
  fclose( fp ) ;


  /* precompute sines and cosines appropriate for making the
     wigners */

  SinEvalPts( n, sinPts ) ;
  CosEvalPts( n, cosPts ) ;
  SinEvalPts2( n, sinPts2 ) ;
  CosEvalPts2( n, cosPts2 ) ;

  /* now make the wigners */
  genWig_L2( m1, m2, bw,
	     sinPts, cosPts,
	     sinPts2, cosPts2,
	     wigners, scratch ) ;


  /* now analyze */
  wigNaiveAnalysis( m1, m2, bw, signal,
		    wigners, result,
		    scratch ) ;

  fp = fopen( argv[5], "w" );
  for ( i = 0 ; i <  (bw-m) ; i++ )
    fprintf(fp, "%.15f\n", result[i]);
  fclose( fp ) ;

  free( workspace ) ;
  free( wigners ) ;
  free( result ) ;
  free( signal ) ;

  return 0 ;
}
