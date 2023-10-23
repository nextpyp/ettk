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

  a routine to see if the C function genWig_L2 actually works.

  for orders m1, m2, let l = max(|m1|, |m2|)
  let bw = bandwidth

  this routine should generate all the Wigner functions

  d_{m1,m2}^l, d_{m1,m2}^{l+1}, ..., d_{m1,m2}^{bw-1}

  evaluated at the n = 2*bw many points

  pi*(2*[0..n-1] + 1)/(2 n)


  input: orders m1, m2
         bandwidth bw,
	 name of output file (to write function values)

  example: test_genWig m1 m2 bw outputFileName

  example: test_genWig 3 4 16 d34.dat


*/

#include <stdio.h>
#include <stdlib.h>
#include "utils_so3.h"
#include "makeWigner.h"

int main ( int argc ,
	   char **argv )
{
  int i, m1, m2, bw, n ;
  int m ;
  double *workspace, *scratch ;
  double *sinPts, *cosPts, *result ;
  double *sinPts2, *cosPts2 ;
  FILE *fp ;

  
  if (argc < 5)
    {
      fprintf(stdout,"Usage: test_genWig m1 m2 bw output_file_name\n");
      exit(0);
    }

  m1 = atoi( argv[1] );
  m2 = atoi( argv[2] );
  bw = atoi( argv[3] );

  m = MAX( ABS( m1 ) , ABS( m2 ) ) ;

  n = 2 * bw ;

  result = ( double * ) malloc(sizeof( double ) * n * ( bw - m ) ) ;
  workspace = (double *) malloc(sizeof( double ) * (4 + 8) * n ) ;
  sinPts = workspace ;
  cosPts = sinPts + n ;
  sinPts2 = cosPts + n ;
  cosPts2 = sinPts2 + n ;
  scratch = cosPts2 + n ; /* scratch needs to be of size 8*n */

  /* note that the definition of wigSpec requires that instead of
     evaluating at beta, I need to evaluate at beta/2; ergo I call
     SinEvalPts2 instead of SinEvalPts, etc etc
  */

  SinEvalPts( n, sinPts ) ;
  CosEvalPts( n, cosPts ) ;
  SinEvalPts2( n, sinPts2 ) ;
  CosEvalPts2( n, cosPts2 ) ;


  genWig_L2( m1, m2, bw,
	     sinPts, cosPts,
	     sinPts2, cosPts2,
	     result, scratch ) ;

  fp = fopen( argv[4], "w" );
  for ( i = 0 ; i < n*(bw-m) ; i++ )
    fprintf(fp, "%.15f\n", result[i]);
  fclose( fp ) ;

  free( workspace ) ;
  free( result ) ;

  return 0 ;
}
