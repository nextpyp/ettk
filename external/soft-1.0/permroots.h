
/***************************************************************************
  **************************************************************************
  
                Spherical Harmonic Transform Kit 2.6
  
   Sean Moore, Dennis Healy, Dan Rockmore, Peter Kostelec
   smoore@bbn.com, {healy,rockmore,geelong}@cs.dartmouth.edu
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
  
   Copyright 1997-2003  Sean Moore, Dennis Healy,
                        Dan Rockmore, Peter Kostelec
  
  
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


/* permroots.h - externally available structures and procedures */
/* These is a table of the roots of unity in bit-reversed order.
   These are used by the FFT evaluation and interpolation routines */
/* Size doesn't matter.  Here are the 4096 roots in bit-reverse order,
   but the first 2048 entries are the 2048 roots in bit-reverse order,
   and so on. */
/* This is an array of doubles in structure, in virtual reality they
   are complex numbers, where the first number is real part, second number is
   imaginary part */

#ifndef _PERMROOTS_H
#define _PERMROOTS_H

extern const double r4096[8192];

#endif /* _PERMROOTS_H */
