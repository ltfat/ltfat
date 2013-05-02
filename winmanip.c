#include <stdlib.h>
#include "config.h"
#include "ltfat.h"


/*  This routines changes the center of a vector from the beginning
 *  to the middle.
 *
 *  f  : Real valued input array.
 *  L  : Length of input array
 *  W  : Number of arrays to transform.
 *  h : Output, same size as input.
 *
 *  For the typical use of a single vector, set W=1.
 *
 *  This function works in the exact same way as the Matlab command
 *
 *  h = fftshift(f,1);
 *
 */

LTFAT_EXTERN void
LTFAT_NAME(fftshift_r)(const LTFAT_REAL *f, const int L, LTFAT_REAL *h)
{
  
   int ii;

   const div_t domod=div(L,2);

   for (ii=0; ii<domod.quot; ii++)
   {
      h[ii]=f[ii+domod.quot+domod.rem];
   }
   for (ii=0; ii<domod.quot+domod.rem; ii++)
   {
      h[ii+domod.quot]=f[ii];
   }
}


/*  This does the reverse of fftshift. When L is even, this function is
 * identical to fftshift, but for odd L there is a difference.
 *
 *  f  : Real valued input array.
 *  L  : Length of input array
 *  W  : Number of arrays to transform.
 *  h  : Output, same size as input.
 *
 *  For the typical use of a single vector, set W=1.
 *
 *  This function works in the exact same way as the Matlab command
 *
 *  h = ifftshift(f,1);
 *
 */
LTFAT_EXTERN void
LTFAT_NAME(ifftshift_r)(const LTFAT_REAL *f, const int L, LTFAT_REAL *h)
{
  
   int ii;
   div_t domod;

   domod=div(L,2);

   for (ii=0; ii<domod.quot+domod.rem; ii++)
   {
      h[ii]=f[ii+domod.quot];
   }
   for (ii=0; ii<domod.quot; ii++)
   {
      h[ii+domod.quot+domod.rem]=f[ii];
   }
}


/* This routine changes an FIR window to a LONG window.
 *
 * Input parameters:
 *  f     : Real valued input array.
 *  Lfir  : Length of input array
 *  Llong  : Length of output array
 *  h     : Output array
 */ 
LTFAT_EXTERN void
LTFAT_NAME(fir2long_r)(const LTFAT_REAL *f, const int Lfir, const int Llong,
	       LTFAT_REAL *h)
{
   const div_t domod=div(Lfir,2);
   
   /* ---- In the odd case, the additional element is kept in the first half. ---*/
     
   for (int ii=0; ii<domod.quot+domod.rem; ii++)
   {
      h[ii]=f[ii];
   }
   for (int ii=domod.quot+domod.rem; ii<Llong-domod.quot;ii++)
   {
      h[ii]=0.0;
   }
   const int ss=Llong-Lfir;
   for (int ii=domod.quot+domod.rem; ii<Lfir;ii++)
   {
      h[ii+ss]=f[ii];
   }
   
}

/* This routine changes an FIR window to a LONG window.
 *
 * Input parameters:
 *  f     : Complex valued input array.
 *  Lfir  : Length of input array
 *  Llong  : Length of output array
 *  h     : Output array
 */ 
LTFAT_EXTERN void
LTFAT_NAME(fir2long_c)(const LTFAT_COMPLEX *f, const int Lfir, const int Llong,
		       LTFAT_COMPLEX *h)
{
   const div_t domod=div(Lfir,2);
   
   /* ---- In the odd case, the additional element is kept in the first half. ---*/
   
   for (int ii=0; ii<domod.quot+domod.rem; ii++)
   {
      h[ii][0]=f[ii][0];
      h[ii][1]=f[ii][1];
   }
   for (int ii=domod.quot+domod.rem; ii<Llong-domod.quot;ii++)
   {
      h[ii][0]=0.0;
      h[ii][1]=0.0;
   }
   const int ss=Llong-Lfir;
   for (int ii=domod.quot+domod.rem; ii<Lfir;ii++)
   {
      h[ii+ss][0]=f[ii][0];
      h[ii+ss][1]=f[ii][1];
   }
   
}



/* This routine changes a LONG window to an FIR window.
 *
 * Input parameters:
 *  f     : Real valued input array.
 *  Llong  : Length of input array
 *  Lfir  : Length of output array
 *  h     : Output array
 */ 
LTFAT_EXTERN void
LTFAT_NAME(long2fir_r)(const LTFAT_REAL *f, const int Llong, const int Lfir, LTFAT_REAL *h)
{
  const div_t domod=div(Lfir,2);

  /* ---- In the odd case, the additional element is kept in the first half. ---*/
  
  for (int ii=0; ii<domod.quot+domod.rem; ii++)
  {
     h[ii]=f[ii];
  }
  const int ss=Llong-Lfir;
  for (int ii=domod.quot+domod.rem; ii<Lfir;ii++)
  {
     h[ii]=f[ii+ss];
  }
	
}



LTFAT_EXTERN void
LTFAT_NAME(long2fir_c)(const LTFAT_COMPLEX *f, const int Llong, const int Lfir, LTFAT_COMPLEX *h)
{
  const div_t domod=div(Lfir,2);

  /* ---- In the odd case, the additional element is kept in the first half. ---*/
     
  for (int ii=0; ii<domod.quot+domod.rem; ii++)
  {
     h[ii][0]=f[ii][0];
     h[ii][1]=f[ii][1];
  }
  const int ss=Llong-Lfir;
  for (int ii=domod.quot+domod.rem; ii<Lfir;ii++)
  {
     h[ii][0]=f[ii+ss][0];
     h[ii][1]=f[ii+ss][1];
  }
  	
}
