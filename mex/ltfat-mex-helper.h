#ifndef LTFAT_MEX_HELPER_H
#define LTFAT_MEX_HELPER_H 1

#include "ltfat.h"
#include "mex.h"

void split2combined(const int L, const mxArray *parg, ltfat_complex *outc)
{
   int ii;
   double *i_r, *i_i;
   

   i_r=mxGetPr(parg);

   if (mxIsComplex(parg))
   {
      i_i=mxGetPi(parg);
      
      for (ii=0;ii<L; ii++)
      {
	 outc[ii][0] = i_r[ii];
	 outc[ii][1] = i_i[ii];
      }      
   }
   else
   {
      /* No imaginary part */
      for (ii=0;ii<L; ii++)
      {
	 outc[ii][0] = i_r[ii];
	 outc[ii][1] = 0.0;
      }
   }
}

void combined2split(const int L, const ltfat_complex *inc, mxArray *parg)
{
   int ii;
   double *outr, *outi;
   
   outr=mxGetPr(parg);
   outi=mxGetPi(parg);

   for (ii=0;ii<L; ii++)
   {
      outr[ii] = inc[ii][0];
      outi[ii] = inc[ii][1];
   }
}

#endif
