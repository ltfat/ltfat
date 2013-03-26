#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include "mex.h"
#include "fftw3.h"

/* Function headers */
void TEMPLATE(preMexFn,LTFAT_REAL)(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]);
void TEMPLATE(postMexFn,LTFAT_REAL)(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]);
void TEMPLATE(split2combined,LTFAT_REAL)(const int L, const mxArray *parg, LTFAT_COMPLEX *outc);
void TEMPLATE(combined2split,LTFAT_REAL)(const int L, const LTFAT_COMPLEX *inc, mxArray *parg);

/* Casting functions */
//mxArray* TEMPLATE(mexSplit2combined(mxArray* prhsEl);


/* Handles data pre-formating */
void TEMPLATE(preMexFn,LTFAT_REAL)(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
#ifdef CHCOMPLEXFORMAT


#endif
}

/* Handles data post-formating */
void TEMPLATE(postMexFn,LTFAT_REAL)(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
#ifdef CHCOMPLEXFORMAT
// TO DO:

#endif
}

void TEMPLATE(split2combined,LTFAT_REAL)(const int L, const mxArray *parg, LTFAT_COMPLEX *outc)
{
   int ii;
   LTFAT_REAL *i_r, *i_i;


   i_r= (LTFAT_REAL*) mxGetPr(parg);

   if (mxIsComplex(parg))
   {
      i_i= (LTFAT_REAL*) mxGetPi(parg);

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

void TEMPLATE(combined2split,LTFAT_REAL)(const int L, const LTFAT_COMPLEX *inc, mxArray *parg)
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
