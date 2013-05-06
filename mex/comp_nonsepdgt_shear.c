#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 7
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS


int ltfat_round(double x)
{
  return (int)(x+.5);
}


#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"

// Calling convention:
// c=comp_nonsepdgt_shear(f,g,a,M,s0,s1,br);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  int a, M, N, L, W, s0, s1, br;
   mwSize ndim;
   mwSize dims[3];

   // Get matrix dimensions.
   L  = mxGetM(prhs[0]);
   W  = mxGetN(prhs[0]);

   a  = (int)mxGetScalar(prhs[2]);
   M  = (int)mxGetScalar(prhs[3]);
   s0 = (int)mxGetScalar(prhs[4]);
   s1 = (int)mxGetScalar(prhs[5]);
   br = (int)mxGetScalar(prhs[6]);

   N  = L/a;

   dims[0]=M;
   dims[1]=N;
   dims[2]=W;
   ndim=3;
   if (W==1)
   {
      ndim=2;
   }

   plhs[0] = ltfatCreateNdimArray(ndim,dims,LTFAT_MX_CLASSID,mxCOMPLEX);
   const LTFAT_REAL _Complex* f_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[0]);
   const LTFAT_REAL _Complex* g_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[1]);
   LTFAT_REAL _Complex* out_combined = (LTFAT_REAL _Complex*) mxGetData(plhs[0]);

   LTFAT_NAME(dgt_shear)(f_combined,g_combined,L,W,a,M,s0,s1,br,out_combined);

  return;
}
#endif

/*
#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

// Calling convention:
// c=comp_nonsepdgt_shear(f,g,a,M,s0,s1,br);


int ltfat_round(double x)
{
  return (int)(x+.5);
}

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
   int a, M, N, L, W, s0, s1, br;
   mwSize ndim;
   mwSize dims[3];

   ltfat_complex *f_combined, *g_combined, *out_combined;


   // Get matrix dimensions.
   L  = mxGetM(prhs[0]);
   W  = mxGetN(prhs[0]);

   a  = (int)mxGetScalar(prhs[2]);
   M  = (int)mxGetScalar(prhs[3]);
   s0 = (int)mxGetScalar(prhs[4]);
   s1 = (int)mxGetScalar(prhs[5]);
   br = (int)mxGetScalar(prhs[6]);

   N  = L/a;

   dims[0]=M;
   dims[1]=N;
   dims[2]=W;
   ndim=3;
   if (W==1)
   {
      ndim=2;
   }

   f_combined=(ltfat_complex*)mxMalloc(L*W*sizeof(ltfat_complex));
   g_combined=mxMalloc(L*sizeof(ltfat_complex));
   out_combined=(ltfat_complex*)mxMalloc(M*N*W*sizeof(ltfat_complex));

   split2combined(L*W, prhs[0], f_combined);
   split2combined(L, prhs[1], g_combined);

   dgt_shear((const ltfat_complex*)f_combined,
	     (const ltfat_complex*)g_combined,
	     L,W,a,M,s0,s1,br,
	     (ltfat_complex*)out_combined);

   mxFree(f_combined);
   mxFree(g_combined);

   plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);

   combined2split(M*N*W, (const ltfat_complex*)out_combined, plhs[0]);

   mxFree(out_combined);
   return;
}
*/
