#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 5
#define TYPEDEPARGS 0
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
//  comp_nonsepwin2multi(g,a,M,lt,L);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   int a, M, L, Lg, lt1, lt2;
   double *lt;

   // Get matrix dimensions.
   Lg = mxGetM(prhs[0]);

   a=(int)mxGetScalar(prhs[1]);
   M=(int)mxGetScalar(prhs[2]);
   L=(int)mxGetScalar(prhs[4]);

   // Read the values of lt and round them to integers.
   lt = mxGetPr(prhs[3]);
   lt1 = ltfat_round(lt[0]);
   lt2 = ltfat_round(lt[1]);

   plhs[0] = ltfatCreateMatrix(L, lt2,LTFAT_MX_CLASSID,mxCOMPLEX);
   const LTFAT_REAL _Complex* g_combined = (const LTFAT_REAL _Complex*) mxGetData(prhs[0]);
   LTFAT_REAL _Complex* out_combined = (LTFAT_REAL _Complex*) mxGetData(plhs[0]);

   LTFAT_NAME(nonsepwin2multi)(g_combined,L,Lg,a,M,lt1,lt2,out_combined);

  return;
}
#endif

/*
#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

// Calling convention:
//  comp_nonsepwin2multi(g,a,M,lt,L);


int ltfat_round(double x)
{
  return (int)(x+.5);
}

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{

   int a, M, L, Lg, lt1, lt2;
   ltfat_complex *g_combined, *out_combined;
   double *lt;

   // Get matrix dimensions.
   Lg = mxGetM(prhs[0]);

   a=(int)mxGetScalar(prhs[1]);
   M=(int)mxGetScalar(prhs[2]);
   L=(int)mxGetScalar(prhs[4]);

   // Read the values of lt and round them to integers.
   lt = mxGetPr(prhs[3]);
   lt1 = ltfat_round(lt[0]);
   lt2 = ltfat_round(lt[1]);

   // Create temporary matrices to convert to correct complex layout.
   g_combined  =mxMalloc(Lg*sizeof(ltfat_complex));
   out_combined=mxMalloc(L*lt2*sizeof(ltfat_complex));

   split2combined(Lg, prhs[0], g_combined);

   nonsepwin2multi((const _Complex double*)g_combined,
		   L,Lg,a,M,lt1,lt2,
		   (_Complex double*)out_combined);

   mxFree(g_combined);

   plhs[0] = mxCreateDoubleMatrix(L, lt2, mxCOMPLEX);

   combined2split(L*lt2, (const ltfat_complex*)out_combined, plhs[0]);

   mxFree(out_combined);

   return;

}
*/
