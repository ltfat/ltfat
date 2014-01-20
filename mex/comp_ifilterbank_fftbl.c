#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 5
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define EXPORTALIAS comp_ifilterbank_fftbl

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include <math.h>
#include "config.h"

static LTFAT_FFTW(plan)* LTFAT_NAME(oldPlans) = 0;
static mwSize* LTFAT_NAME(oldLc) = 0;

static mwSize LTFAT_NAME(oldM) = 0;

void LTFAT_NAME(fftblMexAtExitFnc)()
{
   if(LTFAT_NAME(oldPlans)!=0)
   {
      for(mwIndex m=0;m<LTFAT_NAME(oldM);m++)
      {
         if(LTFAT_NAME(oldPlans)[m]!=0)
         {
            LTFAT_FFTW(destroy_plan)(LTFAT_NAME(oldPlans)[m]);
         }
      }
      ltfat_free(LTFAT_NAME(oldPlans));
      LTFAT_NAME(oldPlans) = 0;
   }

   if(LTFAT_NAME(oldLc)!=0)
   {
      ltfat_free(LTFAT_NAME(oldLc));
      LTFAT_NAME(oldLc) = 0;
   }
   LTFAT_NAME(oldM) = 0;
   #ifdef _DEBUG
   mexPrintf("Exit fnc called: %s\n",__PRETTY_FUNCTION__);
   #endif
}

// Calling convention:
// c = comp_ifilterbank_fftbl(c,G,foff,a,realonly)

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{

   static int atExitFncRegistered = 0;
   if(!atExitFncRegistered)
   {
      LTFAT_NAME(ltfatMexAtExit)(LTFAT_NAME(fftblMexAtExitFnc));
      atExitFncRegistered = 1;
   }

   const mxArray* mxc = prhs[0];
   const mxArray* mxG = prhs[1];
   double* foffDouble = (double*) mxGetData(prhs[2]);
   double* a = (double*) mxGetData(prhs[3]);
   double* realonlyDouble = (double*) mxGetData(prhs[4]);

   // number of channels
   mwSize W = mxGetN(mxGetCell(mxc,0));
   // filter number
   mwSize M = mxGetNumberOfElements(mxG);

   //
   mwSize acols = mxGetN(prhs[3]);

   double afrac[M];
   memcpy(afrac,a,M*sizeof(double));
   if(acols>1)
   {
      for(mwIndex m=0;m<M;m++)
      {
         afrac[m] = afrac[m]/a[M+m];
      }
   }

   // POINTER TO THE FILTERS
   LTFAT_COMPLEX* GPtrs[M];
   // input lengths
   mwSize inLen[M];
   mwSignedIndex foff[M];
   int realonly[M];
   // POINTER TO INPUTS
   LTFAT_COMPLEX* cPtrs[M];
   // filter lengths
   mwSize Gl[M];
   LTFAT_COMPLEX* cbuf[M];

   if(M!=LTFAT_NAME(oldM))
   {
      LTFAT_NAME(fftblMexAtExitFnc)();
      LTFAT_NAME(oldM) = M;
      LTFAT_NAME(oldLc) = (mwSize*) ltfat_calloc(M,sizeof(mwSize));
      LTFAT_NAME(oldPlans) = ltfat_calloc(M,sizeof(LTFAT_FFTW(plan)));
   }


   // over all channels
   for(mwIndex m =0; m<M; m++)
   {
      foff[m] = (mwSignedIndex) foffDouble[m];
      realonly[m] = (realonlyDouble[m]>1e-3);
      cPtrs[m] = (LTFAT_COMPLEX*) mxGetData(mxGetCell(mxc,m));
      inLen[m] = (mwSize) mxGetM(mxGetCell(mxc, m));
      GPtrs[m] = (LTFAT_COMPLEX*) mxGetData(mxGetCell(mxG, m));
      Gl[m] = (mwSize) mxGetNumberOfElements(mxGetCell(mxG, m));
      cbuf[m] = ltfat_malloc(inLen[m]*sizeof(LTFAT_COMPLEX));

      if(LTFAT_NAME(oldLc)[m]!=inLen[m])
      {
         LTFAT_NAME(oldLc)[m] = inLen[m];
         LTFAT_FFTW(plan) ptmp = LTFAT_FFTW(plan_dft_1d)(inLen[m],cbuf[m], cbuf[m], FFTW_FORWARD, FFTW_OPTITYPE);

         if(LTFAT_NAME(oldPlans)[m]!=0)
         {
            LTFAT_FFTW(destroy_plan)(LTFAT_NAME(oldPlans)[m]);
         }
         LTFAT_NAME(oldPlans)[m]=ptmp;
      }
   }

   // output data length
   mwSize L = (mwSize) floor(afrac[0]*inLen[0] + 0.5);
   plhs[0] = ltfatCreateMatrix(L, W, LTFAT_MX_CLASSID, mxCOMPLEX);
   mxArray* mxF = plhs[0];
   LTFAT_COMPLEX* FPtr = (LTFAT_COMPLEX*) mxGetData(mxF);


    LTFAT_NAME(ifilterbank_fftbl_plans)(cPtrs, GPtrs, L, Gl, W, afrac,
                                        M, foff, realonly, FPtr, LTFAT_NAME(oldPlans), cbuf);


    for(mwIndex m =0; m<M; m++)
    {
        ltfat_free(cbuf[m]);
    }


}
#endif




