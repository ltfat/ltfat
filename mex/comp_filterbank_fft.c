#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 3
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define EXPORTALIAS comp_filterbank_fft

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include "math.h"
#include "config.h"

static LTFAT_FFTW(plan)** LTFAT_NAME(oldPlans) = 0;
static mwSize* LTFAT_NAME(oldLc) = 0;
static mwSize LTFAT_NAME(oldM) = 0;

// Calling convention:
// c = comp_filterbank_fft(F,G,a)

void LTFAT_NAME(fftMexAtExitFnc)()
{
   #ifdef _DEBUG
   mexPrintf("Exit fnc called: %s\n",__PRETTY_FUNCTION__);
   #endif
   if(LTFAT_NAME(oldPlans)!=0)
   {
      for(mwIndex m=0;m<LTFAT_NAME(oldM);m++)
	 {
            if(LTFAT_NAME(oldPlans)[m]!=0)
	    {
	       LTFAT_FFTW(destroy_plan)(*LTFAT_NAME(oldPlans)[m]);
	       ltfat_free(LTFAT_NAME(oldPlans)[m]);
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
        
}

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   static int atExitFncRegistered = 0;
   if(!atExitFncRegistered)
   {
      LTFAT_NAME(ltfatMexAtExit)(LTFAT_NAME(fftMexAtExitFnc));
      atExitFncRegistered = 1;
   }

   const mxArray* mxF = prhs[0];
   const mxArray* mxG = prhs[1];
   double* a = (double*) mxGetData(prhs[2]);

   // input data length
   mwSize L = mxGetM(mxF);
   // number of channels
   mwSize W = mxGetN(mxF);
   // filter number
   mwSize M = mxGetNumberOfElements(mxG);

   // output lengths
   mwSize outLen[M];
   // Filter pointer array
   LTFAT_REAL _Complex* GPtrs[M];
   // POINTER TO THE INPUT
   LTFAT_REAL _Complex* FPtr = (LTFAT_REAL _Complex*) mxGetData(prhs[0]);
   // POINTER TO OUTPUTS
   LTFAT_REAL _Complex* cPtrs[M]; // C99 feature
   plhs[0] = mxCreateCellMatrix(M, 1);

   if(M!=LTFAT_NAME(oldM))
   {
      LTFAT_NAME(fftMexAtExitFnc)();
      LTFAT_NAME(oldM) = M;
      LTFAT_NAME(oldLc) = (mwSize*) ltfat_calloc(M,sizeof(mwSize));
      LTFAT_NAME(oldPlans) = (LTFAT_FFTW(plan)**) ltfat_calloc(M,sizeof(LTFAT_FFTW(plan)*));
   }

   for(mwIndex m=0;m<M;++m)
   {
      outLen[m] = (mwSize) ceil( L/a[m] );
      GPtrs[m] = (LTFAT_REAL _Complex*) mxGetPr(mxGetCell(mxG, m));
      mxSetCell(plhs[0], m, ltfatCreateMatrix(outLen[m], W,LTFAT_MX_CLASSID,mxCOMPLEX));
      cPtrs[m] = (LTFAT_REAL _Complex*) mxGetData(mxGetCell(plhs[0],m));
      memset(cPtrs[m],0,outLen[m]*W*sizeof(LTFAT_REAL _Complex)); // Unnecessary memset to 0
       
      if(LTFAT_NAME(oldLc)[m]!=outLen[m])
      {
         LTFAT_NAME(oldLc)[m] = outLen[m];
	 LTFAT_FFTW(plan) ptmp = LTFAT_FFTW(plan_dft_1d)(outLen[m],(LTFAT_REAL (*)[2]) cPtrs[m],(LTFAT_REAL (*)[2]) cPtrs[m], FFTW_BACKWARD, FFTW_OPTITYPE);
		   
	 if(LTFAT_NAME(oldPlans)[m]!=0)
	 {
	    LTFAT_FFTW(destroy_plan)(*LTFAT_NAME(oldPlans)[m]);
	    ltfat_free(LTFAT_NAME(oldPlans)[m]);
	 }
	 LTFAT_NAME(oldPlans)[m] = ltfat_malloc(sizeof(ptmp));
	 memcpy(LTFAT_NAME(oldPlans)[m],&ptmp,sizeof(ptmp));
      }
   }

   // over all channels
   //  #pragma omp parallel for private(m)

   for(mwIndex m =0; m<M; m++)
   {
   // First col of cPtrs[m] is used as a buffer to assure correct memory aligment 
      for(mwIndex w =1; w<W; w++)
      {
         LTFAT_NAME(convsub_fft_plan)(FPtr+w*L,GPtrs[m],L,(size_t)a[m],cPtrs[m],LTFAT_NAME(oldPlans)[m]);
         memcpy(cPtrs[m] + w*outLen[m], cPtrs[m], outLen[m]*sizeof(LTFAT_REAL _Complex));
      }

      LTFAT_NAME(convsub_fft_plan)(FPtr,GPtrs[m],L,(size_t)a[m],cPtrs[m],LTFAT_NAME(oldPlans)[m]);
   }
}
#endif




