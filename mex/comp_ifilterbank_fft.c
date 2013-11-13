#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 3
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define EXPORTALIAS comp_ifilterbank_fft

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include "math.h"
#include "config.h"

// Calling convention:
// c = comp_ifilterbank_fft(c,G,a)

static LTFAT_FFTW(plan)** LTFAT_NAME(oldPlans) = 0;
static mwSize* LTFAT_NAME(oldLc) = 0;

static mwSize LTFAT_NAME(oldM) = 0;

void LTFAT_NAME(ifftMexAtExitFnc)()
{
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
  LTFAT_NAME(oldM) = 0;
  #ifdef _DEBUG
   mexPrintf("Exit fnc called: %s\n",__PRETTY_FUNCTION__);
  #endif
}

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  static int atExitFncRegistered = 0;
  if(!atExitFncRegistered)
  {
     LTFAT_NAME(ltfatMexAtExit)(LTFAT_NAME(ifftMexAtExitFnc));
     atExitFncRegistered = 1;
  }

  const mxArray* mxc = prhs[0];
  const mxArray* mxG = prhs[1];
  double* a = (double*) mxGetData(prhs[2]);

  // input data length
  mwSize L = mxGetM(mxGetCell(mxG,0));
  // number of channels
  mwSize W = mxGetN(mxGetCell(mxc,0));
  // filter number
  mwSize M = mxGetNumberOfElements(mxc);

  // input lengths
  mwSize inLen[M];
 
  plhs[0] = ltfatCreateMatrix(L, W,LTFAT_MX_CLASSID,mxCOMPLEX);
  memset(mxGetData(plhs[0]),0,L*W*sizeof(LTFAT_REAL _Complex));

  // POINTER TO THE OUTPUT
  LTFAT_REAL _Complex* FPtr = (LTFAT_REAL _Complex*) mxGetData(plhs[0]);


  // POINTERS TO THE FILTERS
  LTFAT_REAL _Complex* GPtrs[M];
	 
  // POINTER TO INPUTS
  LTFAT_REAL _Complex* cPtrs[M]; // C99 feature

  for(mwIndex m=0;m<M;m++)
  {
     GPtrs[m] = (LTFAT_REAL _Complex*) mxGetData(mxGetCell(mxG, m));
     cPtrs[m] = (LTFAT_REAL _Complex*) mxGetData(mxGetCell(mxc,m));
     inLen[m] = mxGetM(mxGetCell(mxc,m));
  }

  if(M!=LTFAT_NAME(oldM))
  {
     LTFAT_NAME(ifftMexAtExitFnc)();
     LTFAT_NAME(oldM) = M;
     LTFAT_NAME(oldLc) = (mwSize*) ltfat_calloc(M,sizeof(mwSize));
     LTFAT_NAME(oldPlans) = (LTFAT_FFTW(plan)**) ltfat_calloc(M,sizeof(LTFAT_FFTW(plan)*));
  }

  // over all channels
  //  #pragma omp parallel for private(m)

  for(mwIndex m =0; m<M; m++)
  {
     LTFAT_REAL _Complex* buffer = (LTFAT_REAL _Complex*) ltfat_malloc(inLen[m]*sizeof(LTFAT_REAL _Complex));

     if(LTFAT_NAME(oldLc)[m]!=inLen[m])
     {
        LTFAT_NAME(oldLc)[m] = inLen[m];
        LTFAT_FFTW(plan) ptmp = LTFAT_FFTW(plan_dft_1d)(inLen[m],(LTFAT_REAL (*)[2]) buffer,(LTFAT_REAL (*)[2]) buffer, FFTW_FORWARD, FFTW_OPTITYPE);

        if(LTFAT_NAME(oldPlans)[m]!=0)
        {
           LTFAT_FFTW(destroy_plan)(*LTFAT_NAME(oldPlans)[m]);
           ltfat_free(LTFAT_NAME(oldPlans)[m]);
        }
        LTFAT_NAME(oldPlans)[m] = ltfat_malloc(sizeof(ptmp));
        memcpy(LTFAT_NAME(oldPlans)[m],&ptmp,sizeof(ptmp));
     }
    for(mwIndex w =0; w<W; w++)
    {
       // Obtain pointer to w-th column in output
       LTFAT_REAL _Complex *FPtrCol = FPtr + w*L;
       // Obtaing pointer to w-th column in m-th element of input cell-array
       LTFAT_REAL _Complex *cPtrCol = cPtrs[m] + w*inLen[m];
          
       LTFAT_NAME(upconv_fft_plan)(cPtrCol,inLen[m],GPtrs[m],(size_t)a[m],FPtrCol,LTFAT_NAME(oldPlans)[m],buffer);
    }
    ltfat_free(buffer);
  }
}
#endif




