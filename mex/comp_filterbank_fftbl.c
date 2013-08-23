#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 5
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define EXPORTALIAS comp_filterbank_fftbl
/*
static void (*exitFncPtr[])(void) = {NULL,NULL,NULL};

static bool exitFncRegistered = flase;

void fftblatAtExit(void)
{
   void (**exitFncPtrTmp)(void) = exitFncPtr;
   while(*exitFncPtrTmp!=NULL)
   {
      (**exitFncPtrTmp)();
	  exitFncPtrTmp++;
   }
}
*/
#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include "math.h"


static LTFAT_FFTW(plan)** LTFAT_NAME(oldPlans) = NULL;
static mwSize* LTFAT_NAME(oldLc) = NULL;
static mwSize* LTFAT_NAME(oldM) = NULL;

// Calling convention:
// c = comp_filterbank_fftbl(F,G,foff,a,realonly)

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  const mxArray* mxF = prhs[0];
  const mxArray* mxG = prhs[1];
  double* foff = (double*) mxGetData(prhs[2]);
  double* a = (double*) mxGetData(prhs[3]);
  double* realonly = (double*) mxGetData(prhs[4]);

  // input data length
  mwSize L = mxGetM(mxF);
  // number of channels
  mwSize W = mxGetN(mxF);
  // filter number
  mwSize M = mxGetNumberOfElements(mxG);
  


  //
  mwSize acols = mxGetN(prhs[3]);

  double* afrac = mxMalloc(M*sizeof(double));
  memcpy(afrac,a,M*sizeof(double));
  if(acols>1)
  {
      for(mwIndex m=0;m<M;m++)
      {
         afrac[m] = afrac[m]/a[M+m];
      }
   }


  // output lengths
  mwSize outLen[M];

  //mwSize* outLen = mxMalloc(M*sizeof(mwSize));
  

  for(unsigned int m = 0; m < M; m++)
  {
     outLen[m] = (mwSize) floor( L/afrac[m] +0.5);
  }


     // POINTER TO THE INPUT
     LTFAT_REAL _Complex* FPtr = (LTFAT_REAL _Complex*) mxGetData(prhs[0]);

     // POINTER TO THE FILTERS
     LTFAT_REAL _Complex* GPtrs[M];
     // filter lengths
     mwSize filtLen[M];
     // LTFAT_TYPE** gPtrs = (LTFAT_TYPE**) mxMalloc(M*sizeof(LTFAT_TYPE*));
     for(mwIndex m=0;m<M;m++)
     {
        GPtrs[m] = (LTFAT_REAL _Complex*) mxGetPr(mxGetCell(mxG, m));
        filtLen[m] = mxGetNumberOfElements(mxGetCell(mxG, m));
     }

     // POINTER TO OUTPUTS
     LTFAT_REAL _Complex* cPtrs[M]; // C99 feature
     //LTFAT_TYPE** cPtrs = (LTFAT_TYPE**) mxMalloc(M*sizeof(LTFAT_TYPE*));
     plhs[0] = mxCreateCellMatrix(M, 1);
     for(mwIndex m=0;m<M;++m)
     {
        mxSetCell(plhs[0], m, ltfatCreateMatrix(outLen[m], W,LTFAT_MX_CLASSID,mxCOMPLEX));
        cPtrs[m] = (LTFAT_REAL _Complex*) mxGetData(mxGetCell(plhs[0],m));
        memset(cPtrs[m],0,outLen[m]*W*sizeof(LTFAT_REAL _Complex));
		/*
		if(oldLc!=NULL && oldLc[m]!=outLen[m])
		{
		   oldLc[m] = outLen[m];
		   LTFAT_FFTW(plan) ptmp = LTFAT_FFTW(plan_dft_1d)(outLen[m], cPtrs[m], cPtrs[m], FFTW_BACKWARD, FFTW_ESTIMATE);
		   if(oldPlans[m]!=NULL)
		   {
		      LTFAT_FFTW(destroy_plan)(*oldPlans[m]);,
		   }
		   memcpy(oldPlans[m],&ptmp,sizeof(ptmp));
		}
		*/
     }

	 
     // over all channels
    #pragma omp parallel for 
        for(mwIndex m =0; m<M; m++)
        {
          for(mwIndex w =0; w<W; w++)
          {
           // Obtain pointer to w-th column in input
           LTFAT_REAL _Complex *FPtrCol = FPtr + w*L;
           // Obtaing pointer to w-th column in m-th element of output cell-array
           LTFAT_REAL _Complex *cPtrCol = cPtrs[m] + w*outLen[m];

           LTFAT_NAME(convsub_fftbl)(FPtrCol,L,GPtrs[m],filtLen[m],(int)foff[m],afrac[m],realonly[m],cPtrCol);

          }
        }
}
#endif




