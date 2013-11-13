#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 5
#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define EXPORTALIAS comp_filterbank_fftbl

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include <math.h>
#include "config.h"

static LTFAT_FFTW(plan)** LTFAT_NAME(oldPlans) = 0;
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


// Calling convention:
// c = comp_filterbank_fftbl(F,G,foff,a,realonly)

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{

  static int atExitFncRegistered = 0;
  if(!atExitFncRegistered)
  {
     LTFAT_NAME(ltfatMexAtExit)(LTFAT_NAME(fftblMexAtExitFnc));
     atExitFncRegistered = 1;
  }

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

  double afrac[M];
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
  // POINTER TO THE INPUT
  LTFAT_REAL _Complex* FPtr = (LTFAT_REAL _Complex*) mxGetData(prhs[0]);

  // POINTER TO THE FILTERS
  LTFAT_REAL _Complex* GPtrs[M];
  // filter lengths
  mwSize filtLen[M];
  // POINTER TO OUTPUTS
  LTFAT_REAL _Complex* cPtrs[M]; // C99 feature
  //LTFAT_TYPE** cPtrs = (LTFAT_TYPE**) mxMalloc(M*sizeof(LTFAT_TYPE*));
  plhs[0] = mxCreateCellMatrix(M, 1);
	 
  if(M!=LTFAT_NAME(oldM))
  {
     LTFAT_NAME(fftblMexAtExitFnc)();
     LTFAT_NAME(oldM) = M;
     LTFAT_NAME(oldLc) = (mwSize*) ltfat_calloc(M,sizeof(mwSize));
     LTFAT_NAME(oldPlans) = (LTFAT_FFTW(plan)**) ltfat_calloc(M,sizeof(LTFAT_FFTW(plan)*));
  }
	 
  for(mwIndex m=0;m<M;++m)
  {
     outLen[m] = (mwSize) floor( L/afrac[m] +0.5);
     GPtrs[m] = (LTFAT_REAL _Complex*) mxGetPr(mxGetCell(mxG, m));
     filtLen[m] = mxGetNumberOfElements(mxGetCell(mxG, m));
     mxSetCell(plhs[0], m, ltfatCreateMatrix(outLen[m], W,LTFAT_MX_CLASSID,mxCOMPLEX));
     cPtrs[m] = (LTFAT_REAL _Complex*) mxGetData(mxGetCell(plhs[0],m));
     memset(cPtrs[m],0,outLen[m]*W*sizeof(LTFAT_REAL _Complex));
		
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

	 /*
	 When W>1, the coefficients are stored as matrices in a cell array. The matrix itself is properly memory aligned, but the  
	 adresses of individial cols do not have to be. Since convsub_fftbl_plan uses new-array execute function http://www.fftw.org/doc/New_002darray-Execute-Functions.html
     all arrays have to be properly memory aligned. We go around it by using the first col as a buffer for all other cols.	 
	 
	 Maybe using the plan_many_dft could be used, but it is not clear whether it supports the new-array execution (since we store the plans between mex calls).
	 */

     // over all channels
  //  #pragma omp parallel for
   for(mwIndex m =0; m<M; m++)
   {
      for(mwIndex w =1; w<W; w++)
         {
	    // Using the first col of c as a temp array.
            LTFAT_NAME(convsub_fftbl_plan)(FPtr+w*L,L,GPtrs[m],filtLen[m],(int)foff[m],afrac[m],realonly[m],cPtrs[m],LTFAT_NAME(oldPlans)[m]);

	    // Copy to an appropriate position
	    memcpy(cPtrs[m] + w*outLen[m],cPtrs[m],outLen[m]*sizeof(LTFAT_REAL _Complex));
         }

	  // Working with the first col only.
	  LTFAT_NAME(convsub_fftbl_plan)(FPtr,L,GPtrs[m],filtLen[m],(ptrdiff_t)foff[m],afrac[m],realonly[m],cPtrs[m],LTFAT_NAME(oldPlans)[m]);
   }
}
#endif




