#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 2
#define TYPEDEPARGS 0
#define SINGLEARGS
#define NOCOMPLEXFMTCHANGE


#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include "config.h"

static LTFAT_FFTW(plan)* LTFAT_NAME(p_old) = 0;

void LTFAT_NAME(ifftrealAtExit)()
{
   if(LTFAT_NAME(p_old)!=0)
   {
     LTFAT_FFTW(destroy_plan)(*LTFAT_NAME(p_old));
     free(LTFAT_NAME(p_old));
   }
}


// Calling convention:
//  comp_ifftreal(f,N);

void LTFAT_NAME(ltfatMexFnc)( int UNUSED(nlhs), mxArray *plhs[],
                              int UNUSED(nrhs), const mxArray *prhs[] )
{
  static int atExitRegistered = 0;
  if(!atExitRegistered)
  {
      LTFAT_NAME(ltfatMexAtExit)(LTFAT_NAME(ifftrealAtExit));
      atExitRegistered = 1;
  }

  mwSize ii, L, W, L2;
  LTFAT_FFTW(plan) p;
  LTFAT_REAL *f, s;
  LTFAT_REAL *fin_r, *fin_i;


  L2 = (mwSize) mxGetM(prhs[0]);
  W  = (mwSize) mxGetN(prhs[0]);
  L  = (mwSize) mxGetScalar(prhs[1]);


  if(L/2+1!=L2)
    mexErrMsgTxt("IFFTREAL: Invalid output length.");

  fin_r = (LTFAT_REAL*)mxGetPr(prhs[0]);
  fin_i = (LTFAT_REAL*)mxGetPi(prhs[0]);
  // Case when input is real
  if(!mxIsComplex(prhs[0]))
  {
    mxArray* tmpIn = ltfatCreateMatrix(L2, W, LTFAT_MX_CLASSID , mxCOMPLEX);
    LTFAT_REAL *fin_r_old = (LTFAT_REAL*)mxGetPr(prhs[0]);
    fin_r = (LTFAT_REAL*)mxGetPr(tmpIn);
    fin_i = (LTFAT_REAL*)mxGetPi(tmpIn);
    for(mwIndex jj=0;jj<L2*W;jj++)
    {
        fin_r[jj]= fin_r_old[jj];
        fin_i[jj]= (LTFAT_REAL )0.0;
    }
  }

  // Create output and get pointer
  plhs[0] = ltfatCreateMatrix(L, W, LTFAT_MX_CLASSID , mxREAL);
  f= (LTFAT_REAL*) mxGetPr(plhs[0]);

  LTFAT_FFTW(iodim) dims[1], howmanydims[1];

  // Create plan. Copy data from cin to f.
  dims[0].n = L;
  dims[0].is = 1;
  dims[0].os = 1;

  howmanydims[0].n = W;
  howmanydims[0].is = L2;
  howmanydims[0].os = L;

  // The calling prototype

  // fftw_plan fftw_plan_guru_split_dft_c2r(
  //        int rank, const fftw_iodim *dims,
  //        int howmany_rank, const fftw_iodim *howmany_dims,
  //        double *ri, double *ii, double *out,
  //        unsigned flags);


  p = LTFAT_FFTW(plan_guru_split_dft_c2r)(1, dims,
				   1, howmanydims,
				   fin_r,fin_i,f, FFTW_ESTIMATE);

  LTFAT_NAME(ifftrealAtExit)();
  LTFAT_NAME(p_old) = malloc(sizeof(p));
  memcpy(LTFAT_NAME(p_old),&p,sizeof(p));


  // Real IFFT.
  LTFAT_FFTW(execute)(p);

  //LTFAT_FFTW(destroy_plan)(p);

  // Scale, because FFTW's normalization is different.
  s  = (LTFAT_REAL) (1.0/((LTFAT_REAL)L));
  for (ii=0; ii<L*W; ii++)
    {
      f[ii] *=s;
    }

  return;
}
#endif

