#include "fftw3.h"
#include <stdlib.h>
#include <math.h>

#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 2

#define TYPEDEPARGS 0
#define SINGLEARGS
#define REALARGS
#define COMPLEXARGS
#define NOCOMPLEXFMTCHANGE


static fftw_plan* p_old = 0;

static void fftrealAtExit(void)
{
  if(p_old!=0)
  {
     fftw_destroy_plan(*p_old);
     free(p_old);
  }
}


#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include "config.h"


// Calling convention:
//  comp_dct(f,type);
/*
FFTW_REDFT00 computes an REDFT00 transform, i.e. a DCT-I. (Logical N=2*(n-1), inverse is FFTW_REDFT00.)
FFTW_REDFT10 computes an REDFT10 transform, i.e. a DCT-II (sometimes called “the” DCT). (Logical N=2*n, inverse is FFTW_REDFT01.)
FFTW_REDFT01 computes an REDFT01 transform, i.e. a DCT-III (sometimes called “the” IDCT, being the inverse of DCT-II). (Logical N=2*n, inverse is FFTW_REDFT=10.)
FFTW_REDFT11 computes an REDFT11 transform, i.e. a DCT-IV. (Logical N=2*n, inverse is FFTW_REDFT11.)
*/


void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  // Register exit function only once
  #ifdef LTFAT_DOUBLE
  if(p_old==0)
  {
      mexAtExit(fftrealAtExit);
  }
  #endif

  mwIndex L, W, N, type;
  LTFAT_REAL *f_r, *f_i;
  LTFAT_FFTW(iodim) dims[1], howmanydims[1];
  LTFAT_FFTW(plan) p;
  LTFAT_FFTW(r2r_kind) kind[1];

  L = mxGetM(prhs[0]);
  W = mxGetN(prhs[0]);
  N = 2*L;
  type = (mwIndex) mxGetScalar(prhs[1]);

  // Copy inputs and get pointers
  if( mxIsComplex(prhs[0]))
  {
     plhs[0] = ltfatCreateMatrix(L, W, LTFAT_MX_CLASSID, mxCOMPLEX);
     f_i = (LTFAT_REAL*) mxGetPi(plhs[0]);
     memcpy(f_i,mxGetPi(prhs[0]),W*L*sizeof(LTFAT_REAL));
  }
  else
  {
     plhs[0] = ltfatCreateMatrix(L, W, LTFAT_MX_CLASSID, mxREAL);
  }
  f_r = (LTFAT_REAL*) mxGetPr(plhs[0]);
  memcpy(f_r,mxGetPr(prhs[0]),W*L*sizeof(LTFAT_REAL));

  // Create plan. Copy data from f to cout.
  dims[0].n = L;
  dims[0].is = 1;
  dims[0].os = 1;

  howmanydims[0].n = W;
  howmanydims[0].is = L;
  howmanydims[0].os = L;

  LTFAT_REAL sqrt2 = (LTFAT_REAL) sqrt(2.0);
  LTFAT_REAL postScale = (LTFAT_REAL) 1.0/sqrt2;
  LTFAT_REAL scale = (LTFAT_REAL) sqrt2*(1.0/(double)N)*sqrt((double)L);

  // Re-allocate and prescale input
  if(type==1||type==3)
  {
     for(mwIndex ii=0;ii<W;ii++)
     {
	    f_r[ii*L] *= sqrt2;
     }
     if(mxIsComplex(prhs[0]))
     {
        for(mwIndex ii=0;ii<W;ii++)
        {
           f_i[ii*L] *= sqrt2;
        }
     }
  }

  switch(type)
  {
	 case 1:
		N -= 2;
		for(mwIndex ii=0;ii<W;ii++)
        {
	       f_r[(ii+1)*L-1] *= sqrt2;
        }

        if(mxIsComplex(prhs[0]))
        {
           for(mwIndex ii=0;ii<W;ii++)
           {
              f_i[(ii+1)*L-1] *= sqrt2;
           }
        }

		scale = (LTFAT_REAL) sqrt2*(1.0/((double)N))*sqrt((double)L-1);
        kind[0] = FFTW_REDFT00;
     break;
	 case 2:
        kind[0] = FFTW_REDFT10;
     break;
	 case 3:
        kind[0] = FFTW_REDFT01;
     break;
	 case 4:
        kind[0] = FFTW_REDFT11;
     break;
     default:
		 mexErrMsgTxt("Unknown type.");
  }



  // The calling prototype
  //fftw_plan fftw_plan_guru_split_dft_r2c(
  //        int rank, const fftw_iodim *dims,
  //        int howmany_rank, const fftw_iodim *howmany_dims,
  //        double *in, double *ro, double *io,
  //        unsigned flags);


  p = LTFAT_FFTW(plan_guru_r2r)(1, dims,
				   1, howmanydims,
				   f_r, f_r,
				   kind,
				   FFTW_OPTITYPE);
  /*
  FFTW documentation qote http://www.fftw.org/fftw3_doc/New_002darray-Execute-Functions.html#New_002darray-Execute-Functions:
  ...
  creating a new plan is quick once one exists for a given size
  ...
  so why not to store the old plan..
  */


  if(p_old!=0)
  {
    fftw_destroy_plan(*p_old);
    free(p_old);
  }
  p_old = malloc(sizeof(p));
  memcpy(p_old,&p,sizeof(p));


  // Real FFT.
  LTFAT_FFTW(execute)(p);


  // Do the normalization
  for(int ii=0;ii<L*W;ii++)
  {
	  f_r[ii] *= scale;
  }

  if(type==1||type==2)
  {
     // Scale DC component
     for(int ii=0;ii<W;ii++)
     {
	    f_r[ii*L] *= postScale;
     }
  }

  if(type==1)
  {
     // Scale AC component
     for(int ii=0;ii<W;ii++)
     {
	    f_r[(ii+1)*L-1] *= postScale;
     }
  }

  // If the input is complex, process the imaginary part
  if(mxIsComplex(prhs[0]))
  {
	  LTFAT_FFTW(execute_r2r)(p,f_i,f_i);
	    // Do the normalization
      for(int ii=0;ii<L*W;ii++)
      {
	     f_i[ii] *= scale;
      }
      if(type==1||type==2)
      {
         // Scale DC component
         for(int ii=0;ii<W;ii++)
         {
	        f_i[ii*L] *= postScale;
         }
      }
      if(type==1)
      {
         // Scale AC component
         for(int ii=0;ii<W;ii++)
         {
	        f_i[(ii+1)*L-1] *= postScale;
         }
      }
  }

  // LTFAT_FFTW(destroy_plan)(p);

  return;
}
#endif


/*
#include "mex.h"
#include "config.h"
#include "fftw3.h"


// Calling convention:
//  comp_fftreal(f);


void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{

  int L, W, L2;
  double *f, *cout_r, *cout_i;
  fftw_iodim dims[1], howmanydims[1];
  fftw_plan p;

  L = mxGetM(prhs[0]);
  W = mxGetN(prhs[0]);

  L2 = (L/2)+1;

  // Get pointer to input.
  f=mxGetPr(prhs[0]);

  plhs[0] = mxCreateDoubleMatrix(L2, W, mxCOMPLEX);

  // Get pointer to output.
  cout_r = mxGetPr(plhs[0]);
  cout_i = mxGetPi(plhs[0]);

  // Create plan. Copy data from f to cout.
  dims[0].n = L;
  dims[0].is = 1;
  dims[0].os = 1;

  howmanydims[0].n = W;
  howmanydims[0].is = L;
  howmanydims[0].os = L2;

  // The calling prototype
  //fftw_plan fftw_plan_guru_split_dft_r2c(
  //        int rank, const fftw_iodim *dims,
  //        int howmany_rank, const fftw_iodim *howmany_dims,
  //        double *in, double *ro, double *io,
  //        unsigned flags);


  p = fftw_plan_guru_split_dft_r2c(1, dims,
				   1, howmanydims,
				   f, cout_r, cout_i, FFTW_OPTITYPE);

  // Real FFT.
  fftw_execute(p);

  fftw_destroy_plan(p);

  return;

}
*/
