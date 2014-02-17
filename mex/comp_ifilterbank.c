#include "mex.h"
#include "fftw3.h"
#include <string.h>

/*
  Outside functions
*/
/*
void comp_ifilterbank_td(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] );
void comp_ifilterbank_fft(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] );
void comp_ifilterbank_fftbl(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] );

void comp_ifilterbank_fftbl_atexit();
void comp_ifilterbank_fft_atexit();
void comp_ifilterbank_td_atexit();
*/
/*
  Inside functions
*/
void ifilterbankAtExit();
void addToArray(mxArray* to,const mxArray* from);


static fftw_plan* p_double = NULL;
static fftwf_plan* p_float = NULL;

/*
Since the array is store for the lifetime of the MEX, we introduce limit od the array length.
2^20 ~ 16 MB of complex double
*/
#define MAXARRAYLEN 1048576
// Static pointer for holding the array the FFTW plan uses
static mxArray* mxF = NULL;

// Calling convention:
//  comp_ifilterbank(f,g,a,L);

/*
This will not add the imaginary part of from if to is real.
*/
void addToArray(mxArray* to, const mxArray* from)
{
#define ADDTOARRAY(to,from,L,T)   \
do{                               \
   for(mwIndex ii=0;ii<(L);ii++)    \
       ((to)[ii]) += (T) ((from)[ii]);       \
}while(0)

   mwSize nDimTO = mxGetNumberOfDimensions(to);
   mwSize nDimFROM = mxGetNumberOfDimensions(from);

   if(~mxIsComplex(to)&&mxIsComplex(from))
   {
      mexErrMsgTxt("COMP_IFILTERBANK: Cannot add complex to real.");
   }

   if(nDimFROM!=nDimTO)
   {
       mexErrMsgTxt("COMP_IFILTERBANK: Number of dimensions of arrays are not equal.");
   }

   const mwSize* dimsTO = mxGetDimensions(to);
   const mwSize* dimsFROM = mxGetDimensions(from);
   mwSize L = 1;

   for(mwIndex ii=0;ii<nDimTO;ii++)
   {
       if(dimsTO[ii]!=dimsFROM[ii])
       {
           mexErrMsgTxt("COMP_IFILTERBANK: Dimensions of arrays are not equal.");
       }
       L*=dimsTO[ii];
   }


if(mxIsDouble(to))
{
    if(mxIsDouble(from))
    {
        ADDTOARRAY(mxGetPr(to),mxGetPr(from),L,double);

        if(mxIsComplex(to)&&mxIsComplex(from))
        {
            ADDTOARRAY(mxGetPi(to),mxGetPi(from),L,double);
        }
    }
    else if(mxIsSingle(from))
    {
        ADDTOARRAY(mxGetPr(to),(float*)mxGetPr(from),L,float);
        if(mxIsComplex(to)&&mxIsComplex(from))
        {
            ADDTOARRAY(mxGetPi(to),(float*)mxGetPi(from),L,float);
        }
    }
    else
    {
        mexErrMsgTxt("COMP_IFILTERBANK: Unsupported type.");
    }
}
else
{
    if(mxIsDouble(from))
    {
        ADDTOARRAY((float*)mxGetPr(to),mxGetPr(from),L,float);

        if(mxIsComplex(to)&&mxIsComplex(from))
        {
            ADDTOARRAY((float*)mxGetPi(to),mxGetPi(from),L,float);
        }
    }
    else if(mxIsSingle(from))
    {
        ADDTOARRAY((float*)mxGetPr(to),(float*)mxGetPr(from),L,float);
        if(mxIsComplex(to)&&mxIsComplex(from))
        {
            ADDTOARRAY((float*)mxGetPi(to),(float*)mxGetPi(from),L,float);
        }
    }
    else
    {
        mexErrMsgTxt("COMP_IFILTERBANK: Unsupported type.");
    }
}
#undef ADDTOARRAY
}

/*
  MEX exit fnc. are not called by Matlab
*/
void ifilterbankAtExit()
{
   #ifdef _DEBUG
   mexPrintf("Exit fnc called: %s\n",__PRETTY_FUNCTION__);
   #endif
 /*
   comp_ifilterbank_fftbl_atexit();
   comp_ifilterbank_fft_atexit();
   comp_ifilterbank_td_atexit();
*/
   if(mxF!=NULL)
      mxDestroyArray(mxF);

   if(p_double!=NULL)
   {
       fftw_destroy_plan(*p_double);
       free(p_double);
   }

   if(p_float!=NULL)
   {
       fftwf_destroy_plan(*p_float);
       free(p_float);
   }

}

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   const mxArray* mxc = prhs[0];
   const mxArray* mxg = prhs[1];
   const mxArray* mxa = prhs[2];
   mxArray* tmpF = NULL;
   plhs[0] = NULL;

   // input data length
   const mwSize L = (mwSize) mxGetScalar(prhs[3]);
   const mwSize W = mxGetN(mxGetCell(mxc,0));

   // filter number
   const mwSize M = mxGetNumberOfElements(mxg);

   // a col count
   mwSize acols = mxGetN(mxa);

   // pointer to a
   double *a = (double*) mxGetData(mxa);


   if(acols>1)
   {
      int isOnes = 1;
      for(mwIndex m=0;m<M;m++)
      {
         isOnes = isOnes && a[M+m]==1;
      }

      if(isOnes)
      {
         acols = 1;
      }
   }

/*
   plhs[0] = mxCreateNumericArray(ndim,dims,mxGetClassId(mx),mxCOMPLEX);
   mxArray* mxf = plhs[0];
  */

   // Stuff for sorting the filters
   mwSize tdCount = 0;
   mwSize fftCount = 0;
   mwSize fftblCount = 0;
   mwIndex tdArgsIdx[M];
   mwIndex fftArgsIdx[M];
   mwIndex fftblArgsIdx[M];

   // WALK the filters to determine what has to be done
   for(mwIndex m =0; m<M; m++)
   {
      mxArray * gEl = mxGetCell(mxg, m);
      if(mxGetField(gEl,0,"h")!=NULL)
      {
          tdArgsIdx[tdCount++] = m;
          continue;
      }

      if(mxGetField(gEl,0,"H")!=NULL)
      {
          if(acols==1&&L==mxGetNumberOfElements(mxGetField(gEl,0,"H")))
          {
             fftArgsIdx[fftCount++] = m;
             continue;
          }
          else
          {
             fftblArgsIdx[fftblCount++] = m;
             continue;
          }
      }
   }

   if(tdCount>0)
   {
      /*
         Here, we have to reformat the inputs and pick up results to comply with:
         c=comp_ifilterbank_td(f,g,a,Ls,offset,ext);
         BEWARE OF THE AUTOMATIC DEALLOCATION!! by the Matlab engine.
         Arrays can be very easily freed twice causing segfaults.
         This happends particulary when using mxCreateCell* which stores
         pointers to other mxArray structs. Setting all such pointers to
         NULL after they are used seems to solve it.
      */
      mxArray* plhs_td[1];
      const mxArray* prhs_td[6];
      prhs_td[0] = mxc;
      prhs_td[1] = mxCreateCellMatrix(tdCount,1);
      prhs_td[2] = mxCreateDoubleMatrix(tdCount,1,mxREAL);
	  prhs_td[3] = mxCreateDoubleScalar(L);
      double* aPtr = (double*)mxGetPr(prhs_td[2]);
      prhs_td[4] = mxCreateDoubleMatrix(tdCount,1,mxREAL);
      double* offsetPtr = (double*)mxGetPr(prhs_td[4]);
      prhs_td[5] = mxCreateString("per");

      for(mwIndex m=0;m<tdCount;m++)
      {
          mxArray * gEl = mxGetCell(mxg, tdArgsIdx[m]);
          mxSetCell((mxArray*)prhs_td[1],m,mxGetField(gEl,0,"h"));
          aPtr[m] = a[tdArgsIdx[m]];
          offsetPtr[m] = mxGetScalar(mxGetField(gEl,0,"offset"));
      }

      // Finally call it!
     // comp_ifilterbank_td(1,plhs_td,6, prhs_td);
     mexCallMATLAB(1,plhs_td,6, prhs_td,"comp_ifilterbank_td");


      // Copy pointers to a proper index in the output + unset all duplicate cell elements
      for(mwIndex m=0;m<tdCount;m++)
      {
		  mxSetCell((mxArray*)prhs_td[1],m,NULL);
      }
	  // Copy pointer to output
	  plhs[0] = plhs_td[0];
      mxDestroyArray((mxArray*)prhs_td[1]);
      mxDestroyArray((mxArray*)prhs_td[2]);
      mxDestroyArray((mxArray*)prhs_td[3]);
      mxDestroyArray((mxArray*)prhs_td[4]);
	  mxDestroyArray((mxArray*)prhs_td[5]);

    }

    if(fftCount>0)
    {
        mxArray* plhs_fft[1];
        const mxArray* prhs_fft[3];
        prhs_fft[0] = mxc;
        prhs_fft[1] = mxCreateCellMatrix(fftCount,1);
        prhs_fft[2] = mxCreateDoubleMatrix(fftCount,1,mxREAL);
        double* aPtr = (double*)mxGetPr(prhs_fft[2]);

        for(mwIndex m=0;m<fftCount;m++)
        {
           mxArray * gEl = mxGetCell(mxg, fftArgsIdx[m]);
           mxSetCell((mxArray*)prhs_fft[1],m,mxGetField(gEl,0,"H"));
           // This has overhead
           //mxSetCell((mxArray*)prhs_td[1],m,mxDuplicateArray(mxGetField(gEl,0,"h")));
           aPtr[m] = a[fftArgsIdx[m]];
        }

       // comp_ifilterbank_fft(1,plhs_fft,3, prhs_fft);
         mexCallMATLAB(1,plhs_fft,3, prhs_fft,"comp_ifilterbank_fft");

        for(mwIndex m=0;m<fftCount;m++)
        {
          mxSetCell((mxArray*)prhs_fft[1],m,NULL);
        }
		tmpF = plhs_fft[0];
        mxDestroyArray((mxArray*)prhs_fft[1]);
        mxDestroyArray((mxArray*)prhs_fft[2]);
    }

    if(fftblCount>0)
    {
        mxArray* plhs_fftbl[1];
        const mxArray* prhs_fftbl[5];
        prhs_fftbl[0] = mxc;
        prhs_fftbl[1] = mxCreateCellMatrix(fftblCount,1);
        prhs_fftbl[2] = mxCreateDoubleMatrix(fftblCount,1,mxREAL);
        prhs_fftbl[3] = mxCreateDoubleMatrix(fftblCount,2,mxREAL);
        prhs_fftbl[4] = mxCreateDoubleMatrix(fftblCount,1,mxREAL);
        double* foffPtr = (double*)mxGetPr(prhs_fftbl[2]);
        double* aPtr = (double*)mxGetPr(prhs_fftbl[3]);
        double* realonlyPtr = (double*)mxGetPr(prhs_fftbl[4]);

        for(mwIndex m=0;m<fftblCount;m++)
        {
           mxArray * gEl = mxGetCell(mxg, fftblArgsIdx[m]);
           mxSetCell((mxArray*)prhs_fftbl[1],m,mxGetField(gEl,0,"H"));
           foffPtr[m] = mxGetScalar(mxGetField(gEl,0,"foff"));
           aPtr[m] = a[fftblArgsIdx[m]];

           if(acols>1)
              aPtr[m+fftblCount] = a[fftblArgsIdx[m]+M];
           else
              aPtr[m+fftblCount] = 1;

           realonlyPtr[m] = mxGetScalar(mxGetField(gEl,0,"realonly"));
        }


        //comp_ifilterbank_fftbl(1,plhs_fftbl,5, prhs_fftbl);
        mexCallMATLAB(1,plhs_fftbl,5, prhs_fftbl,"comp_ifilterbank_fftbl");

        for(mwIndex m=0;m<fftblCount;m++)
        {
          mxSetCell((mxArray*)prhs_fftbl[1],m,NULL);
        }

        if(tmpF==NULL)
        {
 		   tmpF = plhs_fftbl[0];
        }
        else
        {
           addToArray(tmpF,plhs_fftbl[0]);
        }


        mxDestroyArray((mxArray*)prhs_fftbl[1]);
        mxDestroyArray((mxArray*)prhs_fftbl[2]);
        mxDestroyArray((mxArray*)prhs_fftbl[3]);
        mxDestroyArray((mxArray*)prhs_fftbl[4]);
    }


	if(fftCount>0 || fftblCount>0)
    {
        // Need to do FFT of mxf
        mwIndex ndim = 2;
        const mwSize dims[] = {L,W};

        if(mxF==NULL || mxGetM(mxF)!=L || mxGetN(mxF)!=W || mxGetClassID(mxF)!=mxGetClassID(tmpF))
        {
            if(mxF!=NULL)
            {
               mxDestroyArray(mxF);
               mxF = NULL;
              // printf("Should be called just once\n");
            }


        if(mxIsDouble(tmpF))
        {
        mxF = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX);
        fftw_iodim fftw_dims[1];
        fftw_iodim howmanydims[1];

        fftw_dims[0].n = L;
        fftw_dims[0].is = 1;
        fftw_dims[0].os = 1;

        howmanydims[0].n = W;
        howmanydims[0].is = L;
        howmanydims[0].os = L;

        if(p_double==NULL)
            p_double = (fftw_plan*) malloc(sizeof(fftw_plan));
         else
            fftw_destroy_plan(*p_double);


        *p_double = fftw_plan_guru_split_dft(
          1, fftw_dims,
          1, howmanydims,
          mxGetPi(mxF), mxGetPr(mxF), mxGetPi(mxF), mxGetPr(mxF),
          FFTW_ESTIMATE);

        }
        else if(mxIsSingle(tmpF))
        {
          mxF = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxCOMPLEX);
         // mexPrintf("M= %i, N= %i\n",mxGetM(mxF),mxGetN(mxF));
          fftwf_iodim fftw_dims[1];
          fftwf_iodim howmanydims[1];

          fftw_dims[0].n = L;
          fftw_dims[0].is = 1;
          fftw_dims[0].os = 1;

          howmanydims[0].n = W;
          howmanydims[0].is = L;
          howmanydims[0].os = L;

          if(p_float==NULL)
             p_float = (fftwf_plan*) malloc(sizeof(fftwf_plan));
          else
             fftwf_destroy_plan(*p_float);

          *p_float = fftwf_plan_guru_split_dft(
          1, fftw_dims,
          1, howmanydims,
          (float*)mxGetPi(mxF), (float*)mxGetPr(mxF),
          (float*) mxGetPi(mxF), (float*)mxGetPr(mxF),
          FFTW_ESTIMATE);

        }


        }

        if(mxIsDouble(tmpF))
        {
		  double* mxFPtr = mxGetPr(mxF);
		  double* mxFPti = mxGetPi(mxF);
          memcpy(mxFPtr,mxGetPr(tmpF),L*W*sizeof(double));
          memcpy(mxFPti,mxGetPi(tmpF),L*W*sizeof(double));

          fftw_execute(*p_double);

		  for(mwIndex ii=0;ii<L*W;ii++)
		  {
		     mxFPtr[ii]/=(double) L;
			 mxFPti[ii]/=(double) L;
		  }
        }
        else if(mxIsSingle(tmpF))
        {
		  float* mxFPtr = (float*) mxGetPr(mxF);
		  float* mxFPti = (float*) mxGetPi(mxF);
          memcpy(mxFPtr,mxGetPr(tmpF),L*W*sizeof(float));
          memcpy(mxFPti,mxGetPi(tmpF),L*W*sizeof(float));

          fftwf_execute(*p_float);

		  for(mwIndex ii=0;ii<L*W;ii++)
		  {
		     mxFPtr[ii]/=(float) L;
			 mxFPti[ii]/=(float) L;
		  }
        }

        if(mxIsDouble(mxF))
        {
           memcpy(mxGetPr(tmpF),mxGetPr(mxF),L*W*sizeof(double));
           memcpy(mxGetPi(tmpF),mxGetPi(mxF),L*W*sizeof(double));
        }
        else if(mxIsSingle(mxF))
        {
           memcpy(mxGetPr(tmpF),mxGetPr(mxF),L*W*sizeof(float));
           memcpy(mxGetPi(tmpF),mxGetPi(mxF),L*W*sizeof(float));
        }

		if(plhs[0]!=NULL)
		{
            addToArray(tmpF,plhs[0]);
        }
        plhs[0] = tmpF;
    }


    /* This should overwrite function registered by mexAtExit in any of the previously
    called MEX files */
   mexAtExit(ifilterbankAtExit);

   if(mxF!=NULL)
      mexMakeArrayPersistent(mxF);

   if(L*W>MAXARRAYLEN && mxF!=NULL)
   {
       //printf("Damn. Should not get here\n");
       mxDestroyArray(mxF);
       mxF = NULL;
   }


//int prd = 0;
}
