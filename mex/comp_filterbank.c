#include "mex.h"
#include "fftw3.h"
#include <string.h>

void comp_filterbank_td(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] );
void comp_filterbank_fft(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] );
void comp_filterbank_fftbl(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] );

// Calling convention:
//  comp_filterbank(f,g,a);

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   const mxArray* mxf = prhs[0];
   mxArray* mxF = NULL;
   const mxArray* mxg = prhs[1];
   const mxArray* mxa = prhs[2];

   // input data length
   const mwSize L = mxGetM(mxf);
   const mwSize W = mxGetN(mxf);

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

   // Cell output
   plhs[0] = mxCreateCellMatrix(M, 1);

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
         c=comp_filterbank_td(f,g,a,offset,ext);
         BEWARE OF THE AUTOMATIC DEALLOCATION!! by the Matlab engine.
         Arrays can be very easily freed twice causing segfaults.
         This happends particulary when using mxCreateCell* which stores
         pointers to other mxArray structs. Setting all such pointers to
         NULL after they are used seems to solve it.
      */
      mxArray* plhs_td[1];
      const mxArray* prhs_td[5];
      prhs_td[0] = mxf;
      prhs_td[1] = mxCreateCellMatrix(tdCount,1);
      prhs_td[2] = mxCreateDoubleMatrix(tdCount,1,mxREAL);
      double* aPtr = (double*)mxGetPr(prhs_td[2]);
      prhs_td[3] = mxCreateDoubleMatrix(tdCount,1,mxREAL);
      double* offsetPtr = (double*)mxGetPr(prhs_td[3]);
      prhs_td[4] = mxCreateString("per");

      for(mwIndex m=0;m<tdCount;m++)
      {
          mxArray * gEl = mxGetCell(mxg, tdArgsIdx[m]);
          mxSetCell((mxArray*)prhs_td[1],m,mxGetField(gEl,0,"h"));
          // This has overhead
          //mxSetCell((mxArray*)prhs_td[1],m,mxDuplicateArray(mxGetField(gEl,0,"h")));

          aPtr[m] = a[tdArgsIdx[m]];
          offsetPtr[m] = mxGetScalar(mxGetField(gEl,0,"offset"));
      }

      // Finally call it!
      comp_filterbank_td(1,plhs_td,5, prhs_td);
      // This has overhead:
      // mexCallMATLAB(1,plhs_td,5, prhs_td,"comp_filterbank_td");

      // Copy pointers to a proper index in the output + unset all duplicate cell elements
      for(mwIndex m=0;m<tdCount;m++)
      {
          mxSetCell(plhs[0],tdArgsIdx[m],mxGetCell(plhs_td[0],m));
          mxSetCell(plhs_td[0],m,NULL);
          mxSetCell((mxArray*)prhs_td[1],m,NULL);
      }
    }


    if(fftCount>0 || fftblCount>0)
    {
        // Need to do FFT of mxf
        mwIndex ndim = 2;
        const mwSize dims[] = {L,W};

        if(mxIsDouble(mxf))
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

        fftw_plan p = fftw_plan_guru_split_dft(
          1, fftw_dims,
          1, howmanydims,
          mxGetPr(mxF), mxGetPi(mxF), mxGetPr(mxF), mxGetPi(mxF),
          FFTW_ESTIMATE);
          memcpy(mxGetPr(mxF),mxGetPr(mxf),L*W*sizeof(double));
          if(mxIsComplex(mxf))
            memcpy(mxGetPi(mxF),mxGetPi(mxf),L*W*sizeof(double));

          fftw_execute(p);
          fftw_destroy_plan(p);
        }
        else if(mxIsSingle(mxf))
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

          fftwf_plan p = fftwf_plan_guru_split_dft(
          1, fftw_dims,
          1, howmanydims,
          (float*)mxGetPr(mxF), (float*)mxGetPi(mxF),
          (float*) mxGetPr(mxF), (float*)mxGetPi(mxF),
          FFTW_ESTIMATE);
          memcpy(mxGetPr(mxF),mxGetPr(mxf),L*W*sizeof(float));
          if(mxIsComplex(mxf))
            memcpy(mxGetPi(mxF),mxGetPi(mxf),L*W*sizeof(float));

          fftwf_execute(p);
          fftwf_destroy_plan(p);
        }

    }

    if(fftCount>0)
    {
        mxArray* plhs_fft[1];
        const mxArray* prhs_fft[3];
        prhs_fft[0] = mxF;
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

        comp_filterbank_fft(1,plhs_fft,3, prhs_fft);

        for(mwIndex m=0;m<fftCount;m++)
        {
          mxSetCell(plhs[0],fftArgsIdx[m],mxGetCell(plhs_fft[0],m));
          mxSetCell(plhs_fft[0],m,NULL);
          mxSetCell((mxArray*)prhs_fft[1],m,NULL);
        }

    }

    if(fftblCount>0)
    {
        mxArray* plhs_fftbl[1];
        const mxArray* prhs_fftbl[5];
        prhs_fftbl[0] = mxF;
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

        comp_filterbank_fftbl(1,plhs_fftbl,5, prhs_fftbl);

        for(mwIndex m=0;m<fftblCount;m++)
        {
          mxSetCell(plhs[0],fftblArgsIdx[m],mxGetCell(plhs_fftbl[0],m));
          mxSetCell(plhs_fftbl[0],m,NULL);
          mxSetCell((mxArray*)prhs_fftbl[1],m,NULL);
        }

    }


}
