#include "mex.h"

void comp_filterbank_td(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] );
void comp_filterbank_fft(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] );
void comp_filterbank_fftbl(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] );

// Calling convention:
//  comp_filterbank(f,g,a);

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
   const mxArray* mxf = prhs[0];
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
          if(acols==1&&L==mxGetNumberOfElement(mxGetField(gEl,0,"H")))
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


    if(fftCount>0)
    {

    }


}
