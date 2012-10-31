#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"
#include "wavelets.h"

/*
Calling convention:
   c = comp_fwt_all( in, filts, J, type, ext)
*/
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 

	double * in; int inLen; int inChan;
	double ** out; int* outLen;
	double ** filts; int fLen;
    int ext = 0;
	int J;
	int noOfFilts = 2;

	// pointer to input data (column vector or matrix with input channels as collumns)
	in = mxGetPr(prhs[0]);
	// input length
	inLen = mxGetM(prhs[0]);
	// number of channels
	inChan = mxGetN(prhs[0]);
	// depth of decomposition
	J = (int) mxGetScalar(prhs[2]);

	// reading ext string
	int extStrLen = mxGetN(prhs[4]);
	char* extStr = (char*) mxMalloc(extStrLen);
	mxGetString(prhs[4],extStr,extStrLen+1);

	// reading flags string
	int flagsStrLen = mxGetN(prhs[3]);
	char* flagsStr = (char*) mxMalloc(flagsStrLen);
	mxGetString(prhs[3],flagsStr,flagsStrLen+1);
	
	// reading filters
	filts = (double **) mxMalloc(noOfFilts*sizeof(double *));
	fLen = mxGetNumberOfElements(mxGetCell(prhs[1], 0));
	for(int ii=0;ii<noOfFilts;ii++)
	    filts[ii] = mxGetPr(mxGetCell(prhs[1], ii));


	// prepare output poiner arrays 
	out = (double **) mxMalloc((J+1)*sizeof(double *));
	outLen = (int *) mxMalloc((J+1)*sizeof(int));
	plhs[0] = mxCreateCellMatrix(J+1, 1);
	
	bool doNonexp = false;

	// deciding which extension, 'per' extension also do nonexpansive representation
	if(!strcmp(extStr,"zpd")) // zero padding
		ext = 0;
	else if(!strcmp(extStr,"sym")) // half-point symmetric extension
		ext = 1;
	else if(!strcmp(extStr,"symw")) // whole-point symmetric extension
		ext = 2;
	else if(!strcmp(extStr,"asym")) // half-point anti-symmetric extension
		ext = 3;
	else if(!strcmp(extStr,"asymw")) // whole-point anti-symmetric extension
		ext = 4;
	else if(!strcmp(extStr,"ppd")) // periodic padding for expansive representation
		ext = 5;
	else if(!strcmp(extStr,"sp0")) // repeating last sample
		ext = 6;
	else
		doNonexp = true; // 'per' or otherwise: periodic extension with nonexpansive representation



	// decimated wavelet transform
	if(!strcmp(flagsStr,"dec"))
	{
      if(doNonexp)
	  {
		  // determine coefficient vectors lengths
		 int toNext = modPow2(pow2(J)-modPow2(inLen,pow2(J)),pow2(J));
		  for(int ii=1;ii<J+1;ii++)
	      {
	         outLen[ii] = ( inLen+pow2(J-ii+1)-1)/pow2(J-ii+1);
	      }
	      outLen[0] = outLen[1];

		  // do output arrays memory allocations
		  for(int ii=0;ii<J+1;ii++)
	      {
		      mxSetCell(plhs[0],ii, mxCreateDoubleMatrix(outLen[ii],inChan,mxREAL));
	      }


		for(int ch=0;ch<inChan;ch++)
		{
		   for(int ii=0;ii<J+1;ii++)
	       {
		      out[ii] = mxGetPr(mxGetCell(plhs[0], ii)) + ch*outLen[ii];
		   }
		      // call computation routine
		      dyadic_dwt_per(in+ch*inLen, inLen, out, outLen, filts, fLen, J);
		}
	  }
	  else
	  {
		  // determine coefficient vectors lengths
	      for(int ii=1;ii<J+1;ii++)
	      {
	         outLen[ii] = (int)( ((double)inLen)/pow2(J-ii+1) + (1.0-1.0/pow2(J-ii+1))*(fLen-1) );
	      }
	      outLen[0] = (int)( ((double)inLen)/pow2(J) + (1.0-1.0/pow2(J))*(fLen-1) );

		     // do output arrays memory allocations
	         for(int ii=0;ii<J+1;ii++)
	         {
		         mxSetCell(plhs[0],ii, mxCreateDoubleMatrix(outLen[ii],inChan,mxREAL));
	         }

			 // call computation routine
		     for(int ch=0;ch<inChan;ch++)
		     {
	            for(int ii=0;ii<J+1;ii++)
	            {
		          out[ii] = mxGetPr(mxGetCell(plhs[0], ii)) + ch*outLen[ii];
			    }
				dyadic_dwt_exp(in+ch*inLen, inLen, out, outLen, filts, fLen, J, ext);
		     }
	  }
	}
	// undecimated wavelet transform
	else if(!strcmp(flagsStr,"undec"))
	{
	   if(doNonexp)
	   {
		  // do output arrays memory allocations

		     for(int ii=0;ii<J+1;ii++)
	         {
		         mxSetCell(plhs[0],ii, mxCreateDoubleMatrix(inLen,inChan,mxREAL));
	         }

		  for(int ch=0;ch<inChan;ch++)
		  {
			 for(int ii=0;ii<J+1;ii++)
	         {
		         out[ii] = mxGetPr(mxGetCell(plhs[0], ii)) + ch*inLen;
	         }
             // call computation routine
		     undec_dwt_per(in+ch*inLen, inLen, out, inLen, filts, fLen, J);
		  }
	   }
	   else
	   {
		  // determine coefficient vectors lengths
		  outLen[J] = inLen + fLen -1;
	      for(int ii=J-1;ii>0;ii--)
	      {
			  outLen[ii] = outLen[ii+1] + pow2(J-ii)*fLen-(pow2(J-ii)-1) -1;
	      }
	      outLen[0] = outLen[1];

		  // do output arrays memory allocations

	         for(int ii=0;ii<J+1;ii++)
	         {
		         mxSetCell(plhs[0],ii, mxCreateDoubleMatrix(outLen[ii],inChan,mxREAL));
	         }

		  for(int ch=0;ch<inChan;ch++)
		  {
			 for(int ii=0;ii<J+1;ii++)
	         {
		         out[ii] = mxGetPr(mxGetCell(plhs[0], ii)) + ch*outLen[ii];
	         }
			 // call computation routine
		     undec_dwt_exp(in+ch*inLen, inLen, out, outLen, filts, fLen, J,ext);
		  }
	   } 
	
	
	}


	mxFree(filts);
	mxFree(out);
	mxFree(extStr);
	mxFree(flagsStr);
}


