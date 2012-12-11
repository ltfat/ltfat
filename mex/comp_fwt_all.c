#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"
#include "wavelets.h"

bool checkArrayEqual(double* in,int inLen)
{
for(int aIdx=0;aIdx<inLen;aIdx++)
	{
		if(in[aIdx]!=in[0])
		{
		   return false;
		}
	}
return true;
}

/*
Calling convention:
                      0      1  2  3     4    5
   c = comp_fwt_all( in, filts, J, a, type, ext)
*/
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{ 

	double * in; int inLen; int inChan;
	double ** out; int* outLen;
	double ** filts; int fLen;
	double *a; int aLen;
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
	a = mxGetPr(prhs[3]);
	aLen = mxGetNumberOfElements(prhs[3]);
	
	if(!checkArrayEqual(a,aLen))
	{
		mexErrMsgTxt("Not supported non-equal subsampling factors.");	
	}

	// reading ext string
	int extStrLen = mxGetN(prhs[5]);
	char* extStr = (char*) mxMalloc(extStrLen);
	mxGetString(prhs[5],extStr,extStrLen+1);

	// reading flags string
	int flagsStrLen = mxGetN(prhs[4]);
	char* flagsStr = (char*) mxMalloc(flagsStrLen);
	mxGetString(prhs[4],flagsStr,flagsStrLen+1);
	
	// reading filters
	noOfFilts = mxGetNumberOfElements(prhs[1]);
	filts = (double **) mxMalloc(noOfFilts*sizeof(double *));
	fLen = mxGetNumberOfElements(mxGetCell(prhs[1], 0));

	for(int ii=0;ii<noOfFilts;ii++)
	    filts[ii] = mxGetPr(mxGetCell(prhs[1], ii));


	// prepare output poiner arrays 
	int outArrays = J*(noOfFilts-1)+1;
	out = (double **) mxMalloc(outArrays*sizeof(double *));
	outLen = (int *) mxMalloc(outArrays*sizeof(int));
	plhs[0] = mxCreateCellMatrix(outArrays, 1);
	
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
	{
		ext = 7;
		doNonexp = true; // 'per' or otherwise: periodic extension with nonexpansive representation
	}

	
	
	int sub = a[0];
	// decimated wavelet transform
	if(!strcmp(flagsStr,"dec"))
	{
      if(doNonexp)
	  {
		  // determine coefficient vectors lengths
		  // int toNext = modPow2(pow2(J)-modPow2(inLen,pow2(J)),pow2(J));
		  for(int ii=0;ii<J;ii++)
	      {
			 for(int ff=0;ff<noOfFilts-1;ff++)
	         {
				int tmpDiv = ipow(a[ff+1],J-ii);
	            outLen[ii*(noOfFilts-1)+1+ff] = ( inLen+tmpDiv-1)/tmpDiv;
			 }
	      }
	      outLen[0] = outLen[1];

		  // do output arrays memory allocations
		  for(int ii=0;ii<outArrays;ii++)
	      {
		      mxSetCell(plhs[0],ii, mxCreateDoubleMatrix(outLen[ii],inChan,mxREAL));
	      }


		for(int ch=0;ch<inChan;ch++)
		{
		   for(int ii=0;ii<outArrays;ii++)
	       {
		      out[ii] = mxGetPr(mxGetCell(plhs[0], ii)) + ch*outLen[ii];
		   }
		      // call computation routine
		     dyadic_dwt_per(in+ch*inLen, inLen, out, outLen, filts, fLen,noOfFilts,sub ,J);
		}
	  }
	  else
	  {
		  // determine coefficient vectors lengths
		  int tmpcLen = inLen;
	      for(int ii=J;ii>=1;ii--)
	      {
			  int tmpOutLen = (int)( (tmpcLen+fLen-1 -1 +sub-1)/sub);
			  for(int ff=0;ff<noOfFilts-1;ff++)
	          {
				 outLen[(ii-1)*(noOfFilts-1)+1+ff] = tmpOutLen;
				 //int tmpDiv = ipow(a[ff+1],J-ii+1);
				 //outLen[(ii-1)*(noOfFilts-1)+1+ff] = (int)( ((double)inLen)/tmpDiv + (1.0-1.0/tmpDiv)*(fLen-1) );
	             // outLen[ii*(noOfFilts-1)+1+ff] = (int)( ((double)inLen)/pow2(J-ii+1) + (1.0-1.0/pow2(J-ii+1))*(fLen-1) );
			  }
			  tmpcLen = tmpOutLen;
	      }
	      outLen[0] = outLen[1];

		     // do output arrays memory allocations
	         for(int ii=0;ii<outArrays;ii++)
	         {
		         mxSetCell(plhs[0],ii, mxCreateDoubleMatrix(outLen[ii],inChan,mxREAL));
	         }

			 // call computation routine
		     for(int ch=0;ch<inChan;ch++)
		     {
	            for(int ii=0;ii<outArrays;ii++)
	            {
		          out[ii] = mxGetPr(mxGetCell(plhs[0], ii)) + ch*outLen[ii];
			    }
				dyadic_dwt_exp(in+ch*inLen, inLen, out, outLen, filts, fLen, noOfFilts, sub, J, ext);
		     }
	  }
	}
	// undecimated wavelet transform
	else if(!strcmp(flagsStr,"undec"))
	{
	   if(doNonexp)
	   {
		  // do output arrays memory allocations

		     for(int ii=0;ii<outArrays;ii++)
	         {
		         mxSetCell(plhs[0],ii, mxCreateDoubleMatrix(inLen,inChan,mxREAL));
	         }

		  for(int ch=0;ch<inChan;ch++)
		  {
			 for(int ii=0;ii<outArrays;ii++)
	         {
		         out[ii] = mxGetPr(mxGetCell(plhs[0], ii)) + ch*inLen;
	         }
             // call computation routine
		     undec_dwt_per(in+ch*inLen, inLen, out, inLen, filts, fLen, noOfFilts, sub, J);
		  }
	   }
	   else
	   {
		  // determine coefficient vectors lengths
		  int tmpcLen = inLen;
		  // outLen[outArrays-1] = inLen + fLen -1;
	      for(int ii=J;ii>0;ii--)
	      {
			 int tmpDiv = ipow(sub,J-ii)*fLen - ipow(sub,J-ii);
			 for(int ff=0;ff<noOfFilts-1;ff++)
	         {
			    //outLen[ii] = tmpcLen + pow2(J-ii)*fLen-(pow2(J-ii)-1) -1;
				outLen[(ii-1)*(noOfFilts-1)+1+ff] = tmpcLen + tmpDiv;
			 }
			 tmpcLen = tmpcLen + tmpDiv;
	      }
	      outLen[0] = outLen[1];

		  // do output arrays memory allocations

	         for(int ii=0;ii<outArrays;ii++)
	         {
		         mxSetCell(plhs[0],ii, mxCreateDoubleMatrix(outLen[ii],inChan,mxREAL));
	         }

		  for(int ch=0;ch<inChan;ch++)
		  {
			 for(int ii=0;ii<outArrays;ii++)
	         {
		         out[ii] = mxGetPr(mxGetCell(plhs[0], ii)) + ch*outLen[ii];
	         }
			 // call computation routine
		     undec_dwt_exp(in+ch*inLen, inLen, out, outLen, filts, fLen, noOfFilts, sub, J,ext);
		  }
	   } 
	
	
	}


	mxFree(filts);
	mxFree(out);
	mxFree(extStr);
	mxFree(flagsStr);
}


