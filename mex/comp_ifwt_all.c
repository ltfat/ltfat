#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"
#include "wavelets.h"

/*
Calling convention:
   f =  comp_ifwt_all(c, filts, J, Ls, flags,ext)
*/
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
	double ** in; int* inLen; int inChan;
	double * out; int outLen;
	double ** filts; int fLen;
    int ext = 0;
	int J;
	int noOfFilts = 2;
		J = (int) mxGetScalar(prhs[2]);

	// prepare output poiner arrays 
	inChan = mxGetN(mxGetCell(prhs[0], 0));
	in = (double **) mxMalloc((J+1)*inChan*sizeof(double *));
	inLen = (int *) mxMalloc((J+1)*sizeof(int));
	for(int ii=0;ii<J+1;ii++)
	    {
		   inLen[ii] = mxGetM(mxGetCell(prhs[0], ii)); 
	    }

	for(int ch=0;ch<inChan;ch++)
	{
	  for(int ii=0;ii<J+1;ii++)
	    {
	       in[ii+ch*(J+1)] = mxGetPr(mxGetCell(prhs[0], ii)) + ch*inLen[ii];
	    }
	 }


	filts = (double **) mxMalloc(noOfFilts*sizeof(double *));
	fLen = mxGetNumberOfElements(mxGetCell(prhs[1], 0));
	for(int ii=0;ii<noOfFilts;ii++)
	    filts[ii] = mxGetPr(mxGetCell(prhs[1], ii));



	// reading ext string
	int extStrLen = mxGetN(prhs[5]);
	char* extStr = (char*) mxMalloc(extStrLen);
	mxGetString(prhs[5],extStr,extStrLen+1);

	// reading flags string
	int flagsStrLen = mxGetN(prhs[4]);
	char* flagsStr = (char*) mxMalloc(flagsStrLen);
	mxGetString(prhs[4],flagsStr,flagsStrLen+1);

	// desired output signal length
	outLen = mxGetScalar(prhs[3]);
	plhs[0] = mxCreateDoubleMatrix(outLen, inChan, mxREAL);
	out = mxGetPr(plhs[0]);
	
	bool doNonexp = false;

	if(!strcmp(extStr,"zpd"))
		ext = 0;
	else if(!strcmp(extStr,"sym"))
		ext = 1;
	else if(!strcmp(extStr,"symw"))
		ext = 2;
	else if(!strcmp(extStr,"asym"))
		ext = 3;
	else if(!strcmp(extStr,"asymw"))
		ext = 4;
	else if(!strcmp(extStr,"ppd"))
		ext = 5;
	else if(!strcmp(extStr,"sp0"))
		ext = 6;
	else
		doNonexp = true;

	/* "per" for nonexpansive */


	if(!strcmp(flagsStr,"dec"))
	{
      if(doNonexp)
	  {
		  for(int ch=0;ch<inChan;ch++)
	      {
		    dyadic_idwt_per(&in[ch*(J+1)], inLen, out+ch*outLen, outLen, filts, fLen, J);
		  }
	    
	  }
	  else
	  {
		  for(int ch=0;ch<inChan;ch++)
	      {
			  dyadic_idwt_exp(&in[ch*(J+1)], inLen, out+ch*outLen, outLen, filts, fLen, J);
		  }
	  }
	}
	else if(!strcmp(flagsStr,"undec"))
	{
	   if(doNonexp)
	   {
		 for(int ch=0;ch<inChan;ch++)
	      {
		   undec_idwt_per(&in[ch*(J+1)], inLen[0], out+ch*outLen, outLen, filts, fLen, J);
		 }
	   }
	   else
	   {
		 for(int ch=0;ch<inChan;ch++)
	      {
		  undec_idwt_exp(&in[ch*(J+1)], inLen, out+ch*outLen, outLen, filts, fLen, J);
		 }
	   } 
	
	
	}


	mxFree(filts);
	mxFree(in);
	mxFree(extStr);
	mxFree(flagsStr);

}


