#include <string.h>
#include <math.h>
#include "wavelets.h"


LTFAT_EXTERN
    void LTFAT_NAME(undec_dwt_per)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *out[], const int outLen, const LTFAT_REAL *filts[], int fLen, int noOfFilts, int subFac, const int J)
{
	/************************* STEP 1: Prepare ***************************/
	LTFAT_REAL **filtsNorm = (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	LTFAT_REAL **tmpOut = (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	for(int ff=0;ff<noOfFilts;ff++)
	{
	   filtsNorm[ff] = (LTFAT_REAL*) ltfat_malloc(fLen*sizeof(LTFAT_REAL));
	   for(int fid=0;fid<fLen;fid++)
	   {
		   filtsNorm[ff][fid] = (LTFAT_REAL) filts[ff][fid]/sqrt((LTFAT_REAL)subFac);
	   }
	}
	int noOfOutputs = J*(noOfFilts-1) + 1;
	/*********************************************************************/
	/************************* STEP 2: CALCULATIONS **********************/
	if(J<=1)
	{
	  int skip = (fLen+1)/2;
	  for(int ff=0;ff<noOfFilts;ff++)
	  {
		  tmpOut[ff] = out[ff];
	  }
	  LTFAT_NAME(conv_td_sub)(in,inLen,tmpOut,outLen,filtsNorm,fLen,noOfFilts,1,skip,5,0);
	}
	else
	{
	  // there is no other way: creating buffer to hold intermediate results
	  LTFAT_REAL *buffer = (LTFAT_REAL *) ltfat_malloc(outLen*sizeof(LTFAT_REAL));
	  tmpOut[0] = buffer;
	  for(int ff=0;ff<noOfFilts-1;ff++)
	  {
		  tmpOut[ff+1] = out[noOfOutputs-(noOfFilts-1)+ff];
	  }

	  int skip = (fLen+1)/2;
	  LTFAT_NAME(conv_td_sub)(in,inLen,tmpOut,outLen,filtsNorm,fLen,noOfFilts,1,skip,5,0);


	  for(int jj=1;jj<J-1;jj++)
	  {
		 skip = (ipow(subFac,jj)*fLen+1)/2;
	     for(int ff=0;ff<noOfFilts-1;ff++)
	     {
	        tmpOut[1+ff] = out[noOfOutputs-(jj+1)*(noOfFilts-1)+ff];
	     }
	     LTFAT_NAME(conv_td_sub)(buffer,outLen,tmpOut,outLen,filtsNorm,fLen,noOfFilts,1,skip,5,ipow(subFac,jj));
	  }
	
	  skip = (ipow(subFac,J-1)*fLen+1)/2;
	  tmpOut[0] = out[0]; 
	  for(int ff=0;ff<noOfFilts-1;ff++)
	  {
	     tmpOut[1+ff] = out[1+ff];
	  }
	  LTFAT_NAME(conv_td_sub)(buffer,outLen,tmpOut,outLen,filtsNorm,fLen,noOfFilts,1,skip,5,ipow(subFac,J-1));
    
	  /****** STEP 3: FREE ALL ******/
	  ltfat_free(buffer);
	}
    for(int ff=0;ff<noOfFilts;ff++)  ltfat_free(filtsNorm[ff]); 
	ltfat_free(filtsNorm);
	ltfat_free(tmpOut);
}

LTFAT_EXTERN
    void LTFAT_NAME(undec_dwt_exp)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *out[], const int outLen[], const LTFAT_REAL *filts[], int fLen, int noOfFilts, int subFac, const int J, int ext)
{
    /************************* STEP 1: Prepare ***************************/
	LTFAT_REAL **filtsNorm = (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	LTFAT_REAL **tmpOut = (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	for(int ff=0;ff<noOfFilts;ff++)
	{
	   filtsNorm[ff] = (LTFAT_REAL*) ltfat_malloc(fLen*sizeof(LTFAT_REAL));
	   for(int fid=0;fid<fLen;fid++)
	   {
		   filtsNorm[ff][fid] = (LTFAT_REAL) filts[ff][fid]/sqrt((LTFAT_REAL)subFac);
	   }
	}
	int noOfOutputs = J*(noOfFilts-1) + 1;
	/*********************************************************************/

	if(J<=1)
	{
	  int skip = 0;
	  for(int ff=0;ff<noOfFilts;ff++)
	  {
		  tmpOut[ff] = out[ff];
	  }
	  LTFAT_NAME(conv_td_sub)(in,inLen,tmpOut,outLen[0],filtsNorm,fLen,noOfFilts,1,skip,ext,0);
	}
	else
	{
	  // there is no other way: creating buffer to hold intermediate results
	  LTFAT_REAL *buffer = (LTFAT_REAL *) ltfat_malloc(outLen[0]*sizeof(LTFAT_REAL));
	  tmpOut[0] = buffer;
	  for(int ff=0;ff<noOfFilts-1;ff++)
	  {
		  tmpOut[ff+1] = out[noOfOutputs-(noOfFilts-1)+ff];
	  }
	  int skip = 0;//(fLen+1)/2;
	  LTFAT_NAME(conv_td_sub)(in,inLen,tmpOut,outLen[noOfOutputs-1],filtsNorm,fLen,noOfFilts,1,skip,ext,0);

	  tmpOut[0] = buffer;
	  for(int jj=1;jj<J-1;jj++)
	  {
	     for(int ff=0;ff<noOfFilts-1;ff++)
	     {
	        tmpOut[1+ff] = out[noOfOutputs-(jj+1)*(noOfFilts-1)+ff];
	     }
		 int outLenIdx = noOfOutputs-(jj+1)*(noOfFilts-1);
		 int inLenIdx = noOfOutputs-jj*(noOfFilts-1);
	     LTFAT_NAME(conv_td_sub)(buffer,outLen[inLenIdx],tmpOut,outLen[outLenIdx],filtsNorm,fLen,noOfFilts,1,skip,ext,ipow(subFac,jj));
	  }
	
	  skip = 0;
	  tmpOut[0] = out[0]; 
	  for(int ff=0;ff<noOfFilts-1;ff++)
	  {
	     tmpOut[1+ff] = out[1+ff];
	  }
	  LTFAT_NAME(conv_td_sub)(buffer,outLen[noOfFilts],tmpOut,outLen[0],filtsNorm,fLen,noOfFilts,1,skip,ext,ipow(subFac,J-1));

    
	 /****** STEP 3: FREE ALL ******/
	  ltfat_free(buffer);
	}
    for(int ff=0;ff<noOfFilts;ff++)  ltfat_free(filtsNorm[ff]); 
	ltfat_free(filtsNorm);
	ltfat_free(tmpOut);
}

LTFAT_EXTERN
	void LTFAT_NAME(dyadic_dwt_per)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *out[], const int outLen[], const LTFAT_REAL *filts[], int fLen, int noOfFilts, int subFac, const int J)
{
	LTFAT_REAL** tmpOut= (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	int skip = (fLen+1)/2;
	int extType = 7;
	if(J<=1)
	{
      for(int ff=0;ff<noOfFilts;ff++)
	  {
	     tmpOut[ff] = out[ff];
	  }
	  LTFAT_NAME(conv_td_sub)(in,inLen,tmpOut,outLen[0],filts,fLen,noOfFilts,subFac,skip,extType,0);
	}
	else
	{
	  // there is no other way: creating buffer to hold intermediate results
	  int noOfOutputs = J*(noOfFilts-1)+1;
	  LTFAT_REAL *buffer = (LTFAT_REAL *) ltfat_malloc(outLen[noOfOutputs-1]*sizeof(LTFAT_REAL));


	  for(int ff=0;ff<noOfFilts-1;ff++)
	  {
	     tmpOut[1+ff] = out[noOfOutputs-(noOfFilts-1)+ff];
	  }
	  tmpOut[0] = buffer;

	  // tmpOut[1] = out[J];
	  LTFAT_NAME(conv_td_sub)(in,inLen,tmpOut,outLen[noOfOutputs-1],filts,fLen,noOfFilts,subFac,skip,extType,0);

	  for(int jj=1;jj<J-1;jj++)
	  {
		 for(int ff=0;ff<noOfFilts-1;ff++)
	     {
	        tmpOut[1+ff] = out[noOfOutputs-(jj+1)*(noOfFilts-1)+ff];
	     }
		 int outLenIdx = noOfOutputs-(jj+1)*(noOfFilts-1);
		 int inLenIdx = noOfOutputs-jj*(noOfFilts-1);
	     // tmpOut[1] = out[J-jj]; 
	     //LTFAT_NAME(conv_td_sub)(buffer,outLen[J-jj+1],tmpOut,outLen[J-jj],filts,fLen,noOfFilts,subFac,skip,extType,0);
		 LTFAT_NAME(conv_td_sub)(buffer,outLen[inLenIdx],tmpOut,outLen[outLenIdx],filts,fLen,noOfFilts,subFac,skip,extType,0);
	  }
	
	  tmpOut[0] = out[0]; 
	  for(int ff=0;ff<noOfFilts-1;ff++)
	  {
	     tmpOut[1+ff] = out[1+ff];
	  }
	  //tmpOut[1] = out[1]; 
	  LTFAT_NAME(conv_td_sub)(buffer,outLen[noOfFilts],tmpOut,outLen[0],filts,fLen,noOfFilts,subFac,skip,extType,0);
    
	  ltfat_free(buffer);
	}

	ltfat_free(tmpOut);
}

LTFAT_EXTERN
	void LTFAT_NAME(dyadic_dwt_exp)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *out[], const int outLen[], const LTFAT_REAL *filts[], int fLen, int noOfFilts, int subFac, const int J, int ext)
{
	int noOfOutputs = J*(noOfFilts-1)+1;
	LTFAT_REAL** tmpOut= (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	int skip=1;
	
	if(J<=1)
	{
	  for(int ff=0;ff<noOfFilts;ff++)
	  {
	     tmpOut[ff] = out[ff];
	  }
	  LTFAT_NAME(conv_td_sub)(in,inLen,tmpOut,outLen[0],filts,fLen,noOfFilts,subFac,1,ext,0);
	}
	else
	{
		// determine maximum length of the coefficient vectors to be used as buffer
		int maxOutLen = outLen[noOfOutputs-1];
		for(int xx=0;xx<noOfOutputs;xx++)
		{
		  if(outLen[xx]>maxOutLen) maxOutLen = outLen[xx];
		}

	  LTFAT_REAL *buffer = (LTFAT_REAL *) ltfat_malloc(maxOutLen*sizeof(LTFAT_REAL));

	  for(int ff=0;ff<noOfFilts-1;ff++)
	  {
	     tmpOut[1+ff] = out[noOfOutputs-(noOfFilts-1)+ff];
	  }
	  tmpOut[0] = buffer;
	  LTFAT_NAME(conv_td_sub)(in,inLen,tmpOut,outLen[noOfOutputs-1],filts,fLen,noOfFilts,subFac,skip,ext,0);

	  //tmpOut[0] = buffer;
	  for(int jj=1;jj<J-1;jj++)
	  {
		  for(int ff=0;ff<noOfFilts-1;ff++)
	     {
	        tmpOut[1+ff] = out[noOfOutputs-(jj+1)*(noOfFilts-1)+ff];
	     }
		 int outLenIdx = noOfOutputs-(jj+1)*(noOfFilts-1);
		 int inLenIdx = noOfOutputs-jj*(noOfFilts-1);

	     LTFAT_NAME(conv_td_sub)(buffer,outLen[inLenIdx],tmpOut,outLen[outLenIdx],filts,fLen,noOfFilts,subFac,skip,ext,0);
	  }
	
	  tmpOut[0] = out[0]; 
	  for(int ff=0;ff<noOfFilts-1;ff++)
	  {
	     tmpOut[1+ff] = out[1+ff];
	  }
	  LTFAT_NAME(conv_td_sub)(buffer,outLen[noOfFilts],tmpOut,outLen[0],filts,fLen,noOfFilts,subFac,skip,ext,0);
    
	  ltfat_free(buffer);
	}

	ltfat_free(tmpOut);
}




LTFAT_EXTERN
void LTFAT_NAME(dyadic_idwt_exp)(LTFAT_REAL *in[], int inLen[], LTFAT_REAL *out, const int outLen, const LTFAT_REAL *filts[], int fLen, int noOfFilts, int up, const int J)
{
	int filtUps = 1;
	int ups = up;
	int noOfInputs = J*(noOfFilts-1)+1;
	//const LTFAT_REAL* tmpIn[2];
	LTFAT_REAL** tmpIn = (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	//int skip = fLen-1-(pow2(ups)-1);
	int skip = fLen-2;
	if(J<=1)
	{
	  tmpIn[0] = in[0];
	  for(int ff=0;ff<noOfFilts-1;ff++)
	  {
	     tmpIn[ff+1] = in[ff+1]; 
	  }
	  LTFAT_NAME(up_conv_td)(tmpIn,inLen[0],out,outLen,filts,fLen,noOfFilts,ups,skip,0,filtUps);
	}
	else
	{
		int maxLen = outLen;
		for(int xx=0;xx<noOfInputs;xx++)
		{
		  if(inLen[xx]>maxLen) maxLen = inLen[xx];
		}

		LTFAT_REAL *buffer2 = out;
		if(maxLen>outLen)
		{
			buffer2  = (LTFAT_REAL *) ltfat_malloc(maxLen*sizeof(LTFAT_REAL));
		}

		LTFAT_REAL *buffer = (LTFAT_REAL *) ltfat_malloc(maxLen*sizeof(LTFAT_REAL));
	    tmpIn[0] = in[0]; 
		for(int ff=0;ff<noOfFilts-1;ff++)
	    {
	       tmpIn[ff+1] = in[ff+1]; 
	    }

		if(modPow2(J,2)==0)
		{
			LTFAT_NAME(up_conv_td)(tmpIn,inLen[0],buffer,inLen[noOfFilts],filts, fLen,noOfFilts,ups,skip,0,filtUps);
		}
		else
		{
	       LTFAT_NAME(up_conv_td)(tmpIn,inLen[0],buffer2,inLen[noOfFilts],filts, fLen,noOfFilts,ups,skip,0,filtUps);
		}

		for(int jj=J-1;jj>1;jj--)
		{
		   // tmpIn[1] = in[(J-1)-jj+2]; 
		   for(int ff=0;ff<noOfFilts-1;ff++)
	       {
	           tmpIn[ff+1] = in[noOfInputs-jj*(noOfFilts-1)+ff]; 
	       }

		   //tmpIn[1] = in[(J-1)-jj+2]; 
		   int actInLenIdx = noOfInputs-jj*(noOfFilts-1);
		   int nextInLenIdx = noOfInputs-(jj-1)*(noOfFilts-1);

		   if(modPow2(jj,2)==0)
		   {
			  tmpIn[0] = buffer2;
	         // LTFAT_NAME(up_conv_td)(tmpIn,inLen[(J-1)-jj+2],buffer,inLen[(J-1)-jj+3],filts,fLen,noOfFilts,ups,skip,0,filtUps);
			   LTFAT_NAME(up_conv_td)(tmpIn,inLen[actInLenIdx],buffer,inLen[nextInLenIdx],filts,fLen,noOfFilts,ups,skip,0,filtUps);
		   }
		   else
		   {
			  tmpIn[0] = buffer;
	          //LTFAT_NAME(up_conv_td)(tmpIn,inLen[(J-1)-jj+2],buffer2,inLen[(J-1)-jj+3],filts,fLen,noOfFilts,ups,skip,0,filtUps);
			  LTFAT_NAME(up_conv_td)(tmpIn,inLen[actInLenIdx],buffer2,inLen[nextInLenIdx],filts,fLen,noOfFilts,ups,skip,0,filtUps);
		   }
		}

		for(int ff=0;ff<noOfFilts-1;ff++)
	    {
			tmpIn[ff+1] = in[noOfInputs-(noOfFilts-1)+ff]; 
	    }
		tmpIn[0] = buffer;
	    LTFAT_NAME(up_conv_td)(tmpIn,inLen[noOfInputs-1],out,outLen,filts,fLen,noOfFilts,ups,skip,0,filtUps);


		ltfat_free(buffer);
		if(buffer2!=out) ltfat_free(buffer2);
	}
	ltfat_free(tmpIn);
}

LTFAT_EXTERN
void LTFAT_NAME(dyadic_idwt_per)(LTFAT_REAL *in[], int inLen[], LTFAT_REAL *out, const int outLen, const LTFAT_REAL *filts[], int fLen, int noOfFilts, int up, const int J)
{
	int filtUps = 1;
	int ups = up;
	int noOfInputs = J*(noOfFilts-1)+1;

	LTFAT_REAL** tmpIn = (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	int skip = (fLen+1)/2 -1;

	if(J<=1)
	{
	  tmpIn[0] = in[0];
	  for(int ff=0;ff<noOfFilts-1;ff++)
	  {
	     tmpIn[ff+1] = in[ff+1]; 
	  }
	  LTFAT_NAME(up_conv_td)(tmpIn,inLen[0],out,outLen,filts,fLen,noOfFilts,ups,skip,1,filtUps);
	}
	else
	{
		LTFAT_REAL *buffer = (LTFAT_REAL *) ltfat_malloc(outLen*sizeof(LTFAT_REAL));
	    tmpIn[0] = in[0]; 
		for(int ff=0;ff<noOfFilts-1;ff++)
	    {
	       tmpIn[ff+1] = in[ff+1]; 
	    }

		if(modPow2(J,2)==0)
		{
			LTFAT_NAME(up_conv_td)(tmpIn,inLen[0],buffer,inLen[noOfFilts],filts, fLen,noOfFilts,ups,skip,1,filtUps);
		}
		else
		{
	       LTFAT_NAME(up_conv_td)(tmpIn,inLen[0],out,inLen[noOfFilts],filts, fLen,noOfFilts,ups,skip,1,filtUps);
		}

		for(int jj=J-1;jj>1;jj--)
		{
           for(int ff=0;ff<noOfFilts-1;ff++)
	       {
	           tmpIn[ff+1] = in[noOfInputs-jj*(noOfFilts-1)+ff]; 
	       }

		   //tmpIn[1] = in[(J-1)-jj+2]; 
		   int actInLenIdx = noOfInputs-jj*(noOfFilts-1);
		   int nextInLenIdx = noOfInputs-(jj-1)*(noOfFilts-1);
		   if(modPow2(jj,2)==0)
		   {
			  tmpIn[0] = out;
	          LTFAT_NAME(up_conv_td)(tmpIn,inLen[actInLenIdx],buffer,inLen[nextInLenIdx],filts,fLen,noOfFilts,ups,skip,1,filtUps);
		   }
		   else
		   {
			  tmpIn[0] = buffer;
	          LTFAT_NAME(up_conv_td)(tmpIn,inLen[actInLenIdx],out,inLen[nextInLenIdx],filts,fLen,noOfFilts,ups,skip,1,filtUps);
		   }
		}

		//tmpIn[1] = in[J]; 
		for(int ff=0;ff<noOfFilts-1;ff++)
	    {
			tmpIn[ff+1] = in[noOfInputs-(noOfFilts-1)+ff]; 
	    }
		tmpIn[0] = buffer;
	    LTFAT_NAME(up_conv_td)(tmpIn,inLen[noOfInputs-1],out,outLen,filts,fLen,noOfFilts,ups,skip,1,filtUps);


		ltfat_free(buffer);
	}
	ltfat_free(tmpIn);

}

LTFAT_EXTERN
void LTFAT_NAME(undec_idwt_exp)(LTFAT_REAL *in[], int inLen[], LTFAT_REAL *out, const int outLen, const LTFAT_REAL *filts[], int fLen, int noOfFilts, int up, const int J)
{
	LTFAT_REAL **filtsNorm = (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	LTFAT_REAL **tmpIn = (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	for(int ff=0;ff<noOfFilts;ff++)
	{
	   filtsNorm[ff] = (LTFAT_REAL*) ltfat_malloc(fLen*sizeof(LTFAT_REAL));
	   for(int fid=0;fid<fLen;fid++)
	   {
		   filtsNorm[ff][fid] = (LTFAT_REAL) filts[ff][fid]/sqrt((LTFAT_REAL)up);
	   }
	}

	int filtUps = 1;
	int ups = 1;
	int tmpFlen = fLen;
	int skip = fLen-1;
	int noOfInputs = J*(noOfFilts-1)+1;


	if(J<=1)
	{
	  for(int ff=0;ff<noOfFilts;ff++)
	  {
		  tmpIn[ff] = in[ff];
	  }
	  LTFAT_NAME(up_conv_td)(tmpIn,inLen[0],out,outLen,filtsNorm,fLen,noOfFilts,ups,skip,0,filtUps);
	}
	else
	{
		filtUps = J-1;
		tmpFlen = ipow(up,filtUps)*fLen-(ipow(up,filtUps)-1);
		skip = tmpFlen-1;
		int buffLen = inLen[2];
		LTFAT_REAL *buffer = (LTFAT_REAL *) ltfat_malloc(buffLen*sizeof(LTFAT_REAL));
	    tmpIn[0] = in[0]; 
	    for(int ff=0;ff<noOfFilts;ff++)
	    {
		  tmpIn[ff] = in[ff];
	    }
        LTFAT_NAME(up_conv_td)(tmpIn,inLen[0],buffer,inLen[noOfFilts],filtsNorm,fLen,noOfFilts,ups,skip,0,ipow(up,filtUps));


		tmpIn[0] = buffer;
		for(int jj=J-1;jj>1;jj--)
		{
			filtUps = jj-1;
			tmpFlen = ipow(up,filtUps)*fLen-(ipow(up,filtUps)-1);
			skip = tmpFlen-1;
			//int currjj=(J-1)-jj+2;
		    //tmpIn[1] = in[currjj]; 
			for(int ff=0;ff<noOfFilts-1;ff++)
	       {
	           tmpIn[ff+1] = in[noOfInputs-jj*(noOfFilts-1)+ff]; 
	       }

		   //tmpIn[1] = in[(J-1)-jj+2]; 
		   int actInLenIdx = noOfInputs-jj*(noOfFilts-1);
		   int nextInLenIdx = noOfInputs-(jj-1)*(noOfFilts-1);
			  
	          LTFAT_NAME(up_conv_td)(tmpIn,inLen[actInLenIdx],buffer,inLen[nextInLenIdx],filtsNorm,fLen,noOfFilts,ups,skip,0,ipow(up,filtUps));

		}

		filtUps = 0;
		tmpFlen = fLen;
		skip = tmpFlen-1;
		for(int ff=0;ff<noOfFilts-1;ff++)
	    {
			tmpIn[ff+1] = in[noOfInputs-(noOfFilts-1)+ff]; 
	    }

	    LTFAT_NAME(up_conv_td)(tmpIn,inLen[noOfInputs-1],out,outLen,filtsNorm,fLen,noOfFilts,ups,skip,0,filtUps);



		ltfat_free(buffer);
	}
	for(int ff=0;ff<noOfFilts;ff++)  ltfat_free(filtsNorm[ff]); 
	ltfat_free(filtsNorm);
	ltfat_free(tmpIn);
}

LTFAT_EXTERN
void LTFAT_NAME(undec_idwt_per)(LTFAT_REAL *in[], int inLen, LTFAT_REAL *out, const int outLen, const LTFAT_REAL *filts[], int fLen, int noOfFilts, int up, const int J)
{
	int filtUps = 1;
	int ups = 1;
	int tmpFlen = fLen;
	int skip = (tmpFlen)/2 - 1 ;
	int upLog2 = log2(up);

	LTFAT_REAL **filtsNorm = (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	LTFAT_REAL **tmpIn = (LTFAT_REAL**) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));
	for(int ff=0;ff<noOfFilts;ff++)
	{
	   filtsNorm[ff] = (LTFAT_REAL*) ltfat_malloc(fLen*sizeof(LTFAT_REAL));
	   for(int fid=0;fid<fLen;fid++)
	   {
		   filtsNorm[ff][fid] = (LTFAT_REAL) filts[ff][fid]/sqrt((LTFAT_REAL)up);
	   }
	}
	int noOfInputs = J*(noOfFilts-1) + 1;

	if(J<=1)
	{
	  for(int ff=0;ff<noOfFilts;ff++)
	  {
		  tmpIn[ff] = in[ff];
	  }
	  LTFAT_NAME(up_conv_td)(tmpIn,inLen,out,outLen,filtsNorm,fLen,noOfFilts,ups,skip,1,filtUps);
	}
	else
	{
		filtUps = J-1;
		tmpFlen = ipow(up,filtUps)*fLen-(ipow(up,filtUps)-1);
		skip = (ipow(up,filtUps)*fLen)/2-ipow(up,filtUps);
		int buffLen = inLen;
		LTFAT_REAL *buffer = (LTFAT_REAL *) ltfat_malloc(buffLen*sizeof(LTFAT_REAL));
	    //tmpIn[0] = in[0]; tmpIn[1] = in[1]; 

	    for(int ff=0;ff<noOfFilts;ff++)
	    {
			tmpIn[ff] = in[ff];
	    }

	    LTFAT_NAME(up_conv_td)(tmpIn,inLen,buffer,inLen,filtsNorm,fLen,noOfFilts,ups,skip,1,ipow(up,filtUps));

		tmpIn[0] = buffer;
		for(int jj=J-1;jj>1;jj--)
		{
			filtUps = jj-1;
			tmpFlen = ipow(up,filtUps)*fLen-(ipow(up,filtUps)-1);
		    skip = (ipow(up,filtUps)*fLen)/2-ipow(up,filtUps);
			int currjj=(J-1)-jj+2;
		    
			for(int ff=0;ff<noOfFilts-1;ff++)
	        {
	           tmpIn[ff+1] = in[noOfInputs-jj*(noOfFilts-1)+ff]; 
	        }
			//tmpIn[1] = in[currjj]; 
			  
	        LTFAT_NAME(up_conv_td)(tmpIn,inLen,buffer,inLen,filtsNorm,fLen,noOfFilts,ups,skip,1,ipow(up,filtUps));
		}

		filtUps = 1;
		tmpFlen = fLen;
		skip = (tmpFlen)/2 -1;
		//tmpIn[1] = in[J]; 
		for(int ff=0;ff<noOfFilts-1;ff++)
	    {
			tmpIn[ff+1] = in[noOfInputs-(noOfFilts-1)+ff]; 
	    }
	    LTFAT_NAME(up_conv_td)(tmpIn,inLen,out,outLen,filtsNorm,fLen,noOfFilts,ups,skip,1,filtUps);


		/***** FREE STUF *******/
		ltfat_free(buffer);
	}
        for(int ff=0;ff<noOfFilts;ff++)  ltfat_free(filtsNorm[ff]); 
	    ltfat_free(filtsNorm);
	    ltfat_free(tmpIn);
}



/*
BASIC routine for the disrete wavelet transforms
noOfFilts - filters filterbank followed by subsampling by factor of sub
skip - initial possition of the filters
ext - extension type
filtUp - filters upsampling 


all filters of equal length and identical subsampling for each of them

if filtUps>=1, the aTrous algorithm is used (as if the input filter was upsampled 2^filtUps, but just nonzero samples are used in computations)
else  usual convolution combined with subsampling is used

REMARK: Output subsampling works even for subFilt>1
REMARK 2: if sub >= filtLen, no buffer cycling is necessary...buffer is filled with all new samples every time. Cycling works, but it is slower than it could be
*/

LTFAT_EXTERN
void LTFAT_NAME(conv_td_sub)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *out[], const int outLen, const LTFAT_REAL *filts[], int fLen, int noOfFilts, int sub, int skip, int ext, int filtUps)
{
	if(filtUps<1) filtUps=1;
	int filtUpsPow2 = filtUps;
	filtUps = log2(filtUpsPow2);
	/*if(filtUps<0) 
	{
		filtUps=0;
	}
	else
	{
	  filtUpsPow2 = pow2(filtUps);
	}*/
	int fLenUps = filtUpsPow2*fLen-(filtUpsPow2-1);

	LTFAT_REAL *righExtbuff = 0;
	// number of output samples that can be calculated "painlessly"
    int outLenN = (inLen - skip + sub -1)/sub;
	
   // prepare cyclic buffer of length of power of two (for effective modulo operations)
   int buffLen = nextPow2(imax(fLenUps,sub+1));
   // buffer index
   int buffPtr = 0;
   // pointer for moving in the input data
   const LTFAT_REAL *tmpIn = in;
   // allocating and initializing the cyclic buffer
   LTFAT_REAL *buffer = (LTFAT_REAL *) ltfat_malloc(buffLen*sizeof(LTFAT_REAL));
   memset(buffer,0,buffLen*sizeof(LTFAT_REAL)); 
   
   // fill buffer with the initial values from the input signal according to the boundary treatment
   // last fLenUps buffer samples are filled to keep buffPtr=0
   extend_left(in,inLen,buffer,buffLen,fLenUps,ext,sub);

   if(outLenN<outLen)
   {
   	   // right extension is necessary, additional buffer from where to copy
	   righExtbuff = (LTFAT_REAL *) ltfat_malloc(buffLen*sizeof(LTFAT_REAL)); 
       memset(righExtbuff,0,buffLen*sizeof(LTFAT_REAL)); 
	   // store extension in the buffer (must be done now to avoid errors when inplace calculation is done)
	   extend_right(in,inLen,righExtbuff,fLenUps,ext,sub);
   }


   /*** initial buffer fill ***/ 
   // number of overflowing samples
   int buffOver = imax(buffPtr+(skip+1)-buffLen, 0);
   // copy non-overflowing
   memcpy(buffer+buffPtr,tmpIn,(skip+1-buffOver)*sizeof(LTFAT_REAL));
   tmpIn += (skip+1-buffOver);
   // copy overflowing
   memcpy(buffer,tmpIn,buffOver*sizeof(LTFAT_REAL));
   tmpIn += buffOver;
   // move buffer pointer (can overflow)
   buffPtr = modPow2(buffPtr+=(skip+1),buffLen);

   // accumulators for each output channel
   LTFAT_REAL *outTmp = (LTFAT_REAL *) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL));
   memset(outTmp,0,noOfFilts*sizeof(LTFAT_REAL));
   // pointer to the actual buffer
   LTFAT_REAL *buffTmp; //= (LTFAT_REAL) 0.0;


   /*****************************************************************************/
   int packDtNo = (16/sizeof(LTFAT_REAL));
   int inBuffLen = packDtNo*((fLen+packDtNo-1)/packDtNo);
   LTFAT_REAL *inBuff = (LTFAT_REAL *) ltfat_malloc(inBuffLen*sizeof(LTFAT_REAL));
   /*****************************************************************************/
   
   /*********** STEP 1: FREE LUNCH ( but also a hot-spot) *******************************/
   // Take the smaller value from "painless" output length and the user defined output length
   int iiLoops = imin(outLenN-1,outLen-1);
   // ceil(fLen/subFilt) number of samples of impluse responses actually used
   int jjLoops = fLen; // (fLen+filtUpsPow2-1)/filtUpsPow2;

   // loop trough all output samples, omit the very last one.
   for (int ii = 0; ii < iiLoops; ii++) 
   {
	   
	   // loop trough impulse response samples
	   for (int jj = 0; jj < jjLoops; jj++) 
	   {
          // buffer index modulo buffer length
		  // int idx = modPow2((-(jj<<filtUps)+buffPtr-1),buffLen);
		   int idx = modPow2((-(jj*filtUpsPow2)+buffPtr-1),buffLen);
		
		  // pinter to the actual value
          buffTmp = buffer+idx;
		  // loop trough filters
		  for(int kk = 0;kk<noOfFilts;kk++)
		  {
				outTmp[kk] += *buffTmp * filts[kk][jj];
		  }
       }

	   // write accumulated values to the output
	   for(int kk = 0;kk<noOfFilts;kk++)
	   {
		 out[kk][ii] = outTmp[kk];
	   }

	   /*
	    for (int jj = 0; jj < jjLoops; jj++) 
	   {
          // buffer index modulo filter buffer length
		  int idx = modPow2((-(jj<<filtUps)+buffPtr-1),buffLen);
		  // actual value
          *(inBuff+jj) = *(buffer+idx);
       }

		for(int kk = 0;kk<noOfFilts;kk++)
	   {
		 LTFAT_REAL* outPtr = &out[kk][ii]; 
		 const LTFAT_REAL* filtPtr = &filts[kk][0]; 
		 *outPtr = (LTFAT_REAL) 0.0;
		 for(int ll=0;ll<fLen;ll++)
		 {
			 *outPtr += inBuff[ll]*filtPtr[ll];
		 }
	   }
	   */

	   // clear temp variables
	   memset(outTmp,0,noOfFilts*sizeof(LTFAT_REAL));

	   // number of samples overflowing buffer
	   int buffOver = imax(buffPtr+sub-buffLen, 0);
	   // have to copy just non-overfloving samples
   	   memcpy(buffer + buffPtr, tmpIn, (sub-buffOver)*sizeof(LTFAT_REAL));
	   // copy the rest to the beginning of the buffer
	   memcpy(buffer,tmpIn+sub-buffOver,buffOver*sizeof(LTFAT_REAL));

	   // move buffer index, modulo by buffer length
	   buffPtr = modPow2(buffPtr += sub,buffLen);
	   // move input pointer
	   tmpIn += sub;
   }

  /*********** STEP 2: FINALIZE FREE LUNCH ************************************/
    // calculate last sample of the full response
	   for (int jj = 0; jj < jjLoops; jj++) 
	   {
          // buffer index modulo filter buffer length
		  //int idx = modPow2((-(jj<<filtUps)+buffPtr-1),buffLen);
		   int idx = modPow2((-(jj*filtUpsPow2)+buffPtr-1),buffLen);
		  // actual value
          buffTmp = buffer+idx;
		  // loop trough filters
		  for(int kk = 0;kk<noOfFilts;kk++)
		  {
				outTmp[kk] += *buffTmp * filts[kk][jj];
		  }
       }

	   // write accumulated values to the output
	   for(int kk = 0;kk<noOfFilts;kk++)
	   {
		   out[kk][outLenN-1] = outTmp[kk];
	   }

  /*********** STEP 3: NOW FOR THE TRICKY PART ************************************/
  if(outLenN<outLen)
  {
	  /************ STEP 3a: DEALING WITH THE REMAINING SAMPLES ******************/
	  // CAREFULL NOW! possibly stepping outside of input signal
	  // last index in the input signal for which reading next sub samples reaches outside of the input signal
	  int lastInIdx = (sub*(outLenN-1)+1+skip);
	  // remaining samples: diff<sub
	  int diff = imax(0,(inLen) - lastInIdx);
	  // fill buffer with the remaining values, values outside of the input signal are read from different buffer
	   // number of samples overflowing buffer
	   int buffOver = imax(buffPtr+diff-buffLen, 0);
   	   memcpy(buffer + buffPtr, in + lastInIdx, (diff-buffOver)*sizeof(LTFAT_REAL));
	   memcpy(buffer,in + lastInIdx+diff-buffOver,buffOver*sizeof(LTFAT_REAL));
	   // move buffer index, modulo by buffer length
	   buffPtr = modPow2(buffPtr += diff,buffLen);


	   // now copying samples that are outside
	   buffOver = imax(buffPtr+(sub-diff)-buffLen, 0);
	   memcpy(buffer + buffPtr,righExtbuff, (sub-diff-buffOver)*sizeof(LTFAT_REAL));
	   memcpy(buffer,righExtbuff+(sub-diff)-buffOver,buffOver*sizeof(LTFAT_REAL));


	   buffPtr = modPow2(buffPtr += (sub-diff),buffLen);
	   // index to the right ext. buffer
	   int rightExtBuffIdx = sub-diff;

	/************ STEP 3b: ALL OK, proceed reading input values from righExtbuff ******************/
	// loop for the remaining output samples
	for(int ii=0;ii<outLen-outLenN;ii++)
	{
		 // clear temp variables
	   memset(outTmp,0,noOfFilts*sizeof(LTFAT_REAL));

	   for (int jj = 0; jj < jjLoops; jj++) 
	   {
          // buffer index modulo filter buffer length
		  // int idx = modPow2((-(jj<<filtUps)+buffPtr-1),buffLen);
		   int idx = modPow2((-(jj*filtUpsPow2)+buffPtr-1),buffLen);
		  // actual value
          buffTmp = buffer+idx;
		  // loop trough filters
		  for(int kk = 0;kk<noOfFilts;kk++)
		  {
				outTmp[kk] += *buffTmp * filts[kk][jj];
		  }
       }

	   // write accumulated values to the output
	   for(int kk = 0;kk<noOfFilts;kk++)
	   {
		   out[kk][outLenN+ii] = outTmp[kk];
	   }
	
	   int buffOver = imax(buffPtr+sub-buffLen, 0);
	   memcpy(buffer + buffPtr,righExtbuff+rightExtBuffIdx, (sub-buffOver)*sizeof(LTFAT_REAL));
	   memcpy(buffer,righExtbuff+rightExtBuffIdx+(sub-buffOver),buffOver*sizeof(LTFAT_REAL));

	   buffPtr = modPow2(buffPtr += sub,buffLen);
	   rightExtBuffIdx = modPow2(rightExtBuffIdx += sub,buffLen);
	
	}
	
  }

  if(righExtbuff)ltfat_free(righExtbuff);
   ltfat_free(outTmp);
   ltfat_free(buffer);
   ltfat_free(inBuff);
}

/*
SKIP in odd upsampled input
*/
LTFAT_EXTERN
void LTFAT_NAME(up_conv_td)(LTFAT_REAL *in[], int inLen, LTFAT_REAL *out, const int outLen, const LTFAT_REAL *filts[], int fLen, int noOfFilts, int up, int skip, int ext, int filtUps)
{
   if(filtUps<1) filtUps=1;
   int filtUpsPow2 = filtUps; // filtUps and up are kept separetely for purposes of time-invariant DWT
   filtUps = log2(filtUpsPow2);
   int upPow2 = up;
   int tmpfLen = filtUpsPow2*fLen - (filtUpsPow2-1);
//   int outLenN = inLen*up - skip;  // number of output samples, that can be calculated from input data given skip
   int buffLen = nextPow2(imax(tmpfLen,upPow2+1));
   int buffPtr = 0;
   int inStart = (skip + upPow2 -1)/upPow2;  // number of skipped input samples given skip (in output)
   LTFAT_REAL * tmpOut = out;
   LTFAT_REAL ** buffers = (LTFAT_REAL **) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*)); // as many bufers as input signals
   int iniStoCopy = imin(inStart,buffLen); // how many input samples to pre-load to the buffers

   int skipModUp = skip%upPow2;
   int skipToNextUp = 0;
   if(skipModUp!=0)  skipToNextUp = upPow2-skipModUp;
   int iiLoops = imin(inLen - inStart,(outLen-skipToNextUp)/upPow2); // just in case outLen/up < inLen - inStart
   int remainsOutSamp = outLen - ((iiLoops)*upPow2 + skipToNextUp);
   int remainsInSamp = (inLen - inStart) - (iiLoops);


   // int rightBuffLen = imax((remainsOutSamp>>up)-remainsInSamp,buffLen);
   int rightBuffLen = imax((remainsOutSamp/upPow2)-remainsInSamp,buffLen);
   LTFAT_REAL ** rightBuffers = (LTFAT_REAL **) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL*));

   for(int kk = 0;kk<noOfFilts;kk++)
   {
      buffers[kk] = (LTFAT_REAL *) ltfat_malloc(buffLen*sizeof(LTFAT_REAL));
	  rightBuffers[kk] = (LTFAT_REAL *) ltfat_malloc(rightBuffLen*sizeof(LTFAT_REAL));
      memset(buffers[kk],0,buffLen*sizeof(LTFAT_REAL));
	  memset(rightBuffers[kk],0,rightBuffLen*sizeof(LTFAT_REAL));
	  if(ext) // if periodic extension
      {
         extend_left(in[kk],inLen,buffers[kk],buffLen,tmpfLen,5,upPow2); // extension as a last (tmpfLen-1) samples of the buffer -> pointer dont have to be moved
		 extend_right(in[kk],inLen,rightBuffers[kk],tmpfLen,5,upPow2);
      }
	  memcpy(buffers[kk],in[kk],iniStoCopy*sizeof(LTFAT_REAL)); // TO DO: add some shift in in if iniStoCopy>buffLen
   }
   buffPtr = modPow2(buffPtr += iniStoCopy,buffLen);

   /*
   STEP 1: Deal with the shift - upsampling misaligment
   if skip=0,up,2*up ... aligned -> calculate nothing
   else calculate missing values to next aligment
   */


   for(int ii=0;ii<skipToNextUp;ii++)
   {
		*tmpOut =  (LTFAT_REAL) 0.0;  
		// int jjLoops = ((fLen-(skipModUp-ii)+upPow2-1)>>up); 
		int jjLoops = ((fLen-(skipModUp-ii)+upPow2-1)/upPow2); 

		    for(int jj=0;jj<jjLoops;jj++)
		    {
        		// int idx = modPow2((-(jj<<filtUps)+buffPtr-1), buffLen);
				int idx = modPow2((-(jj*filtUpsPow2)+buffPtr-1), buffLen);
     			for(int kk = 0;kk<noOfFilts;kk++)
		        {
					// *tmpOut += buffers[kk][idx] * filts[kk][(jj<<up)+skipModUp+ii];
					*tmpOut += buffers[kk][idx] * filts[kk][(jj*upPow2)+skipModUp+ii];
		        }
		    }
	    tmpOut++;
   }
   /******************************************************************************/

   /*
   STEP 2: MAIN LOOP. Number of iterations is equal to the number of input samples (minus already loaded ones) or is limited by the number of output samples
   */

	 // loop over input samples
	 for(int ii=0;ii<iiLoops;ii++)
	 {
		 for(int kk = 0;kk<noOfFilts;kk++)
		 {
			 *(buffers[kk] + buffPtr) = *(in[kk] + inStart + ii);
		 }
		 buffPtr = modPow2(++buffPtr,buffLen);

		 // in upsamled input, go over up-skipModUp values
		 for(int uu=0;uu<upPow2;uu++)
	     {
			 LTFAT_REAL tmpOutVal = (LTFAT_REAL) 0.0;
			//tmpOut[uu] = (LTFAT_REAL) 0.0; // clear output
		   // int jjLoops = ((fLen-uu+upPow2-1)>>up);  // jjLoopsMax or jjLoopsMax
			 int jjLoops = ((fLen-uu+upPow2-1)/upPow2);  // jjLoopsMax or jjLoopsMax

		    for(int jj=0;jj<jjLoops;jj++)
		    {
				// int idx = modPow2((-(jj<<filtUps)+buffPtr-1),buffLen);
        		int idx = modPow2((-(jj*filtUpsPow2)+buffPtr-1),buffLen);
     			for(int kk = 0;kk<noOfFilts;kk++)
		        {
					// tmpOutVal += buffers[kk][idx] * filts[kk][(jj<<up)+uu];
					tmpOutVal += buffers[kk][idx] * filts[kk][(jj*upPow2)+uu];
		        }
		    }
			tmpOut[uu] = tmpOutVal;
		 }

		 tmpOut+=upPow2;
	 }
   /*
    END OF STEP 2
    */
	
	 /*
	 STEP 3a: load additional input sample!s! (if there are any) or load zero
	 */
	 int rightBuffPtr = 0;

        if(remainsInSamp>0)
		  {
		     for(int kk = 0;kk<noOfFilts;kk++)
		     {
			    *(buffers[kk] + buffPtr) = *(in[kk] + inStart + iiLoops);
		     }
		     buffPtr = modPow2(++buffPtr,buffLen);
			 remainsInSamp--;
		  }
		  else
		  {
		 	 for(int kk = 0;kk<noOfFilts;kk++)
		     {
				 *(buffers[kk] + buffPtr) = (LTFAT_REAL) rightBuffers[kk][0];
		     }
		     buffPtr = modPow2(++buffPtr,buffLen);
			 rightBuffPtr++;
		  }

     /*
	 STEP 3b: calculate remaining output samples, supplying zeros as new values
	 */
       for(int ii=0;ii<remainsOutSamp;ii++)
       {
		  if(ii!=0 && ii%upPow2==0)
		  {
		    // just supply zeros to the buffers
		    for(int kk = 0;kk<noOfFilts;kk++)
		    {
			   *(buffers[kk] + buffPtr) = (LTFAT_REAL) rightBuffers[kk][rightBuffPtr];
		    }
			rightBuffPtr++;
		    buffPtr = modPow2(++buffPtr,buffLen);
		  }

		   *tmpOut =  (LTFAT_REAL) 0.0;

		   //  int jjLoops = ((fLen-modPow2(ii,upPow2)+upPow2-1)>>up);
		   int jjLoops = ((fLen-ii%upPow2+upPow2-1)/upPow2);

		    for(int jj=0;jj<jjLoops;jj++)
		    {
        		// int idx = modPow2((-(jj<<filtUps)+buffPtr-1),buffLen);
				int idx = modPow2((-(jj*filtUpsPow2)+buffPtr-1),buffLen);
     			for(int kk = 0;kk<noOfFilts;kk++)
		        {
					//*tmpOut += buffers[kk][idx] * filts[kk][(jj<<up) + modPow2(ii,upPow2) ];
					*tmpOut += buffers[kk][idx] * filts[kk][(jj*upPow2) + ii%upPow2 ];
		        }
		    }
		   tmpOut++;
	    }


   for(int kk = 0;kk<noOfFilts;kk++)
   {
	   ltfat_free(buffers[kk]);
	   ltfat_free(rightBuffers[kk]);
   }
   ltfat_free(buffers);
   ltfat_free(rightBuffers);

}



LTFAT_EXTERN
void LTFAT_NAME(up_conv_sub)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *out, int outLen, LTFAT_REAL *filt, int fLen, int up, int sub, int skip, int ext){
  	skip = 0;
    // type-2 polyphase representation
	LTFAT_REAL **pphase = (LTFAT_REAL **) ltfat_malloc(up*sizeof(LTFAT_REAL*));
	int pphaseLen = (fLen+up-1)/up;
	for(int ii=0;ii<up;ii++)
	{
		LTFAT_REAL **tmpPPhase = &pphase[up-1-ii];
		*tmpPPhase = (LTFAT_REAL *) ltfat_malloc(pphaseLen*sizeof(LTFAT_REAL));
		memset(*tmpPPhase,0,pphaseLen*sizeof(LTFAT_REAL));
		for(int jj=0;jj<(fLen-ii+up-1)/up;jj++)
		{
			(*tmpPPhase)[jj] = filt[ii+up*jj]; 
		}
	}


   

	// prepare cyclic buffer of length of power of two (for effective modulo operations)
   int buffLen = nextPow2(imax(fLen,sub+1));
   // buffer index
   int buffPtr = 0;
   // pointer for moving in the input data
   const LTFAT_REAL *tmpIn = in;
   // allocating and initializing the cyclic buffer
   LTFAT_REAL *buffer = (LTFAT_REAL *) ltfat_malloc(buffLen*sizeof(LTFAT_REAL));
   memset(buffer,0,buffLen*sizeof(LTFAT_REAL));

   /*** initial buffer fill ***/ 
   int toRead = skip + (up-1)+ 1;
   // number of overflowing samples
   int buffOver = imax(buffPtr+toRead-buffLen, 0);
   // copy non-overflowing
   memcpy(buffer+buffPtr,tmpIn,(toRead-buffOver)*sizeof(LTFAT_REAL));
   tmpIn += (toRead-buffOver);
   // copy overflowing
   memcpy(buffer,tmpIn,buffOver*sizeof(LTFAT_REAL));
   tmpIn += buffOver;
   // move buffer pointer (can overflow)
   buffPtr = modPow2(buffPtr+=toRead,buffLen);

     
   int outLenN = (inLen - (toRead-1) + sub -1)/sub;
   // Take the smaller value from "painless" output length and the user defined output length
   int iiLoops = imin(outLenN-1,outLen-1);
   int jjLoops = pphaseLen;
   LTFAT_REAL *outTmp = out;
   LTFAT_REAL *buffTmp = buffer;

   // loop trough all output samples, omit the very last one.
   for (int ii = 0; ii < iiLoops; ii++) 
   {
	   // polyphase repr. idx
	   for (int kk = 0; kk < up; kk++) 
	   {
		   LTFAT_REAL* outkkTmp = (outTmp-(up-1-kk));
	      // polyphase repr. length
	      for (int jj = 0; jj < jjLoops; jj++) 
	      {
		     // buffer index modulo buffer length
		     int idx = modPow2((-kk-jj+buffPtr-1),buffLen);
		     // ponter to the actual value
             buffTmp = buffer+idx;
		     // loop trough polyphase samples
		     *outkkTmp += *(buffTmp)*pphase[kk][jj];
          }
	   }
	   outTmp+=up;


	   // number of samples overflowing buffer
	   int buffOver = imax(buffPtr+sub-buffLen, 0);
	   // have to copy just non-overfloving samples
   	   memcpy(buffer + buffPtr, tmpIn, (sub-buffOver)*sizeof(LTFAT_REAL));
	   // copy the rest to the beginning of the buffer
	   memcpy(buffer,tmpIn+sub-buffOver,buffOver*sizeof(LTFAT_REAL));

	   // move buffer index, modulo by buffer length
	   buffPtr = modPow2(buffPtr += sub,buffLen);
	   // move input pointer
	   tmpIn += sub;
   }








   for(int ii=0;ii<up;ii++) 
	   ltfat_free(pphase[ii]);

   ltfat_free(pphase);
   ltfat_free(buffer);
}


LTFAT_EXTERN
void LTFAT_NAME(up_conv_sub_1toN)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *out[], const int outLen, const LTFAT_REAL *filts[], int fLen, int noOfFilts, int sub, int skip, int ext)
{


}


// fills last buffer samples
LTFAT_EXTERN
void LTFAT_NAME(extend_left)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *buffer,int buffLen, int filtLen, int type, int a){
		int legalExtLen = imin(filtLen-1, inLen);
		LTFAT_REAL *buffTmp = buffer + buffLen - legalExtLen;
	switch (type) {
		case 1: // half-point symmetry
           for(int ii=0;ii<legalExtLen;ii++)
			   buffTmp[ii] = in[legalExtLen-ii-1];
           break;
		case 2: // whole-point symmetry
			legalExtLen = imin(filtLen-1, inLen-1);
			buffTmp = buffer + buffLen - legalExtLen;
           for(int ii=0;ii<legalExtLen;ii++)
			   buffTmp[ii] = in[legalExtLen-ii];
           break;
		 case 3: // half-point antisymmetry
           for(int ii=0;ii<legalExtLen;ii++)
			   buffTmp[ii] = -in[legalExtLen-ii-1];
           break;
		 case 4: // whole-point antisymmetry
			legalExtLen = imin(filtLen-1, inLen-1);
			legalExtLen = imin(filtLen-1, inLen-1);
			buffTmp = buffer + buffLen - legalExtLen;
           for(int ii=0;ii<legalExtLen;ii++)
			   buffTmp[ii] = -in[legalExtLen-ii];
           break;
		 case 5: // periodic padding
           for(int ii=0;ii<legalExtLen;ii++)
			   buffTmp[ii] = in[inLen-1-(legalExtLen-1)+ii];
           break;
		 case 6: // constant padding
			buffTmp = buffer + buffLen - (filtLen-1);
           for(int ii=0;ii<filtLen-1;ii++)
			   buffTmp[ii] = in[0];
           break;
		   case 7: // periodic padding with possible last sample repplication
            int rem = inLen%a;
			if(rem==0)
			{
              for(int ii=0;ii<legalExtLen;ii++)
			     buffTmp[ii] = in[inLen-1-(legalExtLen-1)+ii];
			}
			else
			{
			  int remto = a - rem;

			  // replicated
			  for(int ii=0;ii<remto;ii++)
			     buffTmp[legalExtLen-1-ii] = in[inLen-1];
			  
			  // periodic extension
			  for(int ii=0;ii<legalExtLen-remto;ii++)
			     buffTmp[ii] = in[inLen-1-(legalExtLen-1-1)+ii+remto-1];
			}
           break;
		case 0: // zero-padding by default
		default:
			break;
	}
}

void LTFAT_NAME(extend_right)(const LTFAT_REAL *in,int inLen, LTFAT_REAL *buffer, int filtLen, int type, int a){
	int legalExtLen = imin(filtLen-1, inLen);
	switch (type) {
		case 1: // half-point symmetry
           for(int ii=0;ii<legalExtLen;ii++)
			   buffer[ii] = in[legalExtLen-ii];
           break;
		case 2: // whole-point symmetry
			legalExtLen = imin(filtLen-1, inLen-1);
           for(int ii=0;ii<legalExtLen;ii++)
			   buffer[ii] = in[inLen-1-1-ii];
           break;
		 case 3: // half-point antisymmetry
           for(int ii=0;ii<legalExtLen;ii++)
			   buffer[ii] = -in[inLen-1-ii];
           break;
		 case 4: // whole-point antisymmetry
		   legalExtLen = imin(filtLen-1, inLen-1);
           for(int ii=0;ii<legalExtLen;ii++)
			   buffer[ii] = -in[inLen-1-1-ii];
           break;
		 case 5: // periodic padding
           for(int ii=0;ii<legalExtLen;ii++)
			   buffer[ii] = in[ii];
           break;
		 case 6: // constant padding
           for(int ii=0;ii<filtLen;ii++)
			   buffer[ii] = in[inLen-1];
           break;
		 case 7: // periodic padding with possible last sample repplications
			int rem = inLen%a;
			if(rem==0)
			{
              for(int ii=0;ii<legalExtLen;ii++)
			    buffer[ii] = in[ii];
			}
			else
			{
			  int remto = a - rem;
			  // replicated
			  for(int ii=0;ii<remto;ii++)
			     buffer[ii] = in[inLen-1];

			  // periodized
			  for(int ii=0;ii<legalExtLen-remto;ii++)
			     buffer[ii+remto] = in[ii];
			}
           break;
		case 0: // zero-padding by default
		default:
			break;
	}



}





