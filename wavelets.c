#include <string.h>
#include <math.h>
#include "config.h"
#include "ltfat.h"
#include "wavelets.h"

LTFAT_EXTERN
void LTFAT_NAME(convsub_td)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *out, const int outLen, const LTFAT_REAL *filts, int fLen, int sub, int skip, enum ltfatWavExtType ext)
{
    LTFAT_REAL *filtRev = (LTFAT_REAL *) ltfat_malloc(fLen*sizeof(LTFAT_REAL));
    for(unsigned int ii=0;ii<fLen;ii++)
    {
       *(filtRev+ii) = *(filts + fLen-1 - ii);
    }

	LTFAT_REAL *righExtbuff = 0;
	// number of output samples that can be calculated "painlessly"
    int outLenN = imax((inLen - skip + sub -1)/sub,0);

   // prepare cyclic buffer of length of power of two (for effective modulo operations)
   int buffLen = nextPow2(imax(fLen,sub+1));
   // buffer index
   int buffPtr = 0;

   // allocating and initializing the cyclic buffer
   LTFAT_REAL *buffer = (LTFAT_REAL *) ltfat_malloc(buffLen*sizeof(LTFAT_REAL));
   memset(buffer,0,buffLen*sizeof(LTFAT_REAL));

   // pointer for moving in the input data
   const LTFAT_REAL *tmpIn = in;
   LTFAT_REAL *tmpOut = out;
   LTFAT_REAL *tmpFilts = filts;
   LTFAT_REAL *tmpBuffPtr = buffer;

   // fill buffer with the initial values from the input signal according to the boundary treatment
   // last fLenUps buffer samples are filled to keep buffPtr=0
   extend_left(in,inLen,buffer,buffLen,fLen,ext,sub);

   if(outLenN<outLen)
   {
   	   // right extension is necessary, additional buffer from where to copy
	   righExtbuff = (LTFAT_REAL *) ltfat_malloc(buffLen*sizeof(LTFAT_REAL));
       memset(righExtbuff,0,buffLen*sizeof(LTFAT_REAL));
	   // store extension in the buffer (must be done now to avoid errors when inplace calculation is done)
	   extend_right(in,inLen,righExtbuff,fLen,ext,sub);
   }



//int idx = modPow2((-(jj*filtUpsPow2)+buffPtr-1),buffLen);
#define ONEOUTSAMPLE                                                    \
          tmpFilts = filtRev;                                           \
          int revBufPtr = modPow2(buffPtr-fLen,buffLen);                \
          int loop1it = fLen+1;                                         \
	      while(--loop1it)                                              \
	      {                                                             \
		     tmpBuffPtr = buffer + modPow2(revBufPtr++,buffLen);        \
             *tmpOut += *(tmpBuffPtr) * *(tmpFilts++);                  \
          }                                                             \
          tmpOut++;



#define READNEXTDATA(samples,wherePtr)                                              \
	   buffOver = imax(buffPtr+(samples)-buffLen, 0);                               \
   	   memcpy(buffer + buffPtr, wherePtr, ((samples)-buffOver)*sizeof(LTFAT_REAL)); \
	   memcpy(buffer,wherePtr+(samples)-buffOver,buffOver*sizeof(LTFAT_REAL));      \
	   buffPtr = modPow2(buffPtr += (samples),buffLen);


   int buffOver = 0;
   /*** initial buffer fill ***/
   int sampToRead = imin((skip+1),inLen);
   READNEXTDATA(sampToRead,tmpIn);
   tmpIn += sampToRead;

   /*********** STEP 1: FREE LUNCH ( but also a hot-spot) *******************************/
   // Take the smaller value from "painless" output length and the user defined output length
   int iiLoops = imin(outLenN-1,outLen-1);

   // loop trough all output samples, omit the very last one.
   for (int ii = 0; ii < iiLoops; ii++)
   {
       ONEOUTSAMPLE
       READNEXTDATA(sub,tmpIn)
       tmpIn += sub;
   }

  /*********** STEP 2: FINALIZE FREE LUNCH ************************************/
if(outLenN>0)
{
    ONEOUTSAMPLE
}
  /*********** STEP 3: NOW FOR THE TRICKY PART ************************************/
  if(outLenN<outLen)
  {
	  /************ STEP 3a: DEALING WITH THE REMAINING SAMPLES ******************/
	  // CAREFULL NOW! possibly stepping outside of input signal
	  // last index in the input signal for which reading next sub samples reaches outside of the input signal
	  int rightExtBuffIdx = 0;
	  if(outLenN>0)
      {
         int lastInIdx = (sub*(outLenN-1)+1+skip);
         rightExtBuffIdx = lastInIdx + sub - inLen;
         int diff = imax(0,inLen - lastInIdx);
         READNEXTDATA(diff,(in + lastInIdx))
      }
      else
      {
         rightExtBuffIdx = 1+skip - inLen;
      }

	   // now copying samples that are outside
	   READNEXTDATA(rightExtBuffIdx,righExtbuff)

	/************ STEP 3b: ALL OK, proceed reading input values from righExtbuff ******************/
	// loop for the remaining output samples
	for(int ii=0;ii<outLen-outLenN;ii++)
	{
	   ONEOUTSAMPLE
       READNEXTDATA(sub,(righExtbuff+rightExtBuffIdx))
       rightExtBuffIdx = modPow2(rightExtBuffIdx += sub,buffLen);
	}
  }


#undef READNEXTDATA
#undef ONEOUTSAMPLE
  if(righExtbuff)ltfat_free(righExtbuff);
   ltfat_free(buffer);
   ltfat_free(filtRev);
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
void LTFAT_NAME(conv_td_sub)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *out[], const int outLen, const LTFAT_REAL *filts[], int fLen, int noOfFilts, int sub, int skip,enum ltfatWavExtType ext, int filtUps)
{
	if(filtUps<1) filtUps=1;
	int filtUpsPow2 = filtUps;
	filtUps = log2(filtUpsPow2);
	int fLenUps = filtUpsPow2*fLen-(filtUpsPow2-1);

	LTFAT_REAL *righExtbuff = 0;
	// number of output samples that can be calculated "painlessly"
    int outLenN = imax((inLen - skip + sub -1)/sub,0);

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
//int idx = modPow2((-(jj*filtUpsPow2)+buffPtr-1),buffLen);
#define ONEOUTSAMPLE(outIdx)                                         \
	   for (int jj = 0; jj < fLen; jj++)                             \
	   {                                                             \
		  int idx = modPow2((-(jj*filtUpsPow2)+buffPtr-1),buffLen);  \
          LTFAT_REAL *buffTmp = buffer+idx;                          \
		  for(int kk = 0;kk<noOfFilts;kk++)                          \
		  {                                                          \
				outTmp[kk] += *buffTmp * filts[kk][jj];              \
		  }                                                          \
       }                                                             \
                                                                     \
	   for(int kk = 0;kk<noOfFilts;kk++)                             \
	   {                                                             \
		 out[kk][outIdx] = outTmp[kk];                               \
	   }                                                             \
	   memset(outTmp,0,noOfFilts*sizeof(LTFAT_REAL));


#define READNEXTDATA(samples,wherePtr)                                              \
	   buffOver = imax(buffPtr+(samples)-buffLen, 0);                               \
   	   memcpy(buffer + buffPtr, wherePtr, ((samples)-buffOver)*sizeof(LTFAT_REAL)); \
	   memcpy(buffer,wherePtr+(samples)-buffOver,buffOver*sizeof(LTFAT_REAL));      \
	   buffPtr = modPow2(buffPtr += (samples),buffLen);

   int buffOver = 0;
   /*** initial buffer fill ***/
   int sampToRead = imin((skip+1),inLen);
   READNEXTDATA(sampToRead,tmpIn);
   tmpIn += sampToRead;

   // accumulators for each output channels
   LTFAT_REAL *outTmp = (LTFAT_REAL *) ltfat_malloc(noOfFilts*sizeof(LTFAT_REAL));
   memset(outTmp,0,noOfFilts*sizeof(LTFAT_REAL));

   /*********** STEP 1: FREE LUNCH ( but also a hot-spot) *******************************/
   // Take the smaller value from "painless" output length and the user defined output length
   int iiLoops = imin(outLenN-1,outLen-1);

   // loop trough all output samples, omit the very last one.
   for (int ii = 0; ii < iiLoops; ii++)
   {
       ONEOUTSAMPLE(ii);
       READNEXTDATA(sub,tmpIn);
       tmpIn += sub;
   }

  /*********** STEP 2: FINALIZE FREE LUNCH ************************************/
if(outLenN>0)
{
    ONEOUTSAMPLE(outLenN-1);
}
  /*********** STEP 3: NOW FOR THE TRICKY PART ************************************/
  if(outLenN<outLen)
  {
	  /************ STEP 3a: DEALING WITH THE REMAINING SAMPLES ******************/
	  // CAREFULL NOW! possibly stepping outside of input signal
	  // last index in the input signal for which reading next sub samples reaches outside of the input signal
	  int rightExtBuffIdx = 0;
	  if(outLenN>0)
      {
         int lastInIdx = (sub*(outLenN-1)+1+skip);
         rightExtBuffIdx = lastInIdx + sub - inLen;
         int diff = imax(0,inLen - lastInIdx);
         READNEXTDATA(diff,(in + lastInIdx));
      }
      else
      {
         rightExtBuffIdx = 1+skip - inLen;
      }

	   // now copying samples that are outside
	   READNEXTDATA(rightExtBuffIdx,righExtbuff);

	/************ STEP 3b: ALL OK, proceed reading input values from righExtbuff ******************/
	// loop for the remaining output samples
	for(int ii=0;ii<outLen-outLenN;ii++)
	{
	   ONEOUTSAMPLE(outLenN+ii);
       READNEXTDATA(sub,(righExtbuff+rightExtBuffIdx));
       rightExtBuffIdx = modPow2(rightExtBuffIdx += sub,buffLen);
	}
  }


#undef READNEXTDATA
#undef ONEOUTSAMPLE
  if(righExtbuff)ltfat_free(righExtbuff);
   ltfat_free(outTmp);
   ltfat_free(buffer);
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



// fills last buffer samples
LTFAT_EXTERN
void LTFAT_NAME(extend_left)(const LTFAT_REAL *in, int inLen, LTFAT_REAL *buffer,int buffLen, int filtLen, enum ltfatWavExtType ext, int a){
		int legalExtLen = (filtLen-1)%inLen;
		int inLenTimes = (filtLen-1)/inLen;
		LTFAT_REAL *buffTmp = buffer + buffLen - legalExtLen;
	switch (ext) {
		case SYM: // half-point symmetry
		case EVEN:
           for(int ii=0;ii<legalExtLen;ii++)
			   buffTmp[ii] = in[legalExtLen-ii-1];
           break;
		case SYMW: // whole-point symmetry
			legalExtLen = imin(filtLen-1, inLen-1);
			buffTmp = buffer + buffLen - legalExtLen;
           for(int ii=0;ii<legalExtLen;ii++)
			   buffTmp[ii] = in[legalExtLen-ii];
           break;
		 case ASYM: // half-point antisymmetry
		 case ODD:
           for(int ii=0;ii<legalExtLen;ii++)
			   buffTmp[ii] = -in[legalExtLen-ii-1];
           break;
		 case ASYMW: // whole-point antisymmetry
			legalExtLen = imin(filtLen-1, inLen-1);
			legalExtLen = imin(filtLen-1, inLen-1);
			buffTmp = buffer + buffLen - legalExtLen;
           for(int ii=0;ii<legalExtLen;ii++)
			   buffTmp[ii] = -in[legalExtLen-ii];
           break;
		 case PPD: // periodic padding
         case PER:
             {
               LTFAT_REAL *bufferPtr = buffer + buffLen - (filtLen-1);
               for(int ii=0;ii<legalExtLen;ii++)
               {
  			      *(bufferPtr) = in[inLen-1-(legalExtLen-1)+ii];
  			      bufferPtr++;
               }

               for(int ii=0;ii<inLenTimes;ii++)
               {
                  for(int jj=0;jj<inLen;jj++)
                  {
                     *(bufferPtr) = in[jj];
                     bufferPtr++;
                  }
               }

	         }
           break;
		 case SP0: // constant padding
			buffTmp = buffer + buffLen - (filtLen-1);
           for(int ii=0;ii<filtLen-1;ii++)
			   buffTmp[ii] = in[0];
           break;
		   case PERDEC: // periodic padding with possible last sample repplication
		   {
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
           }
           break;
		case ZPD: // zero-padding by default
		case ZERO:
		default:
			break;
	}
}

void LTFAT_NAME(extend_right)(const LTFAT_REAL *in,int inLen, LTFAT_REAL *buffer, int filtLen, enum ltfatWavExtType ext, int a){
	int legalExtLen = (filtLen-1)%inLen;
	int inLenTimes = (filtLen-1)/inLen;
	switch (ext) {
		case SYM: // half-point symmetry
		case EVEN:
           for(int ii=0;ii<legalExtLen;ii++)
			   buffer[ii] = in[legalExtLen-ii];
           break;
		case SYMW: // whole-point symmetry
			legalExtLen = imin(filtLen-1, inLen-1);
           for(int ii=0;ii<legalExtLen;ii++)
			   buffer[ii] = in[inLen-1-1-ii];
           break;
		 case ASYM: // half-point antisymmetry
		 case ODD:
           for(int ii=0;ii<legalExtLen;ii++)
			   buffer[ii] = -in[inLen-1-ii];
           break;
		 case ASYMW: // whole-point antisymmetry
		   legalExtLen = imin(filtLen-1, inLen-1);
           for(int ii=0;ii<legalExtLen;ii++)
			   buffer[ii] = -in[inLen-1-1-ii];
           break;
		 case PPD: // periodic padding
         case PER:
             {
              LTFAT_REAL *bufferPtr = buffer;
              for(int ii=0;ii<inLenTimes;ii++)
              {
                for(int jj=0;jj<inLen;jj++)
                {
                   *(bufferPtr) = in[jj];
                   bufferPtr++;
                }
              }
              for(int ii=0;ii<legalExtLen;ii++)
              {
                 *(bufferPtr) = in[ii];
                 bufferPtr++;
              }
             }
           break;
		 case SP0: // constant padding
           for(int ii=0;ii<filtLen;ii++)
			   buffer[ii] = in[inLen-1];
           break;
		   case PERDEC: // periodic padding with possible last sample repplication
            {
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
            }
		case ZPD: // zero-padding by default
		case ZERO:
		default:
			break;
	}



}





