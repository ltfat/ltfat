/* NOT PROCESSED DIRECTLY, see ltfat_complexindependent.c */
#ifdef LTFAT_TYPE
#include <string.h>
#include <math.h>
#include "config.h"
#include "ltfat.h"

LTFAT_EXTERN
void LTFAT_NAME(atrousconvsub_td)(const LTFAT_TYPE *in, int inLen, LTFAT_TYPE *out, const int outLen, const LTFAT_TYPE *filts, int fLen, int filtUp, int skip, enum ltfatWavExtType ext)
{
    LTFAT_TYPE *filtRev = (LTFAT_TYPE *) ltfat_malloc(fLen*sizeof(LTFAT_TYPE));
    for(int ii=0;ii<fLen;ii++)
    {
       *(filtRev+ii) = *(filts + fLen-1 - ii);
    }
    int fLenUps = filtUp*fLen-(filtUp-1);

	LTFAT_TYPE *righExtbuff = 0;
	// number of output samples that can be calculated "painlessly"
    int outLenN = imax((inLen - skip),0);

   // prepare cyclic buffer of length of power of two (for effective modulo operations)
   int buffLen = nextPow2(fLenUps);
   // buffer index
   int buffPtr = 0;

   // allocating and initializing the cyclic buffer
   LTFAT_TYPE *buffer = (LTFAT_TYPE *) ltfat_malloc(buffLen*sizeof(LTFAT_TYPE));
   memset(buffer,0,buffLen*sizeof(LTFAT_TYPE));

   // pointer for moving in the input data
   const LTFAT_TYPE *tmpIn = in;
   LTFAT_TYPE *tmpOut = out;
   LTFAT_TYPE *tmpFilts = filtRev;
   LTFAT_TYPE *tmpBuffPtr = buffer;

   // fill buffer with the initial values from the input signal according to the boundary treatment
   // last fLenUps buffer samples are filled to keep buffPtr=0
   LTFAT_NAME(extend_left)(in,inLen,buffer,buffLen,fLenUps,ext,1);

   if(outLenN<outLen)
   {
   	   // right extension is necessary, additional buffer from where to copy
	   righExtbuff = (LTFAT_TYPE *) ltfat_malloc(buffLen*sizeof(LTFAT_TYPE));
       memset(righExtbuff,0,buffLen*sizeof(LTFAT_TYPE));
	   // store extension in the buffer (must be done now to avoid errors when inplace calculation is done)
	   LTFAT_NAME(extend_right)(in,inLen,righExtbuff,fLenUps,ext,1);
   }

#define ONEOUTSAMPLE                                                    \
          tmpFilts = filtRev;                                           \
          int revBufPtr = modPow2(buffPtr-fLenUps,buffLen);             \
          int loop1it = fLen+1;                                         \
	      while(--loop1it)                                              \
	      {                                                             \
		     tmpBuffPtr = buffer + modPow2(revBufPtr,buffLen);          \
		     revBufPtr+=filtUp;                                         \
             *tmpOut += *(tmpBuffPtr) * *(tmpFilts++);                  \
          }                                                             \
          tmpOut++;


#define READNEXTDATA(samples,wherePtr)                                              \
	   buffOver = imax(buffPtr+(samples)-buffLen, 0);                               \
   	   memcpy(buffer + buffPtr, wherePtr, ((samples)-buffOver)*sizeof(LTFAT_TYPE)); \
	   memcpy(buffer,wherePtr+(samples)-buffOver,buffOver*sizeof(LTFAT_TYPE));      \
	   buffPtr = modPow2(buffPtr += (samples),buffLen);

#define READNEXTSAMPLE(wherePtr)                               \
   	   *(buffer + buffPtr) = *wherePtr;                        \
	   buffPtr = modPow2(++buffPtr,buffLen);


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
       READNEXTSAMPLE(tmpIn)
       tmpIn++;
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
         int lastInIdx = ((outLenN-1)+1+skip);
         rightExtBuffIdx = lastInIdx + 1 - inLen;
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
       READNEXTSAMPLE((righExtbuff+rightExtBuffIdx))
       ++rightExtBuffIdx;
       //rightExtBuffIdx = modPow2(++rightExtBuffIdx,buffLen);
	}
  }


#undef READNEXTDATA
#undef READNEXTSAMPLE
#undef ONEOUTSAMPLE
  if(righExtbuff)ltfat_free(righExtbuff);
   ltfat_free(buffer);
   ltfat_free(filtRev);
}

LTFAT_EXTERN
void LTFAT_NAME(atrousupconv_td)(const LTFAT_TYPE *in, int inLen, LTFAT_TYPE *out, const int outLen, const LTFAT_TYPE *filts, int fLen, int filtUp, int skip, enum ltfatWavExtType ext)
{
   int fLenUps = filtUp*fLen-(filtUp-1);
   // Running output pointer
   LTFAT_TYPE* tmpOut = out;
   // Running input pointer
   LTFAT_TYPE* tmpIn =  (LTFAT_TYPE*) in;

   /** prepare cyclic buffer */
   int buffLen = nextPow2(fLenUps);
   LTFAT_TYPE* buffer = (LTFAT_TYPE *) ltfat_calloc(buffLen,sizeof(LTFAT_TYPE));
   int buffPtr = 0;

   int iiLoops = 0;
   int remainsOutSamp = outLen;
   int rightBuffPreLoad = 0;

   if(skip >= inLen)
   {
       rightBuffPreLoad = (skip + 1) - inLen;
       skip = inLen;
   }
   else
   {
      iiLoops = imin(inLen - skip,outLen); // just in case outLen < inLen - inSkip
      remainsOutSamp = outLen - (iiLoops-1);
   }

   LTFAT_TYPE *rightBuffer = (LTFAT_TYPE *) ltfat_calloc(buffLen,sizeof(LTFAT_TYPE));
   LTFAT_TYPE *rightBufferTmp = rightBuffer;

   if(ext==PER) // if periodic extension
   {
         LTFAT_NAME(extend_left)(in,inLen,buffer,buffLen,fLenUps,PER,0); // extension as a last (tmpfLen-1) samples of the buffer -> pointer dont have to be moved
		 LTFAT_NAME(extend_right)(in,inLen,rightBuffer,fLenUps,PER,0);
   }

   int iniStoCopy = imin(skip,buffLen);
   int tmpInSkip = imax(0,skip-buffLen);
   memcpy(buffer,tmpIn+tmpInSkip,iniStoCopy*sizeof(LTFAT_TYPE));
   tmpIn += (iniStoCopy+tmpInSkip);
   buffPtr = modPow2(buffPtr += iniStoCopy,buffLen);


 //LTFAT_TYPE* filtTmp = filts;
 #define ONEOUTSAMPLE(filtTmp,jjLoops)                                   \
	    for(int jj=0;jj<(jjLoops);jj++)                                  \
		    {                                                            \
				int idx = modPow2((-jj*filtUp+buffPtr-1), buffLen);      \
				*tmpOut += *(buffer+idx) * *((filtTmp) + jj);            \
		    }                                                            \
	    tmpOut++;

#define READNEXTSAMPLE(wherePtr)                               \
   	   *(buffer + buffPtr) = *(wherePtr);                      \
	   buffPtr = modPow2(++buffPtr,buffLen);


   /** STEP 2: MAIN LOOP */
   if(iiLoops>0)
   {
	 for(int ii=0;ii<iiLoops-1;ii++)
	 {
	     READNEXTSAMPLE(tmpIn)
	     tmpIn++;
		 ONEOUTSAMPLE(filts,fLen)
	 }
	 READNEXTSAMPLE(tmpIn)
	 //tmpIn++;
   }


     /** STEP 3b: load samples from right buffer */
	 while(rightBuffPreLoad--)
     {
		READNEXTSAMPLE((rightBufferTmp))
		rightBufferTmp++;
     }


     /*
	 STEP 3b: calculate remaining output samples,
	 Again, there can be shift/up misaligment thne shift>inLen
	 */

       for(int ii=0;ii<remainsOutSamp;ii++)
       {
          if(ii!=0)
		  {
		    READNEXTSAMPLE((rightBufferTmp))
		    rightBufferTmp++;
		  }
		  ONEOUTSAMPLE((filts),(fLen))
	   }

#undef READNEXTDATA
#undef ONEOUTSAMPLE
	ltfat_free(buffer);
	ltfat_free(rightBuffer);
}


LTFAT_EXTERN
void LTFAT_NAME(convsub_td)(const LTFAT_TYPE *in, int inLen, LTFAT_TYPE *out, const int outLen, const LTFAT_TYPE *filts, int fLen, int sub, int skip, enum ltfatWavExtType ext)
{
    LTFAT_TYPE *filtRev = (LTFAT_TYPE *) ltfat_malloc(fLen*sizeof(LTFAT_TYPE));
    for(int ii=0;ii<fLen;ii++)
    {
       *(filtRev+ii) = *(filts + fLen-1 - ii);
    }

	LTFAT_TYPE *righExtbuff = 0;
	// number of output samples that can be calculated "painlessly"
    int outLenN = imax((inLen - skip + sub -1)/sub,0);

   // prepare cyclic buffer of length of power of two (for effective modulo operations)
   int buffLen = nextPow2(imax(fLen,sub+1));
   // buffer index
   int buffPtr = 0;

   // allocating and initializing the cyclic buffer
   LTFAT_TYPE *buffer = (LTFAT_TYPE *) ltfat_malloc(buffLen*sizeof(LTFAT_TYPE));
   memset(buffer,0,buffLen*sizeof(LTFAT_TYPE));

   // pointer for moving in the input data
   const LTFAT_TYPE * tmpIn = in;
   LTFAT_TYPE * tmpOut = out;
   LTFAT_TYPE *tmpFilts = filtRev;
   LTFAT_TYPE *tmpBuffPtr = buffer;

   // fill buffer with the initial values from the input signal according to the boundary treatment
   // last fLenUps buffer samples are filled to keep buffPtr=0
   LTFAT_NAME(extend_left)(in,inLen,buffer,buffLen,fLen,ext,sub);

   if(outLenN<outLen)
   {
   	   // right extension is necessary, additional buffer from where to copy
	   righExtbuff = (LTFAT_TYPE *) ltfat_malloc(buffLen*sizeof(LTFAT_TYPE));
       memset(righExtbuff,0,buffLen*sizeof(LTFAT_TYPE));
	   // store extension in the buffer (must be done now to avoid errors when inplace calculation is done)
	   LTFAT_NAME(extend_right)(in,inLen,righExtbuff,fLen,ext,sub);
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
   	   memcpy(buffer + buffPtr, wherePtr, ((samples)-buffOver)*sizeof(LTFAT_TYPE)); \
	   memcpy(buffer,wherePtr+(samples)-buffOver,buffOver*sizeof(LTFAT_TYPE));      \
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


LTFAT_EXTERN
void LTFAT_NAME(upconv_td)(const LTFAT_TYPE *in, int inLen, LTFAT_TYPE *out, const int outLen, const LTFAT_TYPE *filts, int fLen, int up, int skip, enum ltfatWavExtType ext)
{
   // Copy, reverse and conjugate the imp resp.
   LTFAT_TYPE *filtsInv = (LTFAT_TYPE *) ltfat_malloc(fLen*sizeof(LTFAT_TYPE));
   memcpy(filtsInv,filts,fLen*sizeof(LTFAT_TYPE));
   LTFAT_NAME(reverse_array)(filtsInv,filtsInv,fLen);
   LTFAT_NAME(conjugate_array)(filtsInv,filtsInv,fLen);
   skip = -(1 - fLen + skip);

   // Running output pointer
   LTFAT_TYPE* tmpOut = out;
   // Running input pointer
   LTFAT_TYPE* tmpIn =  (LTFAT_TYPE*) in;

   /** prepare cyclic buffer */
   int buffLen = nextPow2(fLen);
   LTFAT_TYPE* buffer = (LTFAT_TYPE *) ltfat_calloc(buffLen,sizeof(LTFAT_TYPE));
   int buffPtr = 0;



   int inSkip = (skip + up - 1)/up;
   int skipModUp = skip%up;
   int skipToNextUp = 0;
   if(skipModUp!=0)  skipToNextUp = up-skipModUp;
   int outAlign = 0;

   int iiLoops = 0;
   int uuLoops = 0;
   int remainsOutSamp = outLen;
   int rightBuffPreLoad = 0;

   if(inSkip >= inLen)
   {
       inSkip = inLen;
       outAlign = skipModUp;
       rightBuffPreLoad = (skip + 1 + up - 1)/up - inLen;
   }
   else
   {
      uuLoops = skipToNextUp;
      iiLoops = imin(inLen - inSkip,(outLen-skipToNextUp + up -1)/up); // just in case outLen/up < inLen - inSkip
      remainsOutSamp = outLen - (uuLoops + (iiLoops-1)*up);
   }

   LTFAT_TYPE *rightBuffer = (LTFAT_TYPE *) ltfat_calloc(buffLen,sizeof(LTFAT_TYPE));
   LTFAT_TYPE *rightBufferTmp = rightBuffer;

   if(ext==PER) // if periodic extension
   {
         LTFAT_NAME(extend_left)(in,inLen,buffer,buffLen,fLen,PER,0); // extension as a last (tmpfLen-1) samples of the buffer -> pointer dont have to be moved
		 LTFAT_NAME(extend_right)(in,inLen,rightBuffer,fLen,PER,0);
   }

   int iniStoCopy = imin(inSkip,buffLen);
   int tmpInSkip = imax(0,inSkip-buffLen);
   memcpy(buffer,tmpIn+tmpInSkip,iniStoCopy*sizeof(LTFAT_TYPE));
   tmpIn += (iniStoCopy+tmpInSkip);
   buffPtr = modPow2(buffPtr += iniStoCopy,buffLen);


 //LTFAT_TYPE* filtTmp = filts;
 #define ONEOUTSAMPLE(filtTmp,jjLoops)                                   \
	    for(int jj=0;jj<(jjLoops);jj++)                                  \
		    {                                                            \
				int idx = modPow2((-jj+buffPtr-1), buffLen);             \
				*tmpOut += *(buffer+idx) * *((filtTmp) +(jj*up));        \
		    }                                                            \
	    tmpOut++;

#define READNEXTSAMPLE(wherePtr)                               \
   	   *(buffer + buffPtr) = *(wherePtr);                      \
	   buffPtr = modPow2(++buffPtr,buffLen);


   /** STEP 1: Deal with the shift - upsampling misaligment */
   for(int uu=0;uu<uuLoops;uu++)
   {
       ONEOUTSAMPLE((filtsInv + skipModUp+uu),((fLen-(skipModUp+uu)+up-1)/up))
   }

   /** STEP 2: MAIN LOOP */
   if(iiLoops>0)
   {
	 for(int ii=0;ii<iiLoops-1;ii++)
	 {
	     READNEXTSAMPLE(tmpIn)
	     tmpIn++;
		 for(int uu=0;uu<up;uu++)
	     {
			ONEOUTSAMPLE((filtsInv+uu),((fLen-uu+up-1)/up))
		 }
	 }
	 READNEXTSAMPLE(tmpIn)
	 tmpIn++;
   }


     /** STEP 3b: load samples from right buffer */
	 while(rightBuffPreLoad--)
     {
		READNEXTSAMPLE((rightBufferTmp))
		rightBufferTmp++;
     }


     /*
	 STEP 3b: calculate remaining output samples,
	 Again, there can be shift/up misaligment thne shift>inLen
	 */

       for(int ii=outAlign;ii<remainsOutSamp+outAlign;ii++)
       {
          if(ii!=outAlign&&ii%up==0)
		  {
		    READNEXTSAMPLE((rightBufferTmp))
		    rightBufferTmp++;
		  }
		  ONEOUTSAMPLE((filtsInv+ii%up),((fLen-ii%up+up-1)/up))
	   }

    #undef ONEOUTSAMPLE
    #undef READNEXTSAMPLE
	ltfat_free(buffer);
	ltfat_free(rightBuffer);
	ltfat_free(filtsInv);
}



// fills last buffer samples
LTFAT_EXTERN
void LTFAT_NAME(extend_left)(const LTFAT_TYPE *in, int inLen, LTFAT_TYPE *buffer,int buffLen, int filtLen, enum ltfatWavExtType ext, int a){
		int legalExtLen = (filtLen-1)%inLen;
		int inLenTimes = (filtLen-1)/inLen;
		LTFAT_TYPE *buffTmp = buffer + buffLen - legalExtLen;
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
               LTFAT_TYPE *bufferPtr = buffer + buffLen - (filtLen-1);
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
		case VALID:
		default:
			break;
	}
}

void LTFAT_NAME(extend_right)(const LTFAT_TYPE *in,int inLen, LTFAT_TYPE *buffer, int filtLen, enum ltfatWavExtType ext, int a){
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
              LTFAT_TYPE *bufferPtr = buffer;
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
		case VALID:
		default:
			break;
	}



}



#endif // LTFAT_TYPE
