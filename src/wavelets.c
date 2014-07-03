/* NOT PROCESSED DIRECTLY, see ltfat_complexindependent.c */
#ifdef LTFAT_TYPE
#include "ltfat.h"
#include "ltfat_types.h"


LTFAT_EXTERN void
LTFAT_NAME(atrousfilterbank_td)(const LTFAT_TYPE *f, const LTFAT_TYPE *g[],
                                const ltfatInt L, const ltfatInt gl[],
                                const ltfatInt W, const ltfatInt a[],
                                const ltfatInt skip[], const ltfatInt M,
                                LTFAT_TYPE *c, ltfatExtType ext)
{
    for(ltfatInt m=0; m<M; m++)
    {
        for(ltfatInt w=0; w<W; w++)
        {
            LTFAT_NAME(atrousconvsub_td)(f+w*L, g[m], L, gl[m], a[m],
                                   skip[m],c + w*M*L + m*L, ext);
        }
    }
}

LTFAT_EXTERN void
LTFAT_NAME(iatrousfilterbank_td)(const LTFAT_TYPE *c, const LTFAT_TYPE *g[],
                                 const ltfatInt L, const ltfatInt gl[],
                                 const ltfatInt W, const ltfatInt a[],
                                 const ltfatInt skip[], const ltfatInt M,
                                 LTFAT_TYPE *f, ltfatExtType ext)
{
   // Set output array to zeros, since the array is used as an accumulator
    memset(f,0,L*W*sizeof*f);

    for(ltfatInt m=0; m<M; m++)
    {
        for(ltfatInt w=0; w<W; w++)
        {
            LTFAT_NAME(atrousupconv_td)(c + w*M*L + m*L, g[m], L, gl[m], a[m],
                                  skip[m],f+w*L, ext);
        }
    }

}


LTFAT_EXTERN void
LTFAT_NAME(filterbank_td)(const LTFAT_TYPE *f, const LTFAT_TYPE *g[],
                          const ltfatInt L, const ltfatInt gl[],
                          const ltfatInt W, const ltfatInt a[],
                          const ltfatInt skip[], const ltfatInt M,
                          LTFAT_TYPE *c[], ltfatExtType ext)
{
    for(ltfatInt m=0; m<M; m++)
    {
        const ltfatInt N = filterbank_td_size(L,a[m],gl[m],skip[m],ext);
        for(ltfatInt w=0; w<W; w++)
        {

            LTFAT_NAME(convsub_td)(f+w*L, g[m], L, gl[m], a[m],
                                   skip[m],c[m]+w*N, ext);
        }
    }
}


LTFAT_EXTERN void
LTFAT_NAME(ifilterbank_td)(const LTFAT_TYPE *c[], const LTFAT_TYPE *g[],
                           const ltfatInt L, const ltfatInt gl[],
                           const ltfatInt W, const ltfatInt a[],
                           const ltfatInt skip[], const ltfatInt M,
                           LTFAT_TYPE *f, ltfatExtType ext)
{
    memset(f,0,L*W*sizeof*f);

    for(ltfatInt m=0; m<M; m++)
    {
        const ltfatInt N = filterbank_td_size(L,a[m],gl[m],skip[m],ext);
        for(ltfatInt w=0; w<W; w++)
        {

            LTFAT_NAME(upconv_td)(c[m]+w*N, g[m], L, gl[m], a[m],
                                  skip[m],f+w*L, ext);
        }
    }

}


LTFAT_EXTERN void
LTFAT_NAME(atrousconvsub_td)(const LTFAT_TYPE *f, const LTFAT_TYPE *g,
                             const ltfatInt L, const ltfatInt gl, const ltfatInt ga,
                             ltfatInt skip, LTFAT_TYPE *c, ltfatExtType ext)
{
    memset(c,0,L*sizeof*c);
    ltfatInt skipLoc = -skip;
    LTFAT_TYPE *filtRev = ltfat_malloc(gl*sizeof*filtRev);
    LTFAT_NAME(reverse_array)((LTFAT_TYPE*)g,filtRev,gl);

    ltfatInt glUps = ga*gl-(ga-1);

    LTFAT_TYPE *righExtbuff = 0;
    // number of output samples that can be calculated "painlessly"
    ltfatInt Nsafe = imax((L - skipLoc),0);

    // prepare cyclic buf of length of power of two (for effective modulo operations)
    ltfatInt bufgl = nextPow2(glUps);
    // buf index
    ltfatInt buffPtr = 0;

    // allocating and initializing the cyclic buf
    LTFAT_TYPE *buf = ltfat_calloc(bufgl,sizeof*buf);

    // pointer for moving in the input data
    const LTFAT_TYPE *tmpIn = f;
    LTFAT_TYPE *tmpOut = c;
    LTFAT_TYPE *tmpg = filtRev;
    LTFAT_TYPE *tmpBuffPtr = buf;

    // fill buf with the initial values from the input signal according to the boundary treatment
    // last glUps buf samples are filled to keep buffPtr=0
    LTFAT_NAME(extend_left)(f,L,buf,bufgl,glUps,ext,1);

    if(Nsafe<L)
    {
        // right extension is necessary, additional buf from where to copy
        righExtbuff = ltfat_malloc(bufgl*sizeof(LTFAT_TYPE));
        memset(righExtbuff,0,bufgl*sizeof(LTFAT_TYPE));
        // store extension in the buf (must be done now to avoid errors when inplace calculation is done)
        LTFAT_NAME(extend_right)(f,L,righExtbuff,glUps,ext,1);
    }

#define ONEOUTSAMPLE                                                    \
         tmpg = filtRev;                                           \
         ltfatInt revBufPtr = modPow2(buffPtr-glUps,bufgl);             \
         ltfatInt loop1it = gl+1;                                         \
	      while(--loop1it)                                              \
	      {                                                             \
		     tmpBuffPtr = buf + modPow2(revBufPtr,bufgl);          \
		     revBufPtr+=ga;                                         \
           *tmpOut += *(tmpBuffPtr) * *(tmpg++);                  \
         }                                                             \
         tmpOut++;


#define READNEXTDATA(samples,wherePtr)                                              \
	   buffOver = imax(buffPtr+(samples)-bufgl, 0);                               \
   	memcpy(buf + buffPtr, wherePtr, ((samples)-buffOver)*sizeof(LTFAT_TYPE)); \
	   memcpy(buf,wherePtr+(samples)-buffOver,buffOver*sizeof(LTFAT_TYPE));      \
	   buffPtr = modPow2(buffPtr += (samples),bufgl);

#define READNEXTSAMPLE(wherePtr)                               \
   	   *(buf + buffPtr) = *wherePtr;                        \
	   buffPtr = modPow2(++buffPtr,bufgl);


    ltfatInt buffOver = 0;
    /*** initial buf fill ***/
    ltfatInt sampToRead = imin((skipLoc+1),L);
    READNEXTDATA(sampToRead,tmpIn);
    tmpIn += sampToRead;

    /*********** STEP 1: FREE LUNCH ( but also a hot-spot) *******************************/
    // Take the smaller value from "painless" output length and the user defined output length
    ltfatInt iiLoops = imin(Nsafe-1,L-1);

    // loop trough all output samples, omit the very last one.
    for (ltfatInt ii = 0; ii < iiLoops; ii++)
    {
        ONEOUTSAMPLE
        READNEXTSAMPLE(tmpIn)
        tmpIn++;
    }

    /*********** STEP 2: FINALIZE FREE LUNCH ************************************/
    if(Nsafe>0)
    {
        ONEOUTSAMPLE
    }
    /*********** STEP 3: NOW FOR THE TRICKY PART ************************************/
    if(Nsafe<L)
    {
        /************ STEP 3a: DEALING WITH THE REMAINING SAMPLES ******************/
        // CAREFULL NOW! possibly stepping outside of input signal
        // last index in the input signal for which reading next a samples reaches outside of the input signal
        ltfatInt rightExtBuffIdx = 0;
        if(Nsafe>0)
        {
            ltfatInt lastInIdx = ((Nsafe-1)+1+skipLoc);
            rightExtBuffIdx = lastInIdx + 1 - L;
            ltfatInt diff = imax(0,L - lastInIdx);
            READNEXTDATA(diff,(f + lastInIdx))
        }
        else
        {
            rightExtBuffIdx = 1+skipLoc - L;
        }

        // now copying samples that are outside
        READNEXTDATA(rightExtBuffIdx,righExtbuff)

        /************ STEP 3b: ALL OK, proceed reading input values from righExtbuff ******************/
        // loop for the remaining output samples
        for(ltfatInt ii=0; ii<L-Nsafe; ii++)
        {
            ONEOUTSAMPLE
            READNEXTSAMPLE((righExtbuff+rightExtBuffIdx))
            ++rightExtBuffIdx;
            //rightExtBuffIdx = modPow2(++rightExtBuffIdx,bufgl);
        }
    }


#undef READNEXTDATA
#undef READNEXTSAMPLE
#undef ONEOUTSAMPLE
    LTFAT_SAFEFREEALL(buf,filtRev,righExtbuff);
}

LTFAT_EXTERN void
LTFAT_NAME(atrousupconv_td)(const LTFAT_TYPE *c, const LTFAT_TYPE *g,
                            const ltfatInt L, const ltfatInt gl,
                            const ltfatInt ga, const ltfatInt skip,
                            LTFAT_TYPE *f, ltfatExtType ext)
{
    ltfatInt glUps = ga*gl-(ga-1);
    ltfatInt skipLoc = -(1-glUps-skip);

    // Copy, reverse and conjugate the imp resp.
    LTFAT_TYPE *gInv = ltfat_malloc(gl*sizeof*gInv);
    memcpy(gInv,g,gl*sizeof*gInv);
    LTFAT_NAME(reverse_array)(gInv,gInv,gl);
    LTFAT_NAME(conjugate_array)(gInv,gInv,gl);


    // Running output pointer
    LTFAT_TYPE* tmpOut = f;
    // Running input pointer
    LTFAT_TYPE* tmpIn =  (LTFAT_TYPE*) c;

    /** prepare cyclic buf */
    ltfatInt bufgl = nextPow2(glUps);
    LTFAT_TYPE* buf = ltfat_calloc(bufgl,sizeof*buf);
    ltfatInt buffPtr = 0;

    ltfatInt iiLoops = 0;
    ltfatInt remainsOutSamp = L;
    ltfatInt rightBuffPreLoad = 0;

    if(skipLoc >= L)
    {
        rightBuffPreLoad = (skipLoc + 1) - L;
        skipLoc = L;
    }
    else
    {
        iiLoops = imin(L - skipLoc,L); // just in case L < L - inSkip
        remainsOutSamp = L - (iiLoops-1);
    }

    LTFAT_TYPE *rightbuf = ltfat_calloc(bufgl,sizeof*rightbuf);
    LTFAT_TYPE *rightbufTmp = rightbuf;

    if(ext==PER) // if periodic extension
    {
        LTFAT_NAME(extend_left)(c,L,buf,bufgl,glUps,PER,0); // extension as a last (tmpgl-1) samples of the buf -> pointer dont have to be moved
        LTFAT_NAME(extend_right)(c,L,rightbuf,glUps,PER,0);
    }

    ltfatInt iniStoCopy = imin(skipLoc,bufgl);
    ltfatInt tmpInSkip = imax(0,skipLoc-bufgl);
    memcpy(buf,tmpIn+tmpInSkip,iniStoCopy*sizeof*buf);
    tmpIn += (iniStoCopy+tmpInSkip);
    buffPtr = modPow2(buffPtr += iniStoCopy,bufgl);


//LTFAT_TYPE* filtTmp = g;
#define ONEOUTSAMPLE(filtTmp,jjLoops)                                   \
	    for(ltfatInt jj=0;jj<(jjLoops);jj++)                                  \
		    {                                                            \
				ltfatInt idx = modPow2((-jj*ga+buffPtr-1), bufgl);      \
				*tmpOut += *(buf+idx) * *((filtTmp) + jj);            \
		    }                                                            \
	    tmpOut++;

#define READNEXTSAMPLE(wherePtr)                               \
   	   *(buf + buffPtr) = *(wherePtr);                      \
	   buffPtr = modPow2(++buffPtr,bufgl);


    /** STEP 2: MAIN LOOP */
    if(iiLoops>0)
    {
        for(ltfatInt ii=0; ii<iiLoops-1; ii++)
        {
            READNEXTSAMPLE(tmpIn)
            tmpIn++;
            ONEOUTSAMPLE(gInv,gl)
        }
        READNEXTSAMPLE(tmpIn)
        //tmpIn++;
    }


    /** STEP 3b: load samples from right buf */
    while(rightBuffPreLoad--)
    {
        READNEXTSAMPLE((rightbufTmp))
        rightbufTmp++;
    }


    /*
    STEP 3b: calculate remaining output samples,
    Again, there can be shift/up misaligment thne shift>L
    */

    for(ltfatInt ii=0; ii<remainsOutSamp; ii++)
    {
        if(ii!=0)
        {
            READNEXTSAMPLE((rightbufTmp))
            rightbufTmp++;
        }
        ONEOUTSAMPLE((gInv),(gl))
    }

#undef READNEXTDATA
#undef ONEOUTSAMPLE
    LTFAT_SAFEFREEALL(buf,rightbuf,gInv);
}


LTFAT_EXTERN void
LTFAT_NAME(convsub_td)(const LTFAT_TYPE *f, const LTFAT_TYPE *g, const ltfatInt L,
                       const ltfatInt gl, const ltfatInt a, const ltfatInt skip,
                       LTFAT_TYPE *c, ltfatExtType ext)
{
    const ltfatInt N = filterbank_td_size(L,a,gl,skip,ext);
    // Since c is used as an accu
    memset(c,0,N*sizeof*c);
    // Reverse and conjugate the filter
    LTFAT_TYPE *filtRev = ltfat_malloc(gl*sizeof*filtRev);
    LTFAT_NAME(reverse_array)((LTFAT_TYPE*)g, filtRev,gl);

    LTFAT_TYPE *righExtbuff = 0;
    // number of output samples that can be calculated "painlessly"
    ltfatInt Nsafe = imax((L + skip + a -1)/a,0);

    // prepare cyclic buf of length of power of two (for effective modulo operations)
    ltfatInt bufgl = nextPow2(imax(gl,a+1));
    // buf index
    ltfatInt buffPtr = 0;

    // allocating and initializing the cyclic buf
    LTFAT_TYPE *buf = ltfat_calloc(bufgl,sizeof*buf);

    // pointer for moving in the input data
    const LTFAT_TYPE * tmpIn = f;
    LTFAT_TYPE * tmpOut = c;
    LTFAT_TYPE *tmpg = filtRev;
    LTFAT_TYPE *tmpBuffPtr = buf;

    // fill buf with the initial values from the input signal according to the boundary treatment
    // last glUps buf samples are filled to keep buffPtr=0
    LTFAT_NAME(extend_left)(f,L,buf,bufgl,gl,ext,a);

    if(Nsafe<N)
    {
        // right extension is necessary, additional buf from where to copy
        righExtbuff = ltfat_calloc(bufgl,sizeof*righExtbuff);
        // store extension in the buf (must be done now to avoid errors when inplace calculation is done)
        LTFAT_NAME(extend_right)(f,L,righExtbuff,gl,ext,a);
    }

#define ONEOUTSAMPLE                                                    \
          tmpg = filtRev;                                           \
          ltfatInt revBufPtr = modPow2(buffPtr-gl,bufgl);                \
          ltfatInt loop1it = gl+1;                                         \
	      while(--loop1it)                                              \
	      {                                                             \
		     tmpBuffPtr = buf + modPow2(revBufPtr++,bufgl);        \
             *tmpOut += *(tmpBuffPtr) * *(tmpg++);                  \
          }                                                             \
          tmpOut++;



#define READNEXTDATA(samples,wherePtr)                                              \
	   buffOver = imax(buffPtr+(samples)-bufgl, 0);                               \
   	memcpy(buf + buffPtr, wherePtr, ((samples)-buffOver)*sizeof*buf); \
	   memcpy(buf,wherePtr+(samples)-buffOver,buffOver*sizeof*buf);      \
	   buffPtr = modPow2(buffPtr += (samples),bufgl);


    ltfatInt buffOver = 0;
    /*** initial buf fill ***/
    ltfatInt sampToRead = imin((-skip+1),L);
    READNEXTDATA(sampToRead,tmpIn);
    tmpIn += sampToRead;

    /*********** STEP 1: FREE LUNCH ( but also a hot-spot) *******************************/
    // Take the smaller value from "painless" output length and the user defined output length
    ltfatInt iiLoops = imin(Nsafe-1,N-1);

    // loop trough all output samples, omit the very last one.
    for (ltfatInt ii = 0; ii < iiLoops; ii++)
    {
        ONEOUTSAMPLE
        READNEXTDATA(a,tmpIn)
        tmpIn += a;
    }

    /*********** STEP 2: FINALIZE FREE LUNCH ************************************/
    if(Nsafe>0)
    {
        ONEOUTSAMPLE
    }
    /*********** STEP 3: NOW FOR THE TRICKY PART ************************************/
    if(Nsafe<N)
    {
        /************ STEP 3a: DEALING WITH THE REMAINING SAMPLES ******************/
        // CAREFULL NOW! possibly stepping outside of input signal
        // last index in the input signal for which reading next a samples reaches outside of the input signal
        ltfatInt rightExtBuffIdx = 0;
        if(Nsafe>0)
        {
            ltfatInt lastInIdx = (a*(Nsafe-1)+1-skip);
            rightExtBuffIdx = lastInIdx + a - L;
            ltfatInt diff = imax(0,L - lastInIdx);
            READNEXTDATA(diff,(f + lastInIdx))
        }
        else
        {
            rightExtBuffIdx = 1-skip - L;
        }

        // now copying samples that are outside
        READNEXTDATA(rightExtBuffIdx,righExtbuff)

        /************ STEP 3b: ALL OK, proceed reading input values from righExtbuff ******************/
        // loop for the remaining output samples
        for(ltfatInt ii=0; ii<N-Nsafe; ii++)
        {
            ONEOUTSAMPLE
            READNEXTDATA(a,(righExtbuff+rightExtBuffIdx))
            rightExtBuffIdx = modPow2(rightExtBuffIdx += a,bufgl);
        }
    }


#undef READNEXTDATA
#undef ONEOUTSAMPLE
    LTFAT_SAFEFREEALL(buf,filtRev,righExtbuff);
}


LTFAT_EXTERN void
LTFAT_NAME(upconv_td)(const LTFAT_TYPE *c, const LTFAT_TYPE *g, const ltfatInt L,
                      const ltfatInt gl, const ltfatInt a, const ltfatInt skip,
                      LTFAT_TYPE *f, ltfatExtType ext)
{
    const ltfatInt N = filterbank_td_size(L,a,gl,skip,ext);

    // Copy, reverse and conjugate the imp resp.
    LTFAT_TYPE *gInv = ltfat_malloc(gl*sizeof*gInv);
    memcpy(gInv,g,gl*sizeof*gInv);
    LTFAT_NAME(reverse_array)(gInv,gInv,gl);
    LTFAT_NAME(conjugate_array)(gInv,gInv,gl);
    ltfatInt skipRev = -(1 - gl - skip);

    // Running output pointer
    LTFAT_TYPE* tmpOut = f;
    // Running input pointer
    const LTFAT_TYPE* tmpIn =  c;

    /** prepare cyclic buf */
    ltfatInt bufgl = nextPow2(gl);
    LTFAT_TYPE* buf = ltfat_calloc(bufgl,sizeof*buf);
    ltfatInt buffPtr = 0;

    ltfatInt inSkip = (skipRev + a - 1)/a;
    ltfatInt skipModUp = skipRev%a;
    ltfatInt skipToNextUp = 0;
    if(skipModUp!=0)  skipToNextUp = a-skipModUp;
    ltfatInt outAlign = 0;

    ltfatInt iiLoops = 0;
    ltfatInt uuLoops = 0;
    ltfatInt remainsOutSamp = L;
    ltfatInt rightBuffPreLoad = 0;

    if(inSkip >= N)
    {
        inSkip = N;
        outAlign = skipModUp;
        rightBuffPreLoad = (skipRev + 1 + a - 1)/a - N;
    }
    else
    {
        uuLoops = skipToNextUp;
        iiLoops = imin(N - inSkip,(L-skipToNextUp + a -1)/a); // just in case L/a < N - inSkip
        remainsOutSamp = L - (uuLoops + (iiLoops-1)*a);
    }

    LTFAT_TYPE *rightbuf = ltfat_calloc(bufgl,sizeof(LTFAT_TYPE));
    LTFAT_TYPE *rightbufTmp = rightbuf;

    if(ext==PER) // if periodic extension
    {
        LTFAT_NAME(extend_left)(c,N,buf,bufgl,gl,PER,0); // extension as a last (tmpgl-1) samples of the buf -> pointer dont have to be moved
        LTFAT_NAME(extend_right)(c,N,rightbuf,gl,PER,0);
    }

    ltfatInt iniStoCopy = imin(inSkip,bufgl);
    ltfatInt tmpInSkip = imax(0,inSkip-bufgl);
    memcpy(buf,tmpIn+tmpInSkip,iniStoCopy*sizeof*buf);
    tmpIn += (iniStoCopy+tmpInSkip);
    buffPtr = modPow2(buffPtr += iniStoCopy,bufgl);



#define ONEOUTSAMPLE(filtTmp,jjLoops)                                   \
	    for(ltfatInt jj=0;jj<(jjLoops);jj++)                                  \
		    {                                                            \
				ltfatInt idx = modPow2((-jj+buffPtr-1), bufgl);             \
				*tmpOut += *(buf+idx) * *((filtTmp) +(jj*a));        \
		    }                                                            \
	    tmpOut++;

#define READNEXTSAMPLE(wherePtr)                               \
   	   *(buf + buffPtr) = *(wherePtr);                      \
	   buffPtr = modPow2(++buffPtr,bufgl);


    /** STEP 1: Deal with the shift - upsampling misaligment */
    for(ltfatInt uu=0; uu<uuLoops; uu++)
    {
        ONEOUTSAMPLE((gInv + skipModUp+uu),((gl-(skipModUp+uu)+a-1)/a))
    }

    /** STEP 2: MAIN LOOP */
    if(iiLoops>0)
    {
        for(ltfatInt ii=0; ii<iiLoops-1; ii++)
        {
            READNEXTSAMPLE(tmpIn)
            tmpIn++;
            for(ltfatInt uu=0; uu<a; uu++)
            {
                ONEOUTSAMPLE((gInv+uu),((gl-uu+a-1)/a))
            }
        }
        READNEXTSAMPLE(tmpIn)
        tmpIn++;
    }


    /** STEP 3b: load samples from right buf */
    while(rightBuffPreLoad--)
    {
        READNEXTSAMPLE((rightbufTmp))
        rightbufTmp++;
    }


    /*
    STEP 3b: calculate remaining output samples,
    Again, there can be shift/a misaligment thne shift>L
    */

    for(ltfatInt ii=outAlign; ii<remainsOutSamp+outAlign; ii++)
    {
        if(ii!=outAlign&&ii%a==0)
        {
            READNEXTSAMPLE((rightbufTmp))
            rightbufTmp++;
        }
        ONEOUTSAMPLE((gInv+ii%a),((gl-ii%a+a-1)/a))
    }

#undef ONEOUTSAMPLE
#undef READNEXTSAMPLE
    LTFAT_SAFEFREEALL(buf,rightbuf,gInv);
}





// fills last buf samples
LTFAT_EXTERN
void LTFAT_NAME(extend_left)(const LTFAT_TYPE *in, ltfatInt L, LTFAT_TYPE *buf,ltfatInt bufgl, ltfatInt gl, ltfatExtType ext, ltfatInt a)
{
    ltfatInt legalExtLen = (gl-1)%L;
    ltfatInt LTimes = (gl-1)/L;
    LTFAT_TYPE *buffTmp = buf + bufgl - legalExtLen;
    switch (ext)
    {
    case SYM: // half-point symmetry
    case EVEN:
        for(ltfatInt ii=0; ii<legalExtLen; ii++)
            buffTmp[ii] = in[legalExtLen-ii-1];
        break;
    case SYMW: // whole-point symmetry
        legalExtLen = imin(gl-1, L-1);
        buffTmp = buf + bufgl - legalExtLen;
        for(ltfatInt ii=0; ii<legalExtLen; ii++)
            buffTmp[ii] = in[legalExtLen-ii];
        break;
    case ASYM: // half-point antisymmetry
    case ODD:
        for(ltfatInt ii=0; ii<legalExtLen; ii++)
            buffTmp[ii] = -in[legalExtLen-ii-1];
        break;
    case ASYMW: // whole-point antisymmetry
        legalExtLen = imin(gl-1, L-1);
        legalExtLen = imin(gl-1, L-1);
        buffTmp = buf + bufgl - legalExtLen;
        for(ltfatInt ii=0; ii<legalExtLen; ii++)
            buffTmp[ii] = -in[legalExtLen-ii];
        break;
    case PPD: // periodic padding
    case PER:
    {
        LTFAT_TYPE *bufPtr = buf + bufgl - (gl-1);
        for(ltfatInt ii=0; ii<legalExtLen; ii++)
        {
            *(bufPtr) = in[L-1-(legalExtLen-1)+ii];
            bufPtr++;
        }

        for(ltfatInt ii=0; ii<LTimes; ii++)
        {
            for(ltfatInt jj=0; jj<L; jj++)
            {
                *(bufPtr) = in[jj];
                bufPtr++;
            }
        }

    }
    break;
    case SP0: // constant padding
        buffTmp = buf + bufgl - (gl-1);
        for(ltfatInt ii=0; ii<gl-1; ii++)
            buffTmp[ii] = in[0];
        break;
    case PERDEC: // periodic padding with possible last sample repplication
    {
        ltfatInt rem = L%a;
        if(rem==0)
        {
            for(ltfatInt ii=0; ii<legalExtLen; ii++)
                buffTmp[ii] = in[L-1-(legalExtLen-1)+ii];
        }
        else
        {
            ltfatInt remto = a - rem;

            // replicated
            for(ltfatInt ii=0; ii<remto; ii++)
                buffTmp[legalExtLen-1-ii] = in[L-1];

            // periodic extension
            for(ltfatInt ii=0; ii<legalExtLen-remto; ii++)
                buffTmp[ii] = in[L-1-(legalExtLen-1-1)+ii+remto-1];
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

void LTFAT_NAME(extend_right)(const LTFAT_TYPE *in,ltfatInt L, LTFAT_TYPE *buf, ltfatInt gl, ltfatExtType ext, ltfatInt a)
{
    ltfatInt legalExtLen = (gl-1)%L;
    ltfatInt LTimes = (gl-1)/L;
    switch (ext)
    {
    case SYM: // half-point symmetry
    case EVEN:
        for(ltfatInt ii=0; ii<legalExtLen; ii++)
            buf[ii] = in[legalExtLen-ii];
        break;
    case SYMW: // whole-point symmetry
        legalExtLen = imin(gl-1, L-1);
        for(ltfatInt ii=0; ii<legalExtLen; ii++)
            buf[ii] = in[L-1-1-ii];
        break;
    case ASYM: // half-point antisymmetry
    case ODD:
        for(ltfatInt ii=0; ii<legalExtLen; ii++)
            buf[ii] = -in[L-1-ii];
        break;
    case ASYMW: // whole-point antisymmetry
        legalExtLen = imin(gl-1, L-1);
        for(ltfatInt ii=0; ii<legalExtLen; ii++)
            buf[ii] = -in[L-1-1-ii];
        break;
    case PPD: // periodic padding
    case PER:
    {
        LTFAT_TYPE *bufPtr = buf;
        for(ltfatInt ii=0; ii<LTimes; ii++)
        {
            for(ltfatInt jj=0; jj<L; jj++)
            {
                *(bufPtr) = in[jj];
                bufPtr++;
            }
        }
        for(ltfatInt ii=0; ii<legalExtLen; ii++)
        {
            *(bufPtr) = in[ii];
            bufPtr++;
        }
    }
    break;
    case SP0: // constant padding
        for(ltfatInt ii=0; ii<gl; ii++)
            buf[ii] = in[L-1];
        break;
    case PERDEC: // periodic padding with possible last sample repplication
    {
        ltfatInt rem = L%a;
        if(rem==0)
        {
            for(ltfatInt ii=0; ii<legalExtLen; ii++)
                buf[ii] = in[ii];
        }
        else
        {
            ltfatInt remto = a - rem;
            // replicated
            for(ltfatInt ii=0; ii<remto; ii++)
                buf[ii] = in[L-1];

            // periodized
            for(ltfatInt ii=0; ii<legalExtLen-remto; ii++)
                buf[ii+remto] = in[ii];
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
