#include "ltfat.h"
#include "ltfat_types.h"
#include "reassign_typeconstant.h"


LTFAT_EXTERN void
LTFAT_NAME(gabreassign)(const LTFAT_REAL *s, const LTFAT_REAL *tgrad,
                        const LTFAT_REAL *fgrad, const ltfatInt L, const ltfatInt W,
                        const ltfatInt a, const ltfatInt M, LTFAT_REAL *sr)
{

   ltfatInt ii, posi, posj;


   const ltfatInt N = L / a;
   const ltfatInt b = L / M;

   ltfatInt *timepos = ltfat_malloc(N * sizeof * timepos);
   ltfatInt *freqpos = ltfat_malloc(M * sizeof * freqpos);

   fftindex(N, timepos);
   fftindex(M, freqpos);

   /* Zero the output array. */
   memset(sr, 0, M * N * W * sizeof * sr);

   for (ltfatInt w = 0; w < W; w++)
   {
      for (ii = 0; ii < M; ii++)
      {
         for (ltfatInt jj = 0; jj < N; jj++)
         {
            /* Do a 'round' followed by a 'mod'. 'round' is not
             * present in all libraries, so use trunc(x+.5) instead */
            /*posi=positiverem((ltfatInt)trunc(tgrad[ii+jj*M]/b+freqpos[ii]+.5),M);
              posj=positiverem((ltfatInt)trunc(fgrad[ii+jj*M]/a+timepos[jj]+.5),N);*/
            posi = positiverem(ltfat_round(tgrad[ii + jj * M] / b + freqpos[ii]), M);
            posj = positiverem(ltfat_round(fgrad[ii + jj * M] / a + timepos[jj]), N);

            sr[posi + posj * M] += s[ii + jj * M];
         }
      }
   }

   LTFAT_SAFEFREEALL(freqpos, timepos);
}

LTFAT_EXTERN void
LTFAT_NAME(filterbankphasegrad)(const LTFAT_COMPLEX* c [],
                                const LTFAT_COMPLEX* ch[],
                                const LTFAT_COMPLEX* cd[],
                                const ltfatInt          M,
                                const ltfatInt        N[],
                                const ltfatInt          L,
                                const LTFAT_REAL   minlvl,
                                LTFAT_REAL*        tgrad[],
                                LTFAT_REAL*        fgrad[],
                                LTFAT_REAL*           cs[])
{
#define FOREACHCOEF \
    for(ltfatInt m=0;m<M;++m){\
        for(ltfatInt ii=0;ii<N[m];++ii){

#define ARRAYEL(c) ((c)[m][ii])
#define ENDFOREACHCOEF }}

LTFAT_REAL minlvlAlt = LTFAT_COMPLEXH(cabs)(c[0][0]);

// Compute spectrogram from coefficients
// Keep max value
FOREACHCOEF
LTFAT_REAL en = LTFAT_COMPLEXH(cabs)(ARRAYEL(c))*LTFAT_COMPLEXH(cabs)(ARRAYEL(c));
ARRAYEL(cs) = en;
if(en>minlvlAlt)
    minlvlAlt = en;
ENDFOREACHCOEF

// Adjust minlvl 
minlvlAlt *= minlvl;

// Force spectrogram values less tha minLvlAlt to minlvlAlt
FOREACHCOEF
LTFAT_REAL csEl = ARRAYEL(cs);
if(csEl<minlvlAlt)
    ARRAYEL(cs) = minlvlAlt;
ENDFOREACHCOEF

// Instantaneous frequency
FOREACHCOEF
LTFAT_REAL tgradEl = LTFAT_COMPLEXH(creal)(
                         ARRAYEL(cd)*LTFAT_COMPLEXH(conj)(ARRAYEL(c))/ARRAYEL(cs)
                                          )/L*2;
ARRAYEL(tgrad) = fabs(tgradEl)<=2?tgradEl:0.0f;
ENDFOREACHCOEF


FOREACHCOEF
ARRAYEL(fgrad) = LTFAT_COMPLEXH(cimag)(
                        ARRAYEL(ch)*LTFAT_COMPLEXH(conj)(ARRAYEL(c))/ARRAYEL(cs)
                                      );
ENDFOREACHCOEF

#undef FOREACHCOEF
#undef ENDFOREACHCOEF
#undef ARRAYEL
}

LTFAT_EXTERN void
LTFAT_NAME(filterbankreassign)(const LTFAT_REAL *s[],
                               const LTFAT_REAL *tgrad[],
                               const LTFAT_REAL *fgrad[],
                               const ltfatInt N[], const double a[],
                               const double cfreq[], const ltfatInt M,
                               LTFAT_REAL *sr[],
                               fbreassHints hints,
                               fbreassOptOut  *repos)
{
#define CHECKZEROCROSSINGANDBREAK( CMP, SIGN) \
     { \
        if ( (tmptgrad) CMP 0.0 )\
        {\
           if (fabs(tmptgrad) < fabs(oldtgrad))\
           {\
              tgradIdx[jj] = ii;\
           }\
           else\
           {\
              tgradIdx[jj] = ii SIGN 1;\
           }\
           break;\
        }\
        oldtgrad = tmptgrad;\
     }

   ltfatInt* chan_pos = NULL;

   int doTimeWraparound = !(hints & REASS_NOTIMEWRAPAROUND);

   if (repos)
   {
      chan_pos = ltfat_malloc((M + 1) * sizeof * chan_pos);

      chan_pos[0] = 0.0;
      for (ltfatInt ii = 0; ii < M; ii++)
      {
         chan_pos[ii + 1] = chan_pos[ii] + N[ii];
      }
   }

   /* Limit tgrad? */

   LTFAT_REAL oneover2 = 1.0 / 2.0;

   // This will hold center frequencies modulo 2.0
   double *cfreq2 = ltfat_malloc(M * sizeof * cfreq2);

   for (ltfatInt m = 0; m < M; m++)
   {
      // Zero the output arrays
      memset(sr[m], 0, N[m]*sizeof(LTFAT_REAL));
      // This is effectivelly modulo by 2.0
      cfreq2[m] = cfreq[m] - floor(cfreq[m] * oneover2) * 2.0;
   }

   ltfatInt* tgradIdx = NULL;
   ltfatInt* fgradIdx = NULL;
   ltfatInt Nold = 0;
   for (ltfatInt m = M - 1; m >= 0; m--)
   {
      // Ensure the temporary arrays have proper lengths
      if (N[m] > Nold)
      {
         if (tgradIdx)
         {
            ltfat_free(tgradIdx);
         }
         if (fgradIdx)
         {
            ltfat_free(fgradIdx);
         }

         tgradIdx = ltfat_malloc(N[m] * sizeof * tgradIdx);
         fgradIdx = ltfat_malloc(N[m] * sizeof * fgradIdx);
         Nold = N[m];
      }

      // We will use this repeatedly
      LTFAT_REAL cfreqm = cfreq2[m];

      /************************
       *
       * Calculating frequency reassignment
       *
       * **********************
       */
      for (ltfatInt jj = 0; jj < N[m]; jj++)
      {
         //
         LTFAT_REAL tmptgrad = 0.0;
         LTFAT_REAL tgradmjj = tgrad[m][jj] + cfreqm;
         LTFAT_REAL oldtgrad = 10; // 10 seems to be big enough
         // Zero this in case it falls trough, although it might not happen
         tgradIdx[jj] = 0;

         if (tgrad[m][jj] > 0)
         {
            ltfatInt ii;
            // Search for zero crossing

            // If the gradient is bigger than 0, start from m upward....
            for (ii = m; ii < M; ii++)
            {
               tmptgrad = cfreq2[ii] - tgradmjj;
               CHECKZEROCROSSINGANDBREAK( >= , -)
            }
            // If the previous for does not break, ii == M
            if (ii == M  && tmptgrad < 0.0)
            {
               for (ltfatInt ii = 0; ii < m ; ii++)
               {
                  tmptgrad = cfreq2[ii] - tgradmjj + 2.0;
                  CHECKZEROCROSSINGANDBREAK( >= , -)
               }
            }
            if (tgradIdx[jj] < 0)
            {
               tgradIdx[jj] = M - 1;
            }
         }
         else
         {
            ltfatInt ii;
            for (ii = m; ii >= 0; ii--)
            {
               tmptgrad = cfreq2[ii] - tgradmjj;
               CHECKZEROCROSSINGANDBREAK( <= , +)
            }
            // If the previous for does not break, ii=-1
            if (ii == -1 && tmptgrad > 0.0)
            {
               for (ltfatInt ii = M - 1; ii >= m; ii--)
               {
                  tmptgrad = cfreq2[ii] - tgradmjj - 2.0;
                  CHECKZEROCROSSINGANDBREAK( <= , +)
               }
            }
            if (tgradIdx[jj] >= M)
            {
               tgradIdx[jj] = 0;
            }
         }
      }

      /**********************************
       *                                *
       * Calculating time-reassignment  *
       *                                *
       **********************************/

      for (ltfatInt jj = 0; jj < N[m]; jj++)
      {
         ltfatInt tmpIdx = tgradIdx[jj];
         ltfatInt fgradIdxTmp = ltfat_round( (fgrad[m][jj] + a[m] * jj) / a[tmpIdx]);

         if(doTimeWraparound)
         {
            fgradIdx[jj] = positiverem( fgradIdxTmp, N[tmpIdx]);
         }
         else
         {
            fgradIdx[jj] = rangelimit( fgradIdxTmp, 0, N[tmpIdx]-1);
         }
      }


      for (ltfatInt jj = 0; jj < N[m]; jj++)
      {
         sr[tgradIdx[jj]][fgradIdx[jj]] += s[m][jj];
      }

      if (repos && chan_pos)
      {
         for (ltfatInt jj = 0; jj < N[m]; jj++)
         {
            ltfatInt tmpIdx =  chan_pos[tgradIdx[jj]] + fgradIdx[jj] ;
            ltfatInt* tmpl = &repos->reposl[tmpIdx];
            repos->repos[tmpIdx][*tmpl] = chan_pos[m] + jj;
            (*tmpl)++;
            if (*tmpl >= repos->reposlmax[tmpIdx])
            {
               fbreassOptOut_expand(repos, tmpIdx);
            }
         }
      }

   }


   LTFAT_SAFEFREEALL(tgradIdx, fgradIdx, cfreq2, chan_pos);
#undef CHECKZEROCROSSINGANDBREAK
}



