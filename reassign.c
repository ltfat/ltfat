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
LTFAT_REAL fgradEl = LTFAT_COMPLEXH(creal)(
                         ARRAYEL(ch)*LTFAT_COMPLEXH(conj)(ARRAYEL(c))/ARRAYEL(cs)
                                          )/L*2;
ARRAYEL(fgrad) = fabs(fgradEl)<=2?fgradEl:0.0f;
ENDFOREACHCOEF


FOREACHCOEF
ARRAYEL(tgrad) = LTFAT_COMPLEXH(cimag)(
                        ARRAYEL(cd)*LTFAT_COMPLEXH(conj)(ARRAYEL(c))/ARRAYEL(cs)
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
                               fbreassOptOut  *repos)
{
#define CHECKZEROCROSSINGANDBREAK( CMP, SIGN) \
     { \
        if ( (tmpfgrad) CMP 0.0 )\
        {\
           if (fabs(tmpfgrad) < fabs(oldfgrad))\
           {\
              fgradIdx[jj] = ii;\
           }\
           else\
           {\
              fgradIdx[jj] = ii SIGN 1;\
           }\
           break;\
        }\
        oldfgrad = tmpfgrad;\
     }

   ltfatInt* chan_pos = NULL;

   if (repos)
   {
      chan_pos = ltfat_malloc((M + 1) * sizeof * chan_pos);

      chan_pos[0] = 0.0;
      for (ltfatInt ii = 0; ii < M; ii++)
      {
         chan_pos[ii + 1] = chan_pos[ii] + N[ii];
      }
   }

   /* Limit fgrad? */

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

   ltfatInt* fgradIdx = NULL;
   ltfatInt* tgradIdx = NULL;
   ltfatInt Nold = 0;
   for (ltfatInt m = M - 1; m >= 0; m--)
   {
      // Ensure the temporary arrays have proper lengths
      if (N[m] > Nold)
      {
         if (fgradIdx)
         {
            ltfat_free(fgradIdx);
         }
         if (tgradIdx)
         {
            ltfat_free(tgradIdx);
         }

         fgradIdx = ltfat_malloc(N[m] * sizeof * fgradIdx);
         tgradIdx = ltfat_malloc(N[m] * sizeof * tgradIdx);
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
         LTFAT_REAL tmpfgrad = 0.0;
         LTFAT_REAL fgradmjj = fgrad[m][jj] + cfreqm;
         LTFAT_REAL oldfgrad = 10; // 10 seems to be big enough
         // Zero this in case it falls trough, although it might not happen
         fgradIdx[jj] = 0;

         if (fgrad[m][jj] > 0)
         {
            ltfatInt ii;
            // Search for zero crossing

            // If the gradient is bigger than 0, start from m upward....
            for (ii = m; ii < M; ii++)
            {
               tmpfgrad = cfreq2[ii] - fgradmjj;
               CHECKZEROCROSSINGANDBREAK( >= , -)
            }
            // If the previous for does not break, ii == M
            if (ii == M  && tmpfgrad < 0.0)
            {
               for (ltfatInt ii = 0; ii < m ; ii++)
               {
                  tmpfgrad = cfreq2[ii] - fgradmjj + 2.0;
                  CHECKZEROCROSSINGANDBREAK( >= , -)
               }
            }
            if (fgradIdx[jj] < 0)
            {
               fgradIdx[jj] = M - 1;
            }
         }
         else
         {
            ltfatInt ii;
            for (ii = m; ii >= 0; ii--)
            {
               tmpfgrad = cfreq2[ii] - fgradmjj;
               CHECKZEROCROSSINGANDBREAK( <= , +)
            }
            // If the previous for does not break, ii=-1
            if (ii == -1 && tmpfgrad > 0.0)
            {
               for (ltfatInt ii = M - 1; ii >= m; ii--)
               {
                  tmpfgrad = cfreq2[ii] - fgradmjj - 2.0;
                  CHECKZEROCROSSINGANDBREAK( <= , +)
               }
            }
            if (fgradIdx[jj] >= M)
            {
               fgradIdx[jj] = 0;
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
         ltfatInt tmpIdx = fgradIdx[jj];
         tgradIdx[jj] = positiverem(
                           ltfat_round( (tgrad[m][jj] + a[m] * jj) / a[tmpIdx]), N[tmpIdx]);
      }


      for (ltfatInt jj = 0; jj < N[m]; jj++)
      {
         sr[fgradIdx[jj]][tgradIdx[jj]] += s[m][jj];
      }

      if (repos && chan_pos)
      {
         for (ltfatInt jj = 0; jj < N[m]; jj++)
         {
            ltfatInt tmpIdx =  chan_pos[fgradIdx[jj]] + tgradIdx[jj] ;
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


   LTFAT_SAFEFREEALL(fgradIdx, tgradIdx, cfreq2, chan_pos);
#undef CHECKZEROCROSSINGANDBREAK
}



