#include "ltfat.h"
#include "ltfat_types.h"


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
LTFAT_NAME(filterbankreassign)(const LTFAT_REAL *s[],
                               const LTFAT_REAL *tgrad[],
                               const LTFAT_REAL *fgrad[],
                               const ltfatInt N[], const double a[],
                               const double cfreq[],
                               const ltfatInt M, LTFAT_REAL *sr[])
{
#define CHECKZEROCROSSINGANDBREAK( CMP, SIGN) \
        if ( tmpfgrad CMP 0.0 )\
        {\
           if (abs(tmpfgrad) < abs(oldfgrad))\
           {\
              fgradIdx[jj] = ii;\
           }\
           else\
           {\
              fgradIdx[jj] = ii SIGN 1;\
           }\
           break;\
        }\
        oldfgrad = tmpfgrad;

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
         LTFAT_REAL tmpfgrad;
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

            if (ii == M - 1 && tmpfgrad < 0.0)
            {
               for (ltfatInt ii = 0; ii < m ; ii++)
               {
                  tmpfgrad = cfreq2[ii] - fgradmjj + 2.0;
                  CHECKZEROCROSSINGANDBREAK( >= , -)
               }
               if (fgradIdx[jj] < 0)
               {
                  fgradIdx[jj] = M - 1;
               }
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
            if (ii == 0 && tmpfgrad > 0.0)
            {
               for (ltfatInt ii = M - 1; ii >= m; ii--)
               {
                  tmpfgrad = cfreq2[ii] - fgradmjj - 2.0;
                  CHECKZEROCROSSINGANDBREAK( <= , +)
               }
               if (fgradIdx[jj] >= M)
               {
                  fgradIdx[jj] = 0;
               }
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


      for (ltfatInt ii = 0; ii < N[m]; ii++)
      {
         sr[fgradIdx[ii]][tgradIdx[ii]] += s[m][ii];
      }
   }


   LTFAT_SAFEFREEALL(fgradIdx, tgradIdx, cfreq2);
#undef CHECKZEROCROSSINGANDBREAK
}
