#include "ltfat.h"
#include "ltfat_types.h"
#include <float.h>

LTFAT_EXTERN void
LTFAT_NAME(logarray)(const LTFAT_REAL* in, ltfatInt L, LTFAT_REAL* out)
{
#ifdef LTFAT_DOUBLE
    LTFAT_REAL eps = DBL_MIN;
#elif defined(LTFAT_SINGLE)
    LTFAT_REAL eps = FLT_MIN;
#endif
    for (ltfatInt l = 0; l < L; l++)
        out[l] = log(in[l] + eps);

}

LTFAT_EXTERN void
LTFAT_NAME(fbmagphasegrad)(const LTFAT_REAL* logs,
                              const LTFAT_REAL* sqtfr,
                              ltfatInt N[], double a[], double fc[], ltfatInt M,
                              ltfatInt neigh[], double posInfo[], LTFAT_REAL gderivweight,
                              LTFAT_REAL tgrad[], LTFAT_REAL fgrad[])
{
    LTFAT_REAL L = a[0] * N[0];
    ltfatInt chStart = 0;

    for (ltfatInt m = 0; m < M; m++)
    {
        const LTFAT_REAL* logsCol = logs + chStart;
        LTFAT_REAL* fgradCol = fgrad + chStart;

        for (ltfatInt n = 1; n < N[m] - 1; n++)
        {
            fgradCol[n] = (logsCol[n + 1] - logsCol[n - 1]) / 2.0;
        }

        fgradCol[0]        = (logsCol[1] - logsCol[N[m] - 1]) / 2.0;
        fgradCol[N[m] - 1] = (logsCol[0] - logsCol[N[m] - 2]) / 2.0;

        chStart += N[m];
    }


    chStart = 0;
    for (ltfatInt m = 0; m < M; m++)
    {
        LTFAT_REAL* tgradCol = tgrad + chStart;
        LTFAT_REAL* fgradCol = fgrad + chStart;

        LTFAT_REAL aboveNom = 0, aboveDenom = 1, belowNom = 0, belowDenom = 1;
        LTFAT_REAL denom = sqtfr[m] * sqtfr[m] * (M_PI * L);

        if (m < M - 1)
        {
            aboveNom = gderivweight * (sqtfr[m + 1] - sqtfr[m]) / sqtfr[m];
            aboveDenom = fc[m + 1] - fc[m];
        }
        if ( m > 0)
        {
            belowNom = gderivweight * (sqtfr[m] - sqtfr[m - 1]) / sqtfr[m];
            belowDenom = fc[m] - fc[m - 1];
        }

        for (ltfatInt n = 0; n < N[m]; n++)
        {
            ltfatInt w = chStart + n;
            ltfatInt* neighCol = neigh + 6*w;
            LTFAT_REAL tempValAbove = 0;
            LTFAT_REAL tempValBelow = 0;

            int numNeigh = 0;
            for (int jj = 0; jj < 2; jj++)
            {
                ltfatInt oneneigh = neighCol[4 + jj];
                if (oneneigh >= 0)
                {
                    tempValAbove += logs[oneneigh] - logs[w] -
                                    fgrad[w] * (posInfo[oneneigh * 2 + 1] - posInfo[w * 2 + 1]) / a[m];
                    numNeigh++;
                }
            }

            if (numNeigh)
                tempValAbove /= (LTFAT_REAL) numNeigh;


            numNeigh = 0;
            for (int jj = 0; jj < 2; jj++)
            {
                ltfatInt oneneigh = neighCol[2 + jj];
                if (oneneigh >= 0)
                {
                    tempValBelow += logs[w] - logs[oneneigh] -
                                    fgrad[w] * (posInfo[oneneigh * 2 + 1] - posInfo[w * 2 + 1]) / a[m];
                    numNeigh++;
                }
            }

            if (numNeigh)
                tempValBelow /= (LTFAT_REAL) numNeigh;


            tgradCol[n] = (tempValAbove + aboveNom) / aboveDenom +
                          (tempValBelow + belowNom) / belowDenom;
            tgradCol[n] /= denom;
        }
        chStart += N[m];
    }

    // And adjust fgrad
    chStart = 0;
    for (ltfatInt m = 0; m < M; m++)
    {
        LTFAT_REAL* fgradCol = fgrad + chStart;
        LTFAT_REAL fac = sqtfr[m] * sqtfr[m] * N[m] / (2.0 * M_PI);
        for (ltfatInt n = 0; n < N[m]; n++)
            fgradCol[n] *= fac;
        chStart += N[m];
    }
}
