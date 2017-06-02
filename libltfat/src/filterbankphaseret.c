#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN void
PHASERET_NAME(logarray)(const LTFAT_REAL* in, ltfatInt L, LTFAT_REAL* out)
{
#ifdef LTFAT_DOUBLE
    LTFAT_REAL eps = DBL_MIN;
#elif defined(LTFAT_SINGLE)
    LTFAT_REAL eps = FLT_MIN;
#endif
    for (ltfatInt l = 0; l < L; l++)
        out[l] = log(in[l] + eps);

}

void
PHASERET_NAME(pghitgrad)(const LTFAT_REAL* logs, double gamma, ltfatInt a,
                         ltfatInt M,
                         ltfatInt N,
                         LTFAT_REAL* tgrad)
{
    ltfatInt M2 = M / 2 + 1;

    const LTFAT_REAL tgradmul = (a * M) / (gamma * 2.0);
    const LTFAT_REAL tgradplus = 2.0 * M_PI * a / M;


    for (ltfatInt n = 0; n < N; n++)
    {
        LTFAT_REAL* tgradCol = tgrad + n * M2;
        const LTFAT_REAL* logsCol = logs + n * M2;

        tgradCol[0]      = 0.0;
        tgradCol[M2 - 1] = 0.0;

        for (ltfatInt m = 1; m < M2 - 1; m++)
            tgradCol[m] = tgradmul * (logsCol[m + 1] - logsCol[m - 1]) + tgradplus * m;
    }
}

void
PHASERET_NAME(filterbankphasegrad)(const LTFAT_REAL* logs, const LTFAT_REAL* sqtfr,
                                   ltfatInt N[], double a[], double fc[], ltfatInt M,
                                   ltfatInt neigh[], double posInfo, LTFAT_REAL gderivweight,
                                   LTFAT_REAL tgrad[], LTFAT_REAL fgrad[])
{
    ltfatInt chStart = 0;

    for (ltfatInt m = 0; m < M; m++)
    {
        const LTFAT_REAL* logsCol = logs + chStart;
        LTFAT_REAL* fgradCol = fgradCol + chStart;

        for (ltfatInt n = 1; n < N[m] - 1; n++)
        {
            fgradCol[n] = (logsCol[m + 1] - logsCol[m - 1])/2.0;
        }

        fgradCol[0]        = (logsCol[1] - logsCol[N[m] - 1])/2.0;
        fgradCol[N[m] - 1] = (logsCol[0] - logsCol[N[m] - 2])/2.0;

        chStart += N[m];
    }






    chStart = 0;
    for (ltfatInt m = 0; m < M; m++)
    {
        LTFAT_REAL* fgradCol = fgradCol + chStart;

        for (ltfatInt n = 1; n < N[m] - 1; n++)
            fgradCol[n] *= sqtfr*sqtfr*N[m]/(2.0*pi);
        
        chStart += N[m];
    }
}


