#include "ltfat.h"
#include "ltfat_types.h"

#define FIRWIN_RESETCOUNTER do{ \
            if (ii == domod.quot + domod.rem) \
                posInt = startInt; \
            }while(0)

#define GABDIAGAPPLY(gg) do{ \
    for (ltfatInt ii = 0; ii < domod.quot + domod.rem; ii++) \
        (gg)[ii] = g[ii] * d[ii % a]; \
\
    for (ltfatInt ii = gl - 1, jj = a - 1; ii >= domod.quot + domod.rem; ii--, jj--)\
    {\
        if(jj<0) jj=a-1;\
        (gg)[ii] = g[ii] * d[jj];\
    }\
}while(0)


LTFAT_EXTERN int
LTFAT_NAME(firwin)(LTFAT_FIRWIN win, int gl, LTFAT_TYPE* g)
{
    double step = 1.0 / gl;
    // for gl even
    double startInt = -0.5;
    const div_t domod = div(gl, 2);

    if (domod.rem)
        startInt = -0.5 + step / 2.0;

    double posInt = 0;

    switch (win)
    {
    case HANN:
    {
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;

            g[ii] = 0.5 + 0.5 * cos(2.0 * M_PI * posInt);
            posInt += step;
        }

        break;
    }

    case SQRTHANN:
    case COS:
    case SINE:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;

            g[ii] = sqrt(0.5 + 0.5 * cos(2.0 * M_PI * posInt));
            posInt += step;
        }

        break;

    case HAMMING:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;

            g[ii] = 0.54 + 0.46 * cos(2.0 * M_PI * posInt);
            posInt += step;
        }

        break;

    default:
        return LTFATERR_UNKNOWNWIN;
    };

    // Fix symmetry of windows which are not zero at -0.5
    if (!domod.rem)
        g[domod.quot + domod.rem] = 0.0;

    return LTFATERR_SUCCESS;
}

// Return first dl entries of the frame diagonal.
LTFAT_EXTERN void
LTFAT_NAME(gabframediag)(const LTFAT_TYPE* g, ltfatInt gl,
                         ltfatInt a, ltfatInt M, ltfatInt dl, LTFAT_TYPE* d)
{
    ltfatInt amax = a > dl ? dl : a;
    // d is used as an accumulator
    memset(d, 0, dl * sizeof * d);

    const div_t domod = div(gl, 2);

    // First half
    for (int aIdx = 0; aIdx < amax; aIdx++)
    {
        for (int ii = aIdx; ii < domod.quot + domod.rem; ii += a)
        {
            LTFAT_REAL gabs = fabs(g[ii]);
            d[aIdx] += gabs * gabs;
        }
    }

    // Second half from the back
    for (int aIdx = amax - 1; aIdx >= 0; aIdx--)
    {
        for (int ii = gl - (a - aIdx); ii >= domod.quot + domod.rem; ii -= a)
        {
            LTFAT_REAL gabs = fabs(g[ii]);
            d[aIdx] += gabs * gabs;
        }
    }

    for (int aIdx = 0; aIdx < amax; aIdx++)
        d[aIdx] *= M;

    // frame diagonal is a-periodic
    if (dl > a)
        LTFAT_NAME(periodize_array)(d, a, d, dl);
}

LTFAT_EXTERN int
LTFAT_NAME(gabtight_painless)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a,
                              ltfatInt M, LTFAT_TYPE* gt)
{
    if (M < a || gl < a)
        return LTFATERR_NOTAFRAME;

    if (M < gl)
        return LTFATERR_NOTPAINLESS;

    LTFAT_TYPE* d = ltfat_malloc(a * sizeof * d);
    LTFAT_NAME(gabframediag)(g, gl, a, M, a, d);

    const div_t domod = div(gl, 2);

    // Invert the diagonal
    for (ltfatInt ii = 0; ii < a; ii++)
        d[ii] = ((LTFAT_REAL) 1.0) / sqrt(d[ii]);

    GABDIAGAPPLY(gt);

    ltfat_free(d);
    return LTFATERR_SUCCESS;
}

LTFAT_EXTERN int
LTFAT_NAME(gabdual_painless)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a,
                             ltfatInt M, LTFAT_TYPE* gd)
{
    if (M < a || gl < a)
        return LTFATERR_NOTAFRAME;

    if (M < gl)
        return LTFATERR_NOTPAINLESS;

    // This temporary array cannot be substituted by gd
    LTFAT_TYPE* d = ltfat_malloc(a * sizeof * d);

    if (!d)
        return LTFATERR_MEMERR;

    LTFAT_NAME(gabframediag)(g, gl, a, M, a, d);

    const div_t domod = div(gl, 2);

    // Invert the diagonal
    for (ltfatInt ii = 0; ii < a; ii++)
        d[ii] = ((LTFAT_REAL) 1.0) / d[ii];

    GABDIAGAPPLY(gd);

    ltfat_free(d);
    return LTFATERR_SUCCESS;
}


#undef FIRWIN_RESETCOUNTER
#undef GABDIAGAPPLY
