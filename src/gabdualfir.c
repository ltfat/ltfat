#include "ltfat.h"
#include "ltfat_types.h"

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

#define CHECKGABPAINLESS do{ \
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");\
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");\
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");\
    CHECK(LTFATERR_NOTAFRAME, M > a && gl >= a, "Not form a frame. Check if M > a && gl >= a");\
    CHECK(LTFATERR_NOTPAINLESS, M >= gl, "Not painless. Check if M>=gl");\
}while(0)


// Return first dl entries of the frame diagonal.
LTFAT_EXTERN int
LTFAT_NAME(gabframediag)(const LTFAT_TYPE* g, ltfatInt gl,
                         ltfatInt a, ltfatInt M, ltfatInt dl, LTFAT_TYPE* d)
{
    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");
    CHECK(LTFATERR_NOTPOSARG, dl > 0, "d must be positive");

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

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(gabtight_painless)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a,
                              ltfatInt M, LTFAT_TYPE* gt)
{
    LTFAT_TYPE* d = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKGABPAINLESS;

    CHECKMEM(d = ltfat_malloc(a * sizeof * d));

    CHECKSTATUS(LTFAT_NAME(gabframediag)(g, gl, a, M, a, d),
                "Call to gabframediag failed.");

    const div_t domod = div(gl, 2);

    // Invert and square root the diagonal
    for (ltfatInt ii = 0; ii < a; ii++)
        d[ii] = 1.0 / sqrt(d[ii]);

    GABDIAGAPPLY(gt);

error:
    if (d) ltfat_free(d);
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(gabdual_painless)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a,
                             ltfatInt M, LTFAT_TYPE* gd)
{
    LTFAT_TYPE* d = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKGABPAINLESS;

    // This temporary array cannot be substituted by gd
    CHECKMEM( d = ltfat_malloc(a * sizeof * d));

    CHECKSTATUS(LTFAT_NAME(gabframediag)(g, gl, a, M, a, d),
                "Call to gabramediag failed.");

    const div_t domod = div(gl, 2);

    // Invert the diagonal
    for (ltfatInt ii = 0; ii < a; ii++)
        d[ii] = 1.0 / d[ii];

    GABDIAGAPPLY(gd);

error:
    if (d) ltfat_free(d);
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(gabpu_painless)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a,
                           ltfatInt M, LTFAT_TYPE* gpu)
{
    int status = LTFATERR_SUCCESS;
    CHECKSTATUS(LTFAT_NAME(gabdual_painless)(g, gl, a, M, gpu),
                "Call to gabdual_painless failed");

    for (ltfatInt ii = 0; ii < gl; ii++)
        gpu[ii] *= g[ii];

error:
    return status;
}


#undef GABDIAGAPPLY


#define FIRWIN_RESETCOUNTER do{ \
if (ii == domod.quot + domod.rem) \
                posInt = startInt; \
            }while(0)

LTFAT_EXTERN int
LTFAT_NAME(firwin)(LTFAT_FIRWIN win, int gl, LTFAT_TYPE* g)
{
    int status = LTFATERR_SUCCESS;
    CHECK(LTFATERR_NOTPOSARG, gl > 0, "gl must be positive");

    double step = 1.0 / gl;
    // for gl even
    double startInt = -0.5;
    const div_t domod = div(gl, 2);

    if (domod.rem)
        startInt = -0.5 + step / 2.0;

    double posInt = 0;

    switch (win)
    {
    case LTFAT_HANN:
    case LTFAT_HANNING:
    case LTFAT_NUTTALL10:
    {
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 0.5 + 0.5 * cos(2.0 * M_PI * posInt);
            posInt += step;
        }
        break;
    }

    case LTFAT_SQRTHANN:
    case LTFAT_COSINE:
    case LTFAT_SINE:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = sqrt(0.5 + 0.5 * cos(2.0 * M_PI * posInt));
            posInt += step;
        }
        break;

    case LTFAT_HAMMING:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 0.54 + 0.46 * cos(2.0 * M_PI * posInt);
            posInt += step;
        }
        break;

    case LTFAT_NUTTALL01:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 0.53836 + 0.46164 * cos(2 * M_PI * posInt);
            posInt += step;
        }
        break;

    case LTFAT_SQUARE:
    case LTFAT_RECT:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = fabs(posInt) < 0.5 ? 1.0 : 0.0;
            posInt += step;
        }
        break;

    case LTFAT_TRIA:
    case LTFAT_TRIANGULAR:
    case LTFAT_BARTLETT:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 1.0 - 2.0 * fabs(posInt);
            posInt += step;
        }
        break;

    case LTFAT_SQRTTRIA:
        LTFAT_NAME(firwin)(LTFAT_TRIA, gl, g);
        for (int ii = 0; ii < gl; ii++)
            g[ii] = sqrt(g[ii]);

        break;

    case LTFAT_BLACKMAN:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 0.42 + 0.5 * cos(2 * M_PI * posInt) + 0.08 * cos(4 * M_PI * posInt);
            posInt += step;
        }
        break;

    case LTFAT_BLACKMAN2:
    {
        double denomfac = 1.0 / 18608.0;
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 7938.0 + 9240.0 * cos(2.0 * M_PI * posInt) +
                    1430.0 * cos( 4.0 * M_PI * posInt);
            g[ii] *= denomfac;
            posInt += step;
        }
        break;
    }
    case LTFAT_NUTTALL:
    case LTFAT_NUTTALL12:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 0.355768 + 0.487396 * cos(2.0 * M_PI * posInt) +
                    0.144232 * cos(4.0 * M_PI * posInt) + 0.012604 * cos(6.0 * M_PI * posInt);
            posInt += step;
        }
        break;

    case LTFAT_OGG:
    case LTFAT_ITERSINE:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            double innercos = cos(M_PI * posInt);
            g[ii] = sin(M_PI / 2.0 * innercos * innercos);
            posInt += step;
        }
        break;

    case LTFAT_NUTTALL20:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 3.0 + 4.0 * cos(2.0 * M_PI * posInt) + cos(4.0 * M_PI * posInt);
            g[ii] /= 8.0;
            posInt += step;
        }
        break;

    case LTFAT_NUTTALL11:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 0.40897 + 0.5 * cos(2.0 * M_PI * posInt) +
                    0.09103 * cos( 4.0 * M_PI * posInt);
            posInt += step;
        }
        break;

    case LTFAT_NUTTALL02:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 0.4243801 + 0.4973406 * cos(2.0 * M_PI * posInt) +
                    0.0782793 * cos( 4.0 * M_PI * posInt);
            posInt += step;
        }
        break;

    case LTFAT_NUTTALL30:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 10.0 + 15.0 * cos(2.0 * M_PI * posInt) +
                    6.0 * cos( 4.0 * M_PI * posInt) + cos(6.0 * M_PI * posInt);
            g[ii] /= 32.0;
            posInt += step;
        }
        break;

    case LTFAT_NUTTALL21:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 0.338946 + 0.481973 * cos(2.0 * M_PI * posInt) +
                    0.161054 * cos(4.0 * M_PI * posInt) + 0.018027 * cos(6.0 * M_PI * posInt);
            posInt += step;
        }
        break;

    case LTFAT_NUTTALL03:
        for (int ii = 0; ii < gl; ii++)
        {
            FIRWIN_RESETCOUNTER;
            g[ii] = 0.3635819 + 0.4891775 * cos(2.0 * M_PI * posInt) +
                    0.1365995 * cos( 4.0 * M_PI * posInt) + 0.0106411 * cos(6.0 * M_PI * posInt);
            posInt += step;
        }
        break;

    default:
        CHECKCANTHAPPEN("Unknown window");
    };

    // Fix symmetry of windows which are not zero at -0.5
    if (!domod.rem)
        g[domod.quot + domod.rem] = 0.0;

error:
    return status;
}
#undef FIRWIN_RESETCOUNTER
#undef CHECKGABPAINLESS
