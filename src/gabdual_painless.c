#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

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
    CHECK(LTFATERR_NOTAFRAME, M > a && gl >= a, "Not a frame. Check if M > a && gl >= a");\
    CHECK(LTFATERR_NOTPAINLESS, M >= gl, "Not painless. Check if M>=gl");\
}while(0)


// Return first dl entries of the frame diagonal.
LTFAT_EXTERN int
LTFAT_NAME(gabframediag)(const LTFAT_TYPE* g, ltfatInt gl,
                         ltfatInt a, ltfatInt M, ltfatInt dl, LTFAT_TYPE* d)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(g); CHECKNULL(d);
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
        LTFAT_NAME(periodize_array)(d, a, dl, d);

error:
    return status;
}

LTFAT_EXTERN int
LTFAT_NAME(gabtight_painless)(const LTFAT_TYPE* g, const ltfatInt gl,
                              const ltfatInt a, const ltfatInt M, LTFAT_TYPE* gt)
{
    LTFAT_TYPE* d = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(g); CHECKNULL(gt);
    CHECKGABPAINLESS;

    CHECKMEM( d = ltfat_malloc(a * sizeof * d));

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
LTFAT_NAME(gabdual_painless)(const LTFAT_TYPE* g, const ltfatInt gl,
                             const ltfatInt a,  const ltfatInt M, LTFAT_TYPE* gd)
{
    LTFAT_TYPE* d = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(g); CHECKNULL(gd);
    CHECKGABPAINLESS;

    // This temporary array cannot be substituted by gd
    CHECKMEM( d = ltfat_malloc(a * sizeof * d));

    CHECKSTATUS(LTFAT_NAME(gabframediag)(g, gl, a, M, a, d),
                "gabramediag failed.");

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
#undef CHECKGABPAINLESS
