#include "ltfat.h"
#include "ltfat_types.h"

typedef enum { HANN, SQRTHANN, COS, SIN, HAMMING } LTFAT_FIRWIN;

void LTFAT_NAME(firwin)(LTFAT_FIRWIN win, int gl, LTFAT_REAL* g)
{
    double step = 1.0 / gl;
    // for gl even
    double startInt = -0.5;
    double stopInt = 0.5 - 1 / step;

    if (gl % 2 == 1) {
        startInt = -0.5 + step / 2.0;
        stopInt = 0.5 - step / 2.0;
    }

    switch (win) {
    case SQRTHANN:
        for (int ii = 0; ii < gl; ii++) {
            double posInt = startInt + ii * step;
            g[ii] = sqrt(0.5 + 0.5 * cos(2.0 * M_PI * posInt));
        }
        break;
    case HANN:
    case COS:
    case SINE:
        for (int ii = 0; ii < gl; ii++) {
            double posInt = startInt + ii * step;
            g[ii] = 0.5 + 0.5 * cos(2.0 * M_PI * posInt);
        }
        break;
    case HAMMING:
        for (int ii = 0; ii < gl; ii++) {
            double posInt = startInt + ii * step;
            g[ii] = 0.54 + 0.46 * cos(2.0 * M_PI * posInt);
        }
        break;
    };
}

LTFAT_EXTERN void 
LTFAT_NAME(gabframediag)(const LTFAT_TYPE* g, ltfatInt gl,
                         ltfatInt a, ltfatInt L, LTFAT_TYPE* d)
{
    ltfatInt amax = a>L?L:a;

    for (int aIdx = 0; aIdx < amax; aIdx++) 
    {
        for (int ii = gl / 2; ii < gl; ii += a) 
        {
            d[aIdx] += g[ii] * g[ii];
        }
        for (int ii = gl / 2 - a; ii >= 0; ii -= a) 
        {
            d[aIdx] += g[ii] * g[ii];
        }
    }

    if(L>a)
    {
        LTFAT_NAME(periodize_array)(d,a,d,L);
    }
}

LTFAT_EXTERN void
LTFAT_NAME(gabtightfir)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a, ltfatInt M,
                 LTFAT_TYPE* gt)
{

}

LTFAT_EXTERN void
LTFAT_NAME(gabdualfir)(const LTFAT_TYPE* g, ltfatInt gl, ltfatInt a, ltfatInt M, LTFAT_TYPE* gd)
{
    assert(M/((double)a) > 1.0 && "Frame condition");
    assert(M>=gl && "Painless condition");

    // Store the frame diagonal in the first a entries of gd
    LTFAT_NAME(gabframediag)(g,gl,a,a,dg);

    // Invert the diagonal
    for(ltfatInt ii = 0; ii< a;ii++)
        gd[ii] = LTFAT_REAL(1.0)/gd[ii];

    // Multiply window with the inverted diagonal
    // 1) Do the first half, skip the fist a entries 
    for (ltfatInt ii = a, jj = 0; ii < gl/2; ii++, jj++)
        gd[ii] = g[ii] * gd[jj % a];

    // 2) Do the second half in the reverse order
    for (ltfatInt ii = gl - 1, jj = a - 1; ii >= gl/2; ii--, jj--)
        gd[ii] = g[ii] * gd[jj % a];

    // 3) Finalize by doing the first a entries
    for(ltfatInt ii=0; ii<a;ii++)
        gd[ii] *= g[ii];
}
