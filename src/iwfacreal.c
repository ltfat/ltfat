#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

LTFAT_API void
LTFAT_NAME(iwfacreal)(const LTFAT_COMPLEX* gf, const ltfatInt L,
                      const ltfatInt R,
                      const ltfatInt a, const ltfatInt M, LTFAT_REAL* g)
{

    ltfatInt h_a, h_m;

    LTFAT_FFTW(plan) p_before;

    const ltfatInt b = L / M;
    const ltfatInt c = ltfat_gcd(a, M, &h_a, &h_m);
    const ltfatInt p = a / c;
    const ltfatInt q = M / c;
    const ltfatInt d = b / p;

    /* This is a floor operation. */
    const ltfatInt d2 = d / 2 + 1;

    /* division by d is because of the way FFTW normalizes the transform. */
    const LTFAT_REAL scaling = (const LTFAT_REAL) ( 1.0 / sqrt((double)M) / d );

    LTFAT_REAL*    sbuf = LTFAT_NAME_REAL(malloc)( d);
    LTFAT_COMPLEX* cbuf = LTFAT_NAME_COMPLEX(malloc)( d2);

    /* Create plan. In-place. */
    p_before = LTFAT_FFTW(plan_dft_c2r_1d)((int)d, (LTFAT_FFTW(complex)*) cbuf, sbuf,
                                           FFTW_MEASURE);

    const ltfatInt ld3 = c * p * q * R;

    /* Advancing pointer: Runs through array pointing out the base for the strided operations. */
    const LTFAT_COMPLEX* gfp = gf;

    for (ltfatInt r = 0; r < c; r++)
    {
        for (ltfatInt w = 0; w < R; w++)
        {
            for (ltfatInt l = 0; l < q; l++)
            {
                for (ltfatInt k = 0; k < p; k++)
                {
                    const ltfatInt negrem = ltfat_positiverem(k * M - l * a, L);
                    for (ltfatInt s = 0; s < d2; s++)
                    {
                        cbuf[s] = gfp[s * ld3] * scaling;
                    }

                    LTFAT_FFTW(execute)(p_before);

                    for (ltfatInt s = 0; s < d; s++)
                    {
                        g[r + (negrem + s * p * M) % L + L * w] = sbuf[s];
                    }
                    gfp++;
                }
            }
        }
    }

    /* Clear the work-arrays. */
    LTFAT_SAFEFREEALL(cbuf, sbuf);

    LTFAT_FFTW(destroy_plan)(p_before);
}


