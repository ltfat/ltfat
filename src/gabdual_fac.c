#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"
#include "ltfat/blaslapack.h"

LTFAT_EXTERN void
LTFAT_NAME(gabdual_fac)(const LTFAT_COMPLEX* gf, const ltfatInt L,
                        const ltfatInt R,
                        const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX* gdualf)
{

    ltfatInt h_a, h_m;

    LTFAT_COMPLEX* Sf;

    const LTFAT_COMPLEX zzero = (LTFAT_COMPLEX) 0.0;//{0.0, 0.0 };
    const LTFAT_COMPLEX alpha = (LTFAT_COMPLEX) 1.0; //{1.0, 0.0 };

    const ltfatInt N = L / a;

    const ltfatInt c = ltfat_gcd(a, M, &h_a, &h_m);
    const ltfatInt p = a / c;
    const ltfatInt q = M / c;
    const ltfatInt d = N / q;

    Sf = LTFAT_NAME_COMPLEX(malloc)(p * p);

    /* Copy the contents of gf to gdualf because LAPACK overwrites it input
     * argument
     */
    memcpy(gdualf, gf, L * R * sizeof * gdualf);

    for (ltfatInt rs = 0; rs < c * d; rs++)
    {
        LTFAT_NAME(ltfat_gemm)(CblasNoTrans, CblasConjTrans, p, p, q * R,
                               &alpha,
                               gf + rs * p * q * R, p,
                               gf + rs * p * q * R, p,
                               &zzero, Sf, p);

        LTFAT_NAME(ltfat_posv)(p, q * R, Sf, p,
                               gdualf + rs * p * q * R, p);

    }

    /* Clear the work-array. */
    ltfat_free(Sf);


}


LTFAT_EXTERN void
LTFAT_NAME(gabdualreal_fac)(const LTFAT_COMPLEX* gf, const ltfatInt L,
                            const ltfatInt R,
                            const ltfatInt a, const ltfatInt M,
                            LTFAT_COMPLEX* gdualf)
{

    ltfatInt h_a, h_m;

    LTFAT_COMPLEX* Sf;

    const LTFAT_COMPLEX zzero = (LTFAT_COMPLEX) 0.0;
    const LTFAT_COMPLEX alpha = (LTFAT_COMPLEX) 1.0; //{1.0, 0.0 };

    const ltfatInt N = L / a;

    const ltfatInt c = ltfat_gcd(a, M, &h_a, &h_m);
    const ltfatInt p = a / c;
    const ltfatInt q = M / c;
    const ltfatInt d = N / q;

    /* This is a floor operation. */
    const ltfatInt d2 = d / 2 + 1;

    Sf = LTFAT_NAME_COMPLEX(malloc)(p * p);

    /* Copy the contents of gf to gdualf because LAPACK overwrites it input
     * argument
     */
    memcpy(gdualf, gf, sizeof(LTFAT_COMPLEX)*L * R);

    for (ltfatInt rs = 0; rs < c * d2; rs++)
    {
        LTFAT_NAME(ltfat_gemm)(CblasNoTrans, CblasConjTrans, p, p, q * R,
                               &alpha,
                               gf + rs * p * q * R, p,
                               gf + rs * p * q * R, p,
                               &zzero, Sf, p);

        LTFAT_NAME(ltfat_posv)(p, q * R, Sf, p,
                               gdualf + rs * p * q * R, p);

    }

    /* Clear the work-array. */
    ltfat_free(Sf);


}
