#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

LTFAT_EXTERN int
LTFAT_NAME(pgauss)(const ltfatInt L, const double w, const double c_t,
                   LTFAT_REAL* g)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(g);
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive (passed %d).", L);
    CHECK(LTFATERR_NOTPOSARG, w > 0, "w must be positive (passed %f).", w);

    ltfatInt lr, k, nk;
    double tmp, sqrtL, safe, gnorm;

    sqrtL = sqrt((double)L);
    safe = 4;
    gnorm = 0;

    /* Outside the interval [-safe,safe] then exp(-pi*x.^2) is numerically zero. */
    nk = (ltfatInt)ceil(safe / sqrt((double)L / sqrt(w)));

    for ( lr = 0; lr < L; lr++)
    {
        g[lr] = 0.0;
        for (k = -nk; k <= nk; k++)
        {
            /* Use a tmp variable to calculate squaring */
            tmp = ((double)lr + c_t) / sqrtL - (double)k * sqrtL;
            g[lr] += exp(-M_PI * tmp * tmp / w);
        }
        gnorm += g[lr] * g[lr];
    }

    /* Normalize it exactly. */
    gnorm = sqrt(gnorm);

    for ( lr = 0; lr < L; lr++)
        g[lr] /= gnorm;

error:
    return status;
}


/* does not work correctly. This code does:
%for k=-nk:nk
%  tmp=exp(-pi*((lr+c_t)/sqrtL-k*sqrtL).^2/w)
%  g=g+tmp.*cos(2*pi*c_f*(lr/L-k))+i*tmp.*sin(2*pi*c_f*(lr/L-k));
%end;
*/

LTFAT_EXTERN int
LTFAT_NAME(pgauss_cmplx)(const ltfatInt L, const double w, const double c_t,
                         const double c_f,
                         LTFAT_COMPLEX* g)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(g);
    CHECK(LTFATERR_BADSIZE, L > 0, "L must be positive (passed %d).", L);
    CHECK(LTFATERR_NOTPOSARG, w > 0, "w must be positive (passed %f).", w);

    ltfatInt lr, k, nk;
    double tmp, sqrtL, safe, gnorm;

    sqrtL = sqrt((double)L);
    safe = 4;
    gnorm = 0;

    /* Outside the interval [-safe,safe] then exp(-pi*x.^2) is numerically zero. */
    nk = (ltfatInt)ceil(safe / sqrt((double)L / sqrt(w)));

    for ( lr = 0; lr < L; lr++)
    {
        g[lr] = (LTFAT_COMPLEX) 0.0;
        for (k = -nk; k <= nk; k++)
        {
            /* Use a tmp variable to calculate squaring */
            tmp = ((double)lr + c_t) / sqrtL - (double)k * sqrtL;
            tmp = exp( -M_PI * tmp * tmp / w );
            g[lr] += (LTFAT_REAL)(tmp) *
                     exp(I * (LTFAT_REAL)( 2.0 * M_PI * c_f * ((( double)lr) / ((double)L) - ((
                                               double)k))));

        }
        double gReal = ltfat_real(g[lr]);
        double gImag = ltfat_imag(g[lr]);
        gnorm += (gReal * gReal + gImag * gImag);
    }

    /* Normalize it exactly. */
    gnorm = sqrt(gnorm);

    for ( lr = 0; lr < L; lr++)
        g[lr] /= gnorm;

error:
    return status;
}
