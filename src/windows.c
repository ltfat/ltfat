#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN void
LTFAT_NAME(pgauss)(const ltfatInt L, const double w, const double c_t,
                   LTFAT_REAL *g)
{

    ltfatInt lr,k,nk;
    double tmp,sqrtl, safe, gnorm;

    sqrtl=sqrt((double)L);
    safe=4;
    gnorm=0;

    /* Outside the interval [-safe,safe] then exp(-pi*x.^2) is numerically zero. */
    nk=(ltfatInt)ceil(safe/sqrt((double)L/sqrt(w)));

    for ( lr=0; lr<L; lr++)
    {
        g[lr]=0.0;
        for (k=-nk; k<=nk; k++)
        {
            /* Use a tmp variable to calculate squaring */
            tmp = ((double)lr+c_t)/sqrtl-(double)k*sqrtl;
            g[lr]+=exp(-PI*tmp*tmp/w);
        }
        gnorm +=g[lr]*g[lr];
    }

    /* Normalize it exactly. */
    gnorm=sqrt(gnorm);

    for ( lr=0; lr<L; lr++)
    {
        g[lr] /= gnorm;
    }
}


/* does not work correctly. This code does:
%for k=-nk:nk
%  tmp=exp(-pi*((lr+c_t)/sqrtl-k*sqrtl).^2/w)
%  g=g+tmp.*cos(2*pi*c_f*(lr/L-k))+i*tmp.*sin(2*pi*c_f*(lr/L-k));
%end;
*/

LTFAT_EXTERN void
LTFAT_NAME(pgauss_cmplx)(const ltfatInt L, const double w, const double c_t, const double c_f,
                         LTFAT_COMPLEX *g)
{

    ltfatInt lr,k,nk;
    double tmp,sqrtl, safe, gnorm;

    sqrtl=sqrt((double)L);
    safe=4;
    gnorm=0;

    /* Outside the interval [-safe,safe] then exp(-pi*x.^2) is numerically zero. */
    nk=(ltfatInt)ceil(safe/sqrt((double)L/sqrt(w)));

    for ( lr=0; lr<L; lr++)
    {
        //g[lr][0]=0.0;
        //g[lr][1]=0.0;
        g[lr] = (LTFAT_COMPLEX) 0.0;
        for (k=-nk; k<=nk; k++)
        {
            /* Use a tmp variable to calculate squaring */
            tmp = ((double)lr+c_t)/sqrtl-(double)k*sqrtl;
            tmp = exp(-PI*tmp*tmp/w);
            g[lr]+=tmp*LTFAT_COMPLEXH(cexp)(I*2*PI*c_f*((double)lr/L-(double)k));

        }
        double gReal = LTFAT_COMPLEXH(creal)(g[lr]);
        double gImag = LTFAT_COMPLEXH(cimag)(g[lr]);
        gnorm += (gReal*gReal+gImag*gImag);
    }

    /* Normalize it exactly. */
    gnorm=sqrt(gnorm);

    for ( lr=0; lr<L; lr++)
    {
        g[lr] /= gnorm;
    }
}
