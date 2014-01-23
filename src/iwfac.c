#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(iwfac)(const LTFAT_COMPLEX *gf, const ltfatInt L, const ltfatInt R,
                          const ltfatInt a, const ltfatInt M, LTFAT_COMPLEX *g)
{

    ltfatInt h_a, h_m;

    ltfatInt rem, negrem;

    LTFAT_REAL scaling, *sbuf, *gfp;

    LTFAT_FFTW(plan) p_before;

    const ltfatInt b=L/M;
    const ltfatInt c=gcd(a, M,&h_a, &h_m);
    const ltfatInt p=a/c;
    const ltfatInt q=M/c;
    const ltfatInt d=b/p;

    /* division by d is because of the way FFTW normalizes the transform. */
    scaling=1.0/sqrt(M)/d;

    sbuf = ltfat_malloc(2*d*sizeof*sbuf);

    /* Create plan. In-place. */
    p_before = LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)sbuf, (LTFAT_COMPLEX*)sbuf,
                                       FFTW_BACKWARD, FFTW_MEASURE);



    const ltfatInt ld3=c*p*q*R;
    gfp=(LTFAT_REAL*)gf;

    for (ltfatInt r=0; r<c; r++)
    {
        for (ltfatInt w=0; w<R; w++)
        {
            for (ltfatInt l=0; l<q; l++)
            {
                for (ltfatInt k=0; k<p; k++)
                {
                    negrem=positiverem(k*M-l*a,L);
                    for (ltfatInt s=0; s<2*d; s+=2)
                    {
                        sbuf[s]   = gfp[s*ld3]*scaling;
                        sbuf[s+1] = gfp[s*ld3+1]*scaling;
                    }

                    LTFAT_FFTW(execute)(p_before);

                    for (ltfatInt s=0; s<d; s++)
                    {
                        rem = (negrem+s*p*M)%L;
                        LTFAT_REAL* gTmp = (LTFAT_REAL*) &(g[r+rem+L*w]);
                        gTmp[0] = sbuf[2*s];
                        gTmp[1] = sbuf[2*s+1];
                    }
                    gfp+=2;
                }
            }
        }
    }

    /* Clear the work-array. */
    ltfat_free(sbuf);
    LTFAT_FFTW(destroy_plan)(p_before);
}



LTFAT_EXTERN void
LTFAT_NAME_REAL(iwfac)(const LTFAT_COMPLEX *gf, const ltfatInt L, const ltfatInt R,
                       const ltfatInt a, const ltfatInt M, LTFAT_REAL *g)
{

    ltfatInt h_a, h_m;

    ltfatInt rem, negrem;

    LTFAT_REAL scaling, *sbuf, *gfp;

    LTFAT_FFTW(plan) p_before;

    const ltfatInt b=L/M;
    const ltfatInt c=gcd(a, M,&h_a, &h_m);
    const ltfatInt p=a/c;
    const ltfatInt q=M/c;
    const ltfatInt d=b/p;

    /* division by d is because of the way FFTW normalizes the transform. */
    scaling=1.0/sqrt(M)/d;

    sbuf = ltfat_malloc(2*d*sizeof*sbuf);

    /* Create plan. In-place. */
    p_before = LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)sbuf, (LTFAT_COMPLEX*)sbuf,
                                       FFTW_BACKWARD, FFTW_MEASURE);



    const ltfatInt ld3=c*p*q*R;
    gfp=(LTFAT_REAL*)gf;

    for (ltfatInt r=0; r<c; r++)
    {
        for (ltfatInt w=0; w<R; w++)
        {
            for (ltfatInt l=0; l<q; l++)
            {
                for (ltfatInt k=0; k<p; k++)
                {
                    negrem=positiverem(k*M-l*a,L);
                    for (ltfatInt s=0; s<2*d; s+=2)
                    {
                        sbuf[s]   = gfp[s*ld3]*scaling;
                        sbuf[s+1] = gfp[s*ld3+1]*scaling;
                    }

                    LTFAT_FFTW(execute)(p_before);

                    for (ltfatInt s=0; s<d; s++)
                    {
                        rem = (negrem+s*p*M)%L;
                        g[r+rem+L*w] = sbuf[2*s];
                    }
                    gfp+=2;
                }
            }
        }
    }

    /* Clear the work-array. */
    ltfat_free(sbuf);
    LTFAT_FFTW(destroy_plan)(p_before);
}


LTFAT_EXTERN void
LTFAT_NAME(iwfacreal)(const LTFAT_COMPLEX *gf, const ltfatInt L, const ltfatInt R,
                      const ltfatInt a, const ltfatInt M, LTFAT_REAL *g)
{

    ltfatInt h_a, h_m;

    LTFAT_FFTW(plan) p_before;

    const ltfatInt b=L/M;
    const ltfatInt c=gcd(a, M,&h_a, &h_m);
    const ltfatInt p=a/c;
    const ltfatInt q=M/c;
    const ltfatInt d=b/p;

    /* This is a floor operation. */
    const ltfatInt d2= d/2+1;

    /* division by d is because of the way FFTW normalizes the transform. */
    const LTFAT_REAL scaling=1.0/sqrt(M)/d;

    LTFAT_REAL    *sbuf = ltfat_malloc( d*sizeof*sbuf);
    LTFAT_COMPLEX *cbuf = ltfat_malloc(d2*sizeof*cbuf);

    /* Create plan. In-place. */
    p_before = LTFAT_FFTW(plan_dft_c2r_1d)(d, cbuf, sbuf, FFTW_MEASURE);

    const ltfatInt ld3=c*p*q*R;

    /* Advancing pointer: Runs through array pointing out the base for the strided operations. */
    const LTFAT_COMPLEX *gfp = gf;

    for (ltfatInt r=0; r<c; r++)
    {
        for (ltfatInt w=0; w<R; w++)
        {
            for (ltfatInt l=0; l<q; l++)
            {
                for (ltfatInt k=0; k<p; k++)
                {
                    const ltfatInt negrem=positiverem(k*M-l*a,L);
                    for (ltfatInt s=0; s<d2; s++)
                    {
                        cbuf[s] = gfp[s*ld3]*scaling;
                    }

                    LTFAT_FFTW(execute)(p_before);

                    for (ltfatInt s=0; s<d; s++)
                    {
                        g[r+(negrem+s*p*M)%L+L*w] = sbuf[s];
                    }
                    gfp++;
                }
            }
        }
    }

    /* Clear the work-arrays. */
    LTFAT_SAFEFREEALL(cbuf,sbuf);

    LTFAT_FFTW(destroy_plan)(p_before);
}



