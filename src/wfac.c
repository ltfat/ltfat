#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(wfac)(const LTFAT_COMPLEX *g, const ltfatInt L, const ltfatInt R,
                         const ltfatInt a, const ltfatInt M,
                         LTFAT_COMPLEX *gf)
{
    ltfatInt h_a, h_m, s;

    LTFAT_REAL *sbuf, *gfp;

    ltfatInt rem, negrem;

    LTFAT_FFTW(plan) p_before;

    const ltfatInt b = L/M;

    const ltfatInt c=gcd(a, M,&h_a, &h_m);
    const ltfatInt p=a/c;
    const ltfatInt q=M/c;
    const ltfatInt d=b/p;

    const double sqrtM=sqrt(M);

    sbuf = (LTFAT_REAL*)ltfat_malloc(2*d*sizeof(LTFAT_REAL));

    /* Create plan. In-place. */
    p_before = LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)sbuf, (LTFAT_COMPLEX*)sbuf,
                                       FFTW_FORWARD, FFTW_MEASURE);

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
                    negrem = positiverem(k*M-l*a,L);
                    for (s=0; s<d; s++)
                    {
                        rem = (negrem+s*p*M)%L;
                        sbuf[2*s]   = LTFAT_COMPLEXH(creal)(sqrtM*g[r+rem+L*w]);
                        sbuf[2*s+1] = LTFAT_COMPLEXH(cimag)(sqrtM*g[r+rem+L*w]);
                    }

                    LTFAT_FFTW(execute)(p_before);

                    for (s=0; s<2*d; s+=2)
                    {
                        gfp[s*ld3]  = sbuf[s];
                        gfp[s*ld3+1]= sbuf[s+1];
                    }
                    gfp+=2;
                }
            }
        }
    }

    ltfat_free(sbuf);
    LTFAT_FFTW(destroy_plan)(p_before);
}


/* wfac for real valued input. */
LTFAT_EXTERN void
LTFAT_NAME_REAL(wfac)(const LTFAT_REAL *g, const ltfatInt L, const ltfatInt R,
                      const ltfatInt a, const ltfatInt M,
                      LTFAT_COMPLEX *gf)
{

    ltfatInt h_a, h_m;

    LTFAT_REAL *sbuf, *gfp;

    ltfatInt s;
    ltfatInt rem, negrem;

    LTFAT_FFTW(plan) p_before;

    const ltfatInt b=L/M;
    const ltfatInt c=gcd(a, M,&h_a, &h_m);
    const ltfatInt p=a/c;
    const ltfatInt q=M/c;
    const ltfatInt d=b/p;

    const double sqrtM=sqrt(M);

    sbuf = (LTFAT_REAL*)ltfat_malloc(2*d*sizeof(LTFAT_REAL));

    /* Create plan. In-place. */
    p_before = LTFAT_FFTW(plan_dft_1d)(d, (LTFAT_COMPLEX*)sbuf,
                                       (LTFAT_COMPLEX*)sbuf,
                                       FFTW_FORWARD, FFTW_MEASURE);

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
                    negrem = positiverem(k*M-l*a,L);
                    for (s=0; s<d; s++)
                    {
                        rem = (negrem+s*p*M)%L;
                        sbuf[2*s]   = sqrtM*g[r+rem+L*w];
                        sbuf[2*s+1] = 0.0;
                    }

                    LTFAT_FFTW(execute)(p_before);

                    for (s=0; s<2*d; s+=2)
                    {
                        gfp[s*ld3]  = sbuf[s];
                        gfp[s*ld3+1]= sbuf[s+1];
                    }
                    gfp+=2;
                }
            }
        }
    }

    ltfat_free(sbuf);
    LTFAT_FFTW(destroy_plan)(p_before);
}

/* wfac for real valued input. Produces only half the output coefficients of wfac_r */
LTFAT_EXTERN void
LTFAT_NAME(wfacreal)(const LTFAT_REAL *g, const ltfatInt L, const ltfatInt R,
                     const ltfatInt a, const ltfatInt M,
                     LTFAT_COMPLEX *gf)
{

    ltfatInt h_a, h_m;

    //LTFAT_REAL *gfp;
    LTFAT_COMPLEX *gfp = gf;

    ltfatInt s;
    ltfatInt rem, negrem;

    LTFAT_FFTW(plan) p_before;

    const ltfatInt b=L/M;
    const ltfatInt c=gcd(a, M,&h_a, &h_m);
    const ltfatInt p=a/c;
    const ltfatInt q=M/c;
    const ltfatInt d=b/p;

    /* This is a floor operation. */
    const ltfatInt d2= d/2+1;

    const double sqrtM=sqrt(M);

    LTFAT_REAL *sbuf = (LTFAT_REAL*)ltfat_malloc(d*sizeof(LTFAT_REAL));
    LTFAT_COMPLEX *cbuf = (LTFAT_COMPLEX*)ltfat_malloc(d2*sizeof(LTFAT_COMPLEX));

    /* Create plan. In-place. */
    p_before = LTFAT_FFTW(plan_dft_r2c_1d)(d, sbuf, cbuf, FFTW_MEASURE);

    // const ltfatInt ld3=2*c*p*q*R;
    const ltfatInt ld3=c*p*q*R;
    //gfp=(LTFAT_REAL*)gf;
    for (ltfatInt r=0; r<c; r++)
    {
        for (ltfatInt w=0; w<R; w++)
        {
            for (ltfatInt l=0; l<q; l++)
            {
                for (ltfatInt k=0; k<p; k++)
                {
                    negrem = positiverem(k*M-l*a,L);
                    for (s=0; s<d; s++)
                    {
                        rem = (negrem+s*p*M)%L;
                        sbuf[s]   = sqrtM*g[r+rem+L*w];
                    }

                    LTFAT_FFTW(execute)(p_before);

                    for (s=0; s<d2; s++)
                    {
                        gfp[s*ld3] = cbuf[s];
                    }
                    gfp++;
                }
            }
        }
    }

    LTFAT_SAFEFREEALL(sbuf,cbuf);
    LTFAT_FFTW(destroy_plan)(p_before);
}
