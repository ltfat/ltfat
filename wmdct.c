#include "ltfat.h"
#include "ltfat_types.h"

#define CH(name) LTFAT_COMPLEXH(name)

#define POSTPROC_REAL \
  for (ltfatInt n=0;n<N*W;n+=2) \
  { \
     for (ltfatInt m=0;m<M;m+=2) \
     { \
       pcoef[m]=CH(creal)(pcoef2[m])+CH(cimag)(pcoef2[m]); \
       pcoef[m+M]=CH(creal)(pcoef2[m+M2])-CH(cimag)(pcoef2[m+M2]); \
     } \
     \
     for (ltfatInt m=1;m<M;m+=2) \
     { \
       pcoef[m]=CH(creal)(pcoef2[m])-CH(cimag)(pcoef2[m]); \
       pcoef[m+M]=CH(creal)(pcoef2[m+M2])+CH(cimag)(pcoef2[m+M2]); \
     } \
 \
     pcoef+=M2; \
     pcoef2+=M4; \
  }

#define POSTPROC_COMPLEX \
  for (ltfatInt n=0;n<N*W;n+=2) \
  { \
     for (ltfatInt m=0;m<M;m+=2) \
     { \
         pcoef[m] =   scalconst*(emipi4*pcoef2[m]  +eipi4*pcoef2[M2-1-m]); \
         pcoef[m+M] = scalconst*(eipi4*pcoef2[m+M2]+emipi4*pcoef2[M4-1-m]); \
     } \
 \
     for (ltfatInt m=1;m<M;m+=2) \
     { \
       pcoef[m] = scalconst*(eipi4*pcoef2[m]    +emipi4*pcoef2[M2-1-m]); \
       pcoef[m+M]=scalconst*(emipi4*pcoef2[m+M2]+eipi4*pcoef2[M4-1-m]); \
     } \
 \
     pcoef+=M2; \
     pcoef2+=M4; \
  }

#define PREPROC \
   for(ltfatInt n=0;n<L;n++) \
      f2[n] = (LTFAT_COMPLEX) cexp(-PI*I*n/(2.0*M)); \
   for(ltfatInt w=W-1;w>=0;w--) \
      for(ltfatInt n=0;n<L;n++) \
         f2[n+w*L] = f2[n]*f[n+w*L];


LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(dwiltiii_long)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                                  const ltfatInt L, const ltfatInt W,
                                  const ltfatInt M, LTFAT_COMPLEX *cout)
{
    const ltfatInt N = L / M;
    const ltfatInt M2 = 2 * M;
    const ltfatInt M4 = 4 * M;
    const LTFAT_REAL scalconst = 1.0 / sqrt(2.0);
    const LTFAT_COMPLEX eipi4 = cexp(I * PI / 4.0);
    const LTFAT_COMPLEX emipi4 = cexp(-I * PI / 4.0);

    LTFAT_COMPLEX *coef2 = ltfat_malloc(2 * M * N * W * sizeof * coef2);
    LTFAT_COMPLEX *f2 = ltfat_malloc(L * W * sizeof * f2);

    PREPROC

    LTFAT_NAME_COMPLEX(dgt_long)(f2, g, L, W, M, 2 * M, FREQINV, coef2);

    LTFAT_COMPLEX *pcoef  = cout;
    LTFAT_COMPLEX *pcoef2 = coef2;

    POSTPROC_COMPLEX

    LTFAT_SAFEFREEALL(coef2, f2);
}

LTFAT_EXTERN void
LTFAT_NAME_REAL(dwiltiii_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                               const ltfatInt L, const ltfatInt W,
                               const ltfatInt M, LTFAT_REAL *cout)
{
    const ltfatInt N = L / M;
    const ltfatInt M2 = 2 * M;
    const ltfatInt M4 = 4 * M;

    LTFAT_COMPLEX *coef2 = ltfat_malloc(2 * M * N * W * sizeof * coef2);
    LTFAT_COMPLEX *f2 = ltfat_malloc(L * W * sizeof * f2);
    LTFAT_COMPLEX *g2 = ltfat_malloc(L * sizeof * g2);

    // Real to complex
    for (ltfatInt ii = 0; ii < L; ii++)
        g2[ii] = g[ii];

    PREPROC

    LTFAT_NAME_COMPLEX(dgt_long)(f2, g2, L, W, M, 2 * M, FREQINV, coef2);


    LTFAT_REAL *pcoef  = cout;
    LTFAT_COMPLEX *pcoef2 = coef2;

    POSTPROC_REAL

    LTFAT_SAFEFREEALL(coef2, f2, g2);

}

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(dwiltiii_fb)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                                const ltfatInt L, const ltfatInt gl,
                                const ltfatInt W, const ltfatInt M,
                                LTFAT_COMPLEX *cout)
{
    const ltfatInt N = L / M;
    const ltfatInt M2 = 2 * M;
    const ltfatInt M4 = 4 * M;
    const LTFAT_REAL scalconst = 1.0 / sqrt(2.0);
    const LTFAT_COMPLEX eipi4 = cexp(I * PI / 4.0);
    const LTFAT_COMPLEX emipi4 = cexp(-I * PI / 4.0);

    LTFAT_COMPLEX *coef2 = ltfat_malloc(2 * M * N * W * sizeof * coef2);
    LTFAT_COMPLEX *f2 = ltfat_malloc(L * W * sizeof * f2);

    PREPROC

    /* coef2=comp_dgt(f,g,a,2*M,L); */
    LTFAT_NAME_COMPLEX(dgt_fb)(f2, g, L, gl, W, M, 2 * M, FREQINV, coef2);


    LTFAT_COMPLEX *pcoef  = cout;
    LTFAT_COMPLEX *pcoef2 = coef2;

    POSTPROC_COMPLEX

    LTFAT_SAFEFREEALL(coef2, f2);

}

LTFAT_EXTERN void
LTFAT_NAME_REAL(dwiltiii_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                             const ltfatInt L, const ltfatInt gl,
                             const ltfatInt W, const ltfatInt M,
                             LTFAT_REAL *cout)
{
    const ltfatInt N = L / M;
    const ltfatInt M2 = 2 * M;
    const ltfatInt M4 = 4 * M;

    LTFAT_COMPLEX *coef2 = ltfat_malloc(2 * M * N * W * sizeof * coef2);
    LTFAT_COMPLEX *f2 = ltfat_malloc(L * W * sizeof * f2);
    LTFAT_COMPLEX *g2 = ltfat_malloc(gl * sizeof * g2);

    //Real to complex
    for (ltfatInt ii = 0; ii < gl; ii++)
        g2[ii] = g[ii];

    PREPROC

    LTFAT_NAME_COMPLEX(dgt_fb)(f2, g2, L, gl, W, M, 2 * M, FREQINV, coef2);

    LTFAT_REAL* pcoef  = cout;
    LTFAT_COMPLEX* pcoef2 = coef2;

    POSTPROC_REAL

    LTFAT_SAFEFREEALL(coef2, f2, g2);
}

#undef CH
#undef POSTPROC_REAL
#undef POSTPROC_COMPLEX
#undef PREPROC

