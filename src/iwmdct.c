#include "ltfat.h"
#include "ltfat_types.h"

#define CH(name) LTFAT_COMPLEXH(name)

#define PREPROC_COMPLEX \
  for (ltfatInt n=0;n<N*W;n+=2) \
  { \
     for (ltfatInt m=0;m<M;m+=2) \
     { \
        pcoef2[m] = eipi4*pcoef[m]; \
        pcoef2[M2-1-m] = emipi4*pcoef[m]; \
        pcoef2[m+M2] = emipi4*pcoef[m+M]; \
        pcoef2[M4-1-m] = eipi4*pcoef[m+M]; \
     } \
 \
     for (ltfatInt m=1;m<M;m+=2) \
     { \
        pcoef2[m] = emipi4*pcoef[m]; \
        pcoef2[M2-1-m] = eipi4*pcoef[m]; \
        pcoef2[m+M2] = eipi4*pcoef[m+M];  \
        pcoef2[M4-1-m] = emipi4*pcoef[m+M]; \
     } \
 \
     pcoef+=M2; \
     pcoef2+=M4; \
  }

#define POSTPROC_REAL \
   for(ltfatInt w=0;w<W;w++) \
      for(ltfatInt n=0;n<L;n++) \
         f[n+w*L] = scalconst*CH(creal)(f2[n+w*L]*CH(cexp)(I*PI*n/(2.0*M)));

#define POSTPROC_COMPLEX \
   for(ltfatInt w=0;w<W;w++) \
      for(ltfatInt n=0;n<L;n++) \
         f[n+w*L] = scalconst*f2[n+w*L]*CH(cexp)(I*PI*n/(2.0*M));

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(idwiltiii_long)(const LTFAT_COMPLEX *c, const LTFAT_COMPLEX *g,
                                   const ltfatInt L, const ltfatInt W,
                                   const ltfatInt M, LTFAT_COMPLEX *f)
{
    const ltfatInt N = L / M;
    const ltfatInt M2 = 2 * M;
    const ltfatInt M4 = 4 * M;
    const LTFAT_REAL scalconst = 1.0 / sqrt(2.0);
    const LTFAT_COMPLEX eipi4 = cexp(I * PI / 4.0);
    const LTFAT_COMPLEX emipi4 = cexp(-I * PI / 4.0);

    LTFAT_COMPLEX *coef2 = ltfat_calloc(2 * M * N * W, sizeof * coef2);
    LTFAT_COMPLEX *f2 = ltfat_malloc(L * W * sizeof * f2);


    const LTFAT_COMPLEX *pcoef  = c;
    LTFAT_COMPLEX *pcoef2 = coef2;

    PREPROC_COMPLEX

    LTFAT_NAME(idgt_long)(coef2, g, L, W, M, 2 * M, FREQINV, f2);

    POSTPROC_COMPLEX

    LTFAT_SAFEFREEALL(coef2, f2);
}

LTFAT_EXTERN void
LTFAT_NAME_REAL(idwiltiii_long)(const LTFAT_REAL *c, const LTFAT_REAL *g,
                                const ltfatInt L, const ltfatInt W,
                                const ltfatInt M, LTFAT_REAL *f)
{
    const ltfatInt N = L / M;
    const ltfatInt M2 = 2 * M;
    const ltfatInt M4 = 4 * M;
    const LTFAT_REAL scalconst = 1.0 / sqrt(2.0);
    const LTFAT_COMPLEX eipi4 = cexp(I * PI / 4.0);
    const LTFAT_COMPLEX emipi4 = cexp(-I * PI / 4.0);

    LTFAT_COMPLEX *coef2 = ltfat_calloc(2 * M * N * W, sizeof * coef2);
    LTFAT_COMPLEX *f2 = ltfat_malloc(L * W * sizeof * f2);
    LTFAT_COMPLEX *g2 = ltfat_malloc(L * sizeof * g2);
    for (ltfatInt ii = 0; ii < L; ii++)
        g2[ii] = g[ii];


    const LTFAT_REAL *pcoef  = c;
    LTFAT_COMPLEX *pcoef2 = coef2;

    PREPROC_COMPLEX

    LTFAT_NAME(idgt_long)(coef2, g2, L, W, M, 2 * M, FREQINV, f2);

    POSTPROC_REAL

    LTFAT_SAFEFREEALL(coef2, f2, g2);

}

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(idwiltiii_fb)(const LTFAT_COMPLEX *c, const LTFAT_COMPLEX *g,
                                 const ltfatInt L, const ltfatInt gl,
                                 const ltfatInt W, const ltfatInt M,
                                 LTFAT_COMPLEX *f)
{
    const ltfatInt N = L / M;
    const ltfatInt M2 = 2 * M;
    const ltfatInt M4 = 4 * M;
    const LTFAT_REAL scalconst = 1.0 / sqrt(2.0);
    const LTFAT_COMPLEX eipi4 = cexp(I * PI / 4.0);
    const LTFAT_COMPLEX emipi4 = cexp(-I * PI / 4.0);

    LTFAT_COMPLEX *coef2 = ltfat_calloc(2 * M * N * W, sizeof * coef2);
    LTFAT_COMPLEX *f2 = ltfat_malloc(L * W * sizeof * f2);


    const LTFAT_COMPLEX *pcoef  = c;
    LTFAT_COMPLEX *pcoef2 = coef2;

    PREPROC_COMPLEX

    LTFAT_NAME(idgt_fb)(coef2, g, L, gl, W, M, 2 * M, FREQINV, f2);

    POSTPROC_COMPLEX

    LTFAT_SAFEFREEALL(coef2, f2);

}

LTFAT_EXTERN void
LTFAT_NAME_REAL(idwiltiii_fb)(const LTFAT_REAL *c, const LTFAT_REAL *g,
                              const ltfatInt L, const ltfatInt gl, const ltfatInt W, const ltfatInt M,
                              LTFAT_REAL *f)
{
    const ltfatInt N = L / M;
    const ltfatInt M2 = 2 * M;
    const ltfatInt M4 = 4 * M;
    const LTFAT_REAL scalconst = 1.0 / sqrt(2.0);
    const LTFAT_COMPLEX eipi4 = cexp(I * PI / 4.0);
    const LTFAT_COMPLEX emipi4 = cexp(-I * PI / 4.0);

    LTFAT_COMPLEX *coef2 = ltfat_calloc(2 * M * N * W, sizeof * coef2);
    LTFAT_COMPLEX *f2 = ltfat_malloc(L * W * sizeof * f2);
    LTFAT_COMPLEX *g2 = ltfat_malloc(gl * sizeof * g2);
    for (ltfatInt ii = 0; ii < gl; ii++)
        g2[ii] = g[ii];


    const LTFAT_REAL* pcoef  = c;
    LTFAT_COMPLEX* pcoef2 = coef2;

    PREPROC_COMPLEX

    LTFAT_NAME(idgt_fb)(coef2, g2, L, gl, W, M, 2 * M, FREQINV, f2);

    POSTPROC_REAL

    LTFAT_SAFEFREEALL(coef2, f2, g2);
}

#undef CH
#undef PREPROC_COMPLEX
#undef POSTPROC_REAL
#undef POSTPROC_COMPLEX

