#include "ltfat.h"
#include "ltfat_types.h"

#define CH(name) LTFAT_COMPLEXH(name)

#define POSTPROC_REAL \
  for (ltfatInt n=0;n<N*W;n+=2) \
  { \
     pcoef[0]=CH(creal)(pcoef2[0]); \
\
     for (ltfatInt m=1;m<M;m+=2) \
     { \
       pcoef[m]=-scalconst*CH(cimag)(pcoef2[m]); \
       pcoef[m+M]=scalconst*CH(creal)(pcoef2[m+coef2_ld]); \
     } \
 \
     for (ltfatInt m=2;m<M;m+=2) \
     { \
       pcoef[m]=scalconst*CH(creal)(pcoef2[m]); \
       pcoef[m+M]=-scalconst*CH(cimag)(pcoef2[m+coef2_ld]); \
     } \
 \
     pcoef[M]=CH(creal)(pcoef2[M+nyquestadd]); \
     pcoef+=2*M; \
     pcoef2+=2*coef2_ld; \
  }

#define POSTPROC_COMPLEX \
  for (ltfatInt n=0;n<N*W;n+=2) \
  { \
     pcoef[0] = pcoef2[0]; \
 \
     for (ltfatInt m=1;m<M;m+=2) \
     { \
       pcoef[m] = scalconst*I*(pcoef2[m]-pcoef2[M2-m]); \
       pcoef[m+M]=scalconst*(pcoef2[m+M2]+pcoef2[M4-m]); \
     } \
 \
     for (ltfatInt m=2;m<M;m+=2) \
     { \
         pcoef[m] = scalconst*(pcoef2[m]+pcoef2[M2-m]); \
         pcoef[m+M] = scalconst*I*(pcoef2[m+M2]-pcoef2[M4-m]); \
     } \
 \
     pcoef[M]=pcoef2[M+nyquestadd]; \
     pcoef+=M2; \
     pcoef2+=M4; \
  }


LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(dwilt_long)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                               const ltfatInt L, const ltfatInt W,
                               const ltfatInt M, LTFAT_COMPLEX *cout)
{
    const ltfatInt N=L/M;
    const ltfatInt M2=2*M;
    const ltfatInt M4=4*M;
    const LTFAT_REAL scalconst = 1.0/sqrt(2.0);

    LTFAT_COMPLEX *coef2 = ltfat_malloc(2*M*N*W*sizeof*coef2);

    /* coef2=comp_dgt(f,g,a,2*M,L); */
    LTFAT_NAME_COMPLEX(dgt_long)(f, g, L, W, M, 2*M, FREQINV, coef2);

    const ltfatInt nyquestadd = (M%2)*M2;

    LTFAT_COMPLEX *pcoef  = cout;
    LTFAT_COMPLEX *pcoef2 = coef2;

    POSTPROC_COMPLEX

    ltfat_free(coef2);

}

LTFAT_EXTERN void
LTFAT_NAME_REAL(dwilt_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                            const ltfatInt L, const ltfatInt W,
                            const ltfatInt M, LTFAT_REAL *cout)
{
    const ltfatInt N=L/M;
    const ltfatInt coef2_ld = M+1;
    const LTFAT_REAL scalconst = sqrt(2.0);
    const ltfatInt nyquestadd = (M%2)*coef2_ld;

    LTFAT_COMPLEX *coef2 = ltfat_malloc((M+1)*N*W*sizeof*coef2);

    /* coef2=comp_dgt(f,g,a,2*M,L); */
    LTFAT_NAME(dgtreal_long)(f, g, L, W, M, 2*M, FREQINV,coef2);


    LTFAT_REAL *pcoef  = cout;
    LTFAT_COMPLEX *pcoef2 = coef2;

    POSTPROC_REAL

    ltfat_free(coef2);

}

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(dwilt_fb)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                             const ltfatInt L, const ltfatInt gl,
                             const ltfatInt W, const ltfatInt M,
                             LTFAT_COMPLEX *cout)
{
    const ltfatInt N=L/M;
    const ltfatInt M2=2*M;
    const ltfatInt M4=4*M;
    const LTFAT_REAL scalconst = 1.0/sqrt(2.0);

    LTFAT_COMPLEX *coef2 = ltfat_malloc(2*M*N*W*sizeof*coef2);

    /* coef2=comp_dgt(f,g,a,2*M,L); */
    LTFAT_NAME_COMPLEX(dgt_fb)(f, g, L, gl, W, M, 2*M, FREQINV, coef2);

    const ltfatInt nyquestadd = (M%2)*M2;

    LTFAT_COMPLEX *pcoef  = cout;
    LTFAT_COMPLEX *pcoef2 = coef2;

    POSTPROC_COMPLEX

    ltfat_free(coef2);

}

LTFAT_EXTERN void
LTFAT_NAME_REAL(dwilt_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                          const ltfatInt L, const ltfatInt gl,
                          const ltfatInt W, const ltfatInt M,
                          LTFAT_REAL *cout)
{
    const ltfatInt N = L/M;
    const ltfatInt coef2_ld = M + 1;
    const ltfatInt nyquestadd = (M%2)*coef2_ld;
    const LTFAT_REAL scalconst = (LTFAT_REAL) sqrt(2.0);

    LTFAT_COMPLEX *coef2 = ltfat_malloc((M+1)*N*W*sizeof*coef2);
    LTFAT_NAME(dgtreal_fb)(f, g, L, gl, W, M, 2*M, FREQINV, coef2);

    LTFAT_REAL* pcoef  = cout;
    LTFAT_COMPLEX* pcoef2 = coef2;

    POSTPROC_REAL

    ltfat_free(coef2);
}

#undef CH
#undef POSTPROC_REAL
#undef POSTPROC_COMPLEX

