#include "config.h"
#include "ltfat.h"
#include "ltfat_types.h"

#define CH(name) LTFAT_COMPLEXH_NAME(name)

#define POSTPROC_REAL \
  for (Lint n=0;n<N*W;n+=2) \
  { \
     for (Lint m=0;m<M;m+=2) \
     { \
       pcoef[m]=CH(creal)(pcoef2[m])+CH(cimag)(pcoef2[m]); \
       pcoef[m+M]=CH(creal)(pcoef2[m+M2])-CH(cimag)(pcoef2[m+M2]); \
     } \
     \
     for (Lint m=1;m<M;m+=2) \
     { \
       pcoef[m]=CH(creal)(pcoef2[m])-CH(cimag)(pcoef2[m]); \
       pcoef[m+M]=CH(creal)(pcoef2[m+M2])+CH(cimag)(pcoef2[m+M2]); \
     } \
 \
     pcoef+=M2; \
     pcoef2+=M4; \
  }

#define POSTPROC_COMPLEX \
  for (Lint n=0;n<N*W;n+=2) \
  { \
     for (Lint m=0;m<M;m+=2) \
     { \
         pcoef[m] =   scalconst*(emipi4*pcoef2[m]  +eipi4*pcoef2[M2-1-m]); \
         pcoef[m+M] = scalconst*(eipi4*pcoef2[m+M2]+emipi4*pcoef2[M4-1-m]); \
     } \
 \
     for (Lint m=1;m<M;m+=2) \
     { \
       pcoef[m] = scalconst*(eipi4*pcoef2[m]    +emipi4*pcoef2[M2-1-m]); \
       pcoef[m+M]=scalconst*(emipi4*pcoef2[m+M2]+eipi4*pcoef2[M4-1-m]); \
     } \
 \
     pcoef+=M2; \
     pcoef2+=M4; \
  }

#define PREPROC \
   for(Lint n=0;n<L;n++) \
      f2[n] = (LTFAT_COMPLEX) cexp(-PI*I*n/(2.0*M)); \
   for(Lint w=W-1;w>=0;w--) \
      for(Lint n=0;n<L;n++) \
         f2[n+w*L] = f2[n]*f[n+w*L];


LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(dwiltiii_long)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			   const Lint L, const Lint W, const Lint M,
			   LTFAT_COMPLEX *cout)
{
   const Lint N=L/M;
   const Lint M2=2*M;
   const Lint M4=4*M;
   const LTFAT_REAL scalconst = 1.0/sqrt(2.0);
   const LTFAT_COMPLEX eipi4 = cexp(I*PI/4.0);
   const LTFAT_COMPLEX emipi4 = cexp(-I*PI/4.0);

   LTFAT_COMPLEX *coef2 = (LTFAT_COMPLEX*)ltfat_malloc(2*M*N*W*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *f2 = (LTFAT_COMPLEX*)ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));

   PREPROC

   LTFAT_NAME_COMPLEX(dgt_long)(f2, g, L, W, M, 2*M, coef2);

   LTFAT_COMPLEX *pcoef  = cout;
   LTFAT_COMPLEX *pcoef2 = coef2;

   POSTPROC_COMPLEX

   LTFAT_SAFEFREEALL(coef2,f2);
}

LTFAT_EXTERN void
LTFAT_NAME_REAL(dwiltiii_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
			   const Lint L, const Lint W, const Lint M,
			   LTFAT_REAL *cout)
{
   const Lint N=L/M;
   const Lint M2 = 2*M;
   const Lint M4=4*M;

   LTFAT_COMPLEX *coef2 = (LTFAT_COMPLEX*)ltfat_malloc(2*M*N*W*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *f2 = (LTFAT_COMPLEX*)ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *g2 = (LTFAT_COMPLEX*)ltfat_malloc(L*sizeof(LTFAT_COMPLEX));
   for(Lint ii=0;ii<L;ii++)
       g2[ii]=g[ii];

   PREPROC

   LTFAT_NAME_COMPLEX(dgt_long)(f2, g2, L, W, M, 2*M, coef2);


   LTFAT_REAL *pcoef  = cout;
   LTFAT_COMPLEX *pcoef2 = coef2;

   POSTPROC_REAL

   LTFAT_SAFEFREEALL(coef2,f2,g2);

}

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(dwiltiii_fb)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			    const Lint L, const Lint gl, const Lint W, const Lint M,
			   LTFAT_COMPLEX *cout)
{
   const int N=L/M;
   const int M2=2*M;
   const int M4=4*M;
   const LTFAT_REAL scalconst = 1.0/sqrt(2.0);
   const LTFAT_COMPLEX eipi4 = cexp(I*PI/4.0);
   const LTFAT_COMPLEX emipi4 = cexp(-I*PI/4.0);

   LTFAT_COMPLEX *coef2 = (LTFAT_COMPLEX*)ltfat_malloc(2*M*N*W*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *f2 = (LTFAT_COMPLEX*)ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));

  PREPROC

  /* coef2=comp_dgt(f,g,a,2*M,L); */
  LTFAT_NAME_COMPLEX(dgt_fb)(f2, g, L, gl, W, M, 2*M, coef2);


  LTFAT_COMPLEX *pcoef  = cout;
  LTFAT_COMPLEX *pcoef2 = coef2;

  POSTPROC_COMPLEX

  LTFAT_SAFEFREEALL(coef2,f2);

}

LTFAT_EXTERN void
LTFAT_NAME_REAL(dwiltiii_fb)(const LTFAT_REAL *f, const LTFAT_REAL *g,
			   const Lint L, const Lint gl, const Lint W, const Lint M,
			   LTFAT_REAL *cout)
{
   const Lint N = L/M;
   const Lint M2 = 2*M;
   const Lint M4=4*M;

   LTFAT_COMPLEX *coef2 = (LTFAT_COMPLEX*)ltfat_malloc(2*M*N*W*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *f2 = (LTFAT_COMPLEX*)ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *g2 = (LTFAT_COMPLEX*)ltfat_malloc(gl*sizeof(LTFAT_COMPLEX));
   for(Lint ii=0;ii<gl;ii++)
       g2[ii]=g[ii];

   PREPROC

   LTFAT_NAME_COMPLEX(dgt_fb)(f2, g2, L, gl, W, M, 2*M, coef2);

   LTFAT_REAL* pcoef  = cout;
   LTFAT_COMPLEX* pcoef2 = coef2;

   POSTPROC_REAL

   LTFAT_SAFEFREEALL(coef2,f2,g2);
}

#undef CH
#undef POSTPROC_REAL
#undef POSTPROC_COMPLEX
#undef PREPROC

