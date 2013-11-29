#include "config.h"
#include "ltfat.h"
#include "ltfat_types.h"

#define CH(name) LTFAT_COMPLEXH_NAME(name)

#define PREPROC_COMPLEX \
  for (Lint n=0;n<N*W;n+=2) \
  { \
     for (Lint m=0;m<M;m+=2) \
     { \
        pcoef2[m] = eipi4*pcoef[m]; \
        pcoef2[M2-1-m] = emipi4*pcoef[m]; \
        pcoef2[m+M2] = emipi4*pcoef[m+M]; \
        pcoef2[M4-1-m] = eipi4*pcoef[m+M]; \
     } \
 \
     for (Lint m=1;m<M;m+=2) \
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
   for(Lint w=0;w<W;w++) \
      for(Lint n=0;n<L;n++) \
         f[n+w*L] = scalconst*CH(creal)(f2[n+w*L]*CH(cexp)(I*PI*n/(2.0*M)));

#define POSTPROC_COMPLEX \
   for(Lint w=0;w<W;w++) \
      for(Lint n=0;n<L;n++) \
         f[n+w*L] = scalconst*f2[n+w*L]*CH(cexp)(I*PI*n/(2.0*M));

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(idwiltiii_long)(const LTFAT_COMPLEX *c, const LTFAT_COMPLEX *g,
			   const Lint L, const Lint W, const Lint M,
			   LTFAT_COMPLEX *f)
{
   const Lint N=L/M;
   const Lint M2=2*M;
   const Lint M4=4*M;
   const LTFAT_REAL scalconst = 1.0/sqrt(2.0);
   const LTFAT_COMPLEX eipi4 = cexp(I*PI/4.0);
   const LTFAT_COMPLEX emipi4 = cexp(-I*PI/4.0);

   LTFAT_COMPLEX *coef2 = (LTFAT_COMPLEX*)ltfat_calloc(2*M*N*W,sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *f2 = (LTFAT_COMPLEX*)ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));


   LTFAT_COMPLEX *pcoef  = c;
   LTFAT_COMPLEX *pcoef2 = coef2;

   PREPROC_COMPLEX

   LTFAT_NAME(idgt_long)(coef2, g, L, W, M, 2*M, f2);
   
   POSTPROC_COMPLEX

   LTFAT_SAFEFREEALL(coef2,f2);
}

LTFAT_EXTERN void
LTFAT_NAME_REAL(idwiltiii_long)(const LTFAT_REAL *c, const LTFAT_REAL *g,
			   const Lint L, const Lint W, const Lint M,
			   LTFAT_REAL *f)
{
   const Lint N=L/M;
   const Lint M2 = 2*M;
   const Lint M4=4*M;
   const LTFAT_REAL scalconst = 1.0/sqrt(2.0);
   const LTFAT_COMPLEX eipi4 = cexp(I*PI/4.0);
   const LTFAT_COMPLEX emipi4 = cexp(-I*PI/4.0);

   LTFAT_COMPLEX *coef2 = (LTFAT_COMPLEX*)ltfat_calloc(2*M*N*W,sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *f2 = (LTFAT_COMPLEX*)ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *g2 = (LTFAT_COMPLEX*)ltfat_malloc(L*sizeof(LTFAT_COMPLEX));
   for(Lint ii=0;ii<L;ii++)
       g2[ii]=g[ii];


   LTFAT_REAL *pcoef  = c;
   LTFAT_COMPLEX *pcoef2 = coef2;

   PREPROC_COMPLEX

   LTFAT_NAME(idgt_long)(coef2, g2, L, W, M, 2*M, f2);

   POSTPROC_REAL

   LTFAT_SAFEFREEALL(coef2,f2,g2);

}

LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(idwiltiii_fb)(const LTFAT_COMPLEX *c, const LTFAT_COMPLEX *g,
			    const Lint L, const Lint gl, const Lint W, const Lint M,
			   LTFAT_COMPLEX *f)
{
   const int N=L/M;
   const int M2=2*M;
   const int M4=4*M;
   const LTFAT_REAL scalconst = 1.0/sqrt(2.0);
   const LTFAT_COMPLEX eipi4 = cexp(I*PI/4.0);
   const LTFAT_COMPLEX emipi4 = cexp(-I*PI/4.0);

   LTFAT_COMPLEX *coef2 = (LTFAT_COMPLEX*)ltfat_calloc(2*M*N*W,sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *f2 = (LTFAT_COMPLEX*)ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));


   LTFAT_COMPLEX *pcoef  = c;
   LTFAT_COMPLEX *pcoef2 = coef2;

   PREPROC_COMPLEX
   
   LTFAT_NAME(idgt_fb)(coef2, g, L, gl, W, M, 2*M, f2);
   
   POSTPROC_COMPLEX

   LTFAT_SAFEFREEALL(coef2,f2);

}

LTFAT_EXTERN void
LTFAT_NAME_REAL(idwiltiii_fb)(const LTFAT_REAL *c, const LTFAT_REAL *g,
			   const Lint L, const Lint gl, const Lint W, const Lint M,
			   LTFAT_REAL *f)
{
   const Lint N = L/M;
   const Lint M2 = 2*M;
   const Lint M4=4*M;
   const LTFAT_REAL scalconst = 1.0/sqrt(2.0);
   const LTFAT_COMPLEX eipi4 = cexp(I*PI/4.0);
   const LTFAT_COMPLEX emipi4 = cexp(-I*PI/4.0);

   LTFAT_COMPLEX *coef2 = (LTFAT_COMPLEX*)ltfat_calloc(2*M*N*W,sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *f2 = (LTFAT_COMPLEX*)ltfat_malloc(L*W*sizeof(LTFAT_COMPLEX));
   LTFAT_COMPLEX *g2 = (LTFAT_COMPLEX*)ltfat_malloc(gl*sizeof(LTFAT_COMPLEX));
   for(Lint ii=0;ii<gl;ii++)
       g2[ii]=g[ii];


   LTFAT_REAL* pcoef  = c;
   LTFAT_COMPLEX* pcoef2 = coef2;

   PREPROC_COMPLEX
   
   LTFAT_NAME(idgt_fb)(coef2, g2, L, gl, W, M, 2*M, f2);
   
   POSTPROC_REAL
   
   LTFAT_SAFEFREEALL(coef2,f2,g2);
}

#undef CH
#undef PREPROC_COMPLEX
#undef POSTPROC_REAL
#undef POSTPROC_COMPLEX

