#include "config.h"
#include "math.h"
#include "fftw3.h"
#include "ltfat.h"

LTFAT_EXTERN void
LTFAT_NAME(dwilt_long)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			   const int L, const int W, const int M, 
			   LTFAT_COMPLEX *cout)
{
   const int N=L/M;
   const int M2=2*M;
   const int M4=4*M;
   const LTFAT_REAL scalconst = 1.0/sqrt(2.0);
   int m;
   LTFAT_COMPLEX *pcoef, *pcoef2;

   
   LTFAT_COMPLEX *coef2 = (LTFAT_COMPLEX*)ltfat_malloc(2*M*N*W*sizeof(LTFAT_COMPLEX));

  /* coef2=comp_dgt(f,g,a,2*M,L); */
  LTFAT_NAME(dgt_long)(f, g, L, W, M, 2*M, coef2);

  const int nyquestadd = (M%2)*M2;

  pcoef  = cout;
  pcoef2 = coef2;
  for (int n=0;n<N*W;n+=2)
  {
     /*  Unmodulated case. */
     pcoef[0][0]=pcoef2[0][0];
     pcoef[0][1]=pcoef2[0][1];
     
     /* odd value of m */
     for (m=1;m<M;m+=2)
     {
       pcoef[m][0]=scalconst*(-pcoef2[m][1]+pcoef2[M2-m][1]);
       pcoef[m][1]=scalconst*( pcoef2[m][0]-pcoef2[M2-m][0]);
       
       pcoef[m+M][0]=scalconst*(pcoef2[m+M2][0]+pcoef2[M4-m][0]);
       pcoef[m+M][1]=scalconst*(pcoef2[m+M2][1]+pcoef2[M4-m][1]);

     }
     
     /* even value of m */
     for (m=2;m<M;m+=2)
     {
	pcoef[m][0]=scalconst*(pcoef2[m][0]+pcoef2[M2-m][0]);
	pcoef[m][1]=scalconst*(pcoef2[m][1]+pcoef2[M2-m][1]);
	
	pcoef[m+M][0]=scalconst*(-pcoef2[m+M2][1]+pcoef2[M4-m][1]);
	pcoef[m+M][1]=scalconst*( pcoef2[m+M2][0]-pcoef2[M4-m][0]);	

     }

     /* Nyquest */
     pcoef[M][0]=pcoef2[M+nyquestadd][0];
     pcoef[M][1]=pcoef2[M+nyquestadd][1];

     pcoef+=M2;
     pcoef2+=M4;
  }

  ltfat_free(coef2);

}

LTFAT_EXTERN void
LTFAT_NAME(dwiltreal_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
			   const int L, const int W, const int M, 
			   LTFAT_REAL *cout)
{
   const int N=L/M;
   const int coef2_ld = M+1;
   const LTFAT_REAL scalconst = sqrt(2.0);
   int m;
   LTFAT_COMPLEX *pcoef2;
   LTFAT_REAL *pcoef;

   
   /*coef=zeros(2*M,N/2,W);*/

   LTFAT_COMPLEX *coef2 = (LTFAT_COMPLEX*)ltfat_malloc(2*(M+1)*N*W*sizeof(LTFAT_COMPLEX));

  /* coef2=comp_dgt(f,g,a,2*M,L); */
  LTFAT_NAME(dgtreal_long)(f, g, L, W, M, 2*M, coef2);

  const int nyquestadd = (M%2)*coef2_ld;

  pcoef  = cout;
  pcoef2 = coef2;
  for (int n=0;n<N*W;n+=2)
  {
     /*  Unmodulated case. */
     pcoef[0]=pcoef2[0][0];
     
     /* odd value of m */
     for (m=1;m<M;m+=2)
     {
       pcoef[m]=-scalconst*pcoef2[m][1];       
       pcoef[m+M]=scalconst*pcoef2[m+coef2_ld][0];
     }
     
     /* even value of m */
     for (m=2;m<M;m+=2)
     {
       pcoef[m]=scalconst*pcoef2[m][0];
       pcoef[m+M]=-scalconst*pcoef2[m+coef2_ld][1];
     }

     /* Nyquest */
     pcoef[M]=pcoef2[M+nyquestadd][0];

     pcoef+=2*M;
     pcoef2+=2*coef2_ld;
  }

  ltfat_free(coef2);

}

