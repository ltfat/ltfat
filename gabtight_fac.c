#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ltfat.h"

LTFAT_EXTERN void
LTFAT_NAME(gabtight_fac)(const LTFAT_COMPLEX *gf, const int L,
			      const int a, const int M,
			      LTFAT_COMPLEX *gtightf)
{
  
   int h_a, h_m;
   
   int w, rs;

   LTFAT_COMPLEX *Sf, *U, *VT, *gfwork;
   LTFAT_REAL *S;

   const LTFAT_COMPLEX zzero = {0.0, 0.0 };
   const LTFAT_COMPLEX alpha = {1.0, 0.0 };
   
   const int R = 1;

   const int N=L/a;
   
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=N/q;
      
   S  = (LTFAT_REAL*)ltfat_malloc(p*sizeof(LTFAT_REAL));
   Sf = (LTFAT_COMPLEX*)ltfat_malloc(p*p*sizeof(LTFAT_COMPLEX));
   U  = (LTFAT_COMPLEX*)ltfat_malloc(p*p*sizeof(LTFAT_COMPLEX));
   VT = (LTFAT_COMPLEX*)ltfat_malloc(p*q*R*sizeof(LTFAT_COMPLEX));
   gfwork = (LTFAT_COMPLEX*)ltfat_malloc(L*R*sizeof(LTFAT_COMPLEX));

   /* Copy the contents of gf to gfwork because LAPACK overwrites
    * the input.
    */
   memcpy(gfwork,gf,sizeof(LTFAT_COMPLEX)*L*R);
  
   for (w=0;w<R;w++)
   {
      for (rs=0;rs<c*d;rs++)
      {
	/* Compute the thin SVD */
	LTFAT_NAME(ltfat_gesvd)(p, q*R, gfwork+rs*p*q*R, p,
		     S, U, p, VT, p);

	/* Combine U and V. */
	LTFAT_NAME(ltfat_gemm)(CblasNoTrans,CblasNoTrans,p,q*R,p,
			       &alpha,(const LTFAT_COMPLEX*)U,p,
			       (const LTFAT_COMPLEX*)VT,p,
			       &zzero,gtightf+rs*p*q*R, p);


      }
   }

   ltfat_free(gfwork);
   ltfat_free(Sf);
   ltfat_free(S);
   ltfat_free(U);
   ltfat_free(VT);


}


LTFAT_EXTERN void
LTFAT_NAME(gabtightreal_fac)(const LTFAT_COMPLEX *gf, const int L,
			      const int a, const int M,
			      LTFAT_COMPLEX *gtightf)
{
  
   int h_a, h_m;
   
   int w, rs;

   LTFAT_COMPLEX *Sf, *U, *VT, *gfwork;
   LTFAT_REAL *S;

   const LTFAT_COMPLEX zzero = {0.0, 0.0 };
   const LTFAT_COMPLEX alpha = {1.0, 0.0 };
   
   const int R = 1;

   const int N=L/a;
   
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=N/q;

   /* This is a floor operation. */
   const int d2= d/2+1;
      
   S  = (LTFAT_REAL*)ltfat_malloc(p*sizeof(LTFAT_REAL));
   Sf = (LTFAT_COMPLEX*)ltfat_malloc(p*p*sizeof(LTFAT_COMPLEX));
   U  = (LTFAT_COMPLEX*)ltfat_malloc(p*p*sizeof(LTFAT_COMPLEX));
   VT = (LTFAT_COMPLEX*)ltfat_malloc(p*q*R*sizeof(LTFAT_COMPLEX));
   gfwork = (LTFAT_COMPLEX*)ltfat_malloc(L*R*sizeof(LTFAT_COMPLEX));

   /* Copy the contents of gf to gfwork because LAPACK overwrites
    * the input.
    */
   memcpy(gfwork,gf,sizeof(LTFAT_COMPLEX)*L*R);
  
   for (w=0;w<R;w++)
   {
      for (rs=0;rs<c*d2;rs++)
      {
	/* Compute the thin SVD */
	LTFAT_NAME(ltfat_gesvd)(p, q*R, gfwork+rs*p*q*R, p,
		     S, U, p, VT, p);

	/* Combine U and V. */
	LTFAT_NAME(ltfat_gemm)(CblasNoTrans,CblasNoTrans,p,q*R,p,
			       &alpha,(const LTFAT_COMPLEX*)U,p,
			       (const LTFAT_COMPLEX*)VT,p,
			       &zzero,gtightf+rs*p*q*R, p);

      }
   }

   ltfat_free(gfwork);
   ltfat_free(Sf);
   ltfat_free(S);
   ltfat_free(U);
   ltfat_free(VT);


}

