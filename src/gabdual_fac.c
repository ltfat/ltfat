#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ltfat.h"

LTFAT_EXTERN void
LTFAT_NAME(gabdual_fac)(const LTFAT_COMPLEX *gf, const int L,
			     const int a, const int M,
			     LTFAT_COMPLEX *gdualf)
{
  
   int h_a, h_m;  

   LTFAT_COMPLEX *Sf;
   
   const int R = 1;

   const LTFAT_COMPLEX zzero = {0.0, 0.0 };
   const LTFAT_COMPLEX alpha = {1.0, 0.0 };

   const int N=L/a;
   
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=N/q;
      
   Sf = (LTFAT_COMPLEX*)ltfat_malloc(p*p*sizeof(LTFAT_COMPLEX));

   /* Copy the contents of gf to gdualf because LAPACK overwrites it input
    * argument
    */
   memcpy(gdualf,gf,sizeof(LTFAT_COMPLEX)*L*R);
  
   for (int w=0;w<R;w++)
   {
      for (int rs=0;rs<c*d;rs++)
      {
	LTFAT_NAME(ltfat_gemm)(CblasNoTrans,CblasConjTrans,p,p,q*R,
		    &alpha,
		    gf+rs*p*q*R,p,
		    gf+rs*p*q*R,p,
		    &zzero,Sf,p);
       
	LTFAT_NAME(ltfat_posv)(p, q*R, Sf, p,
		    gdualf+rs*p*q*R, p);

      }
   }

   /* Clear the work-array. */
   ltfat_free(Sf);


}


LTFAT_EXTERN void
LTFAT_NAME(gabdualreal_fac)(const LTFAT_COMPLEX *gf, const int L,
			     const int a, const int M,
			     LTFAT_COMPLEX *gdualf)
{
  
   int h_a, h_m;  

   LTFAT_COMPLEX *Sf;
   
   const int R = 1;

   const LTFAT_COMPLEX zzero = {0.0, 0.0 };
   const LTFAT_COMPLEX alpha = {1.0, 0.0 };

   const int N=L/a;
   
   const int c=gcd(a, M,&h_a, &h_m);
   const int p=a/c;
   const int q=M/c;
   const int d=N/q;
      
   /* This is a floor operation. */
   const int d2= d/2+1;

   Sf = (LTFAT_COMPLEX*)ltfat_malloc(p*p*sizeof(LTFAT_COMPLEX));

   /* Copy the contents of gf to gdualf because LAPACK overwrites it input
    * argument
    */
   memcpy(gdualf,gf,sizeof(LTFAT_COMPLEX)*L*R);
  
   for (int w=0;w<R;w++)
   {
      for (int rs=0;rs<c*d2;rs++)
      {
	LTFAT_NAME(ltfat_gemm)(CblasNoTrans,CblasConjTrans,p,p,q*R,
		    &alpha,
		    gf+rs*p*q*R,p,
		    gf+rs*p*q*R,p,
		    &zzero,Sf,p);
       
	LTFAT_NAME(ltfat_posv)(p, q*R, Sf, p,
		    gdualf+rs*p*q*R, p);

      }
   }

   /* Clear the work-array. */
   ltfat_free(Sf);


}

