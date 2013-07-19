/* NOT PROCESSED DIRECTLY, see ltfat_complexindependent.c */
#ifdef LTFAT_TYPE
#include <math.h>
#include <complex.h>
#include "config.h"
#include "ltfat.h"


#ifndef GGA_UNROLL
#   define GGA_UNROLL 8
#endif


#ifndef PI
#   define PI 3.141592653589793
#endif



LTFAT_EXTERN
void LTFAT_NAME(gga)(const LTFAT_TYPE *fPtr, const double *indVecPtr,
                const int L, const int W, const int M, LTFAT_COMPLEXH *cPtr)
{

double pik_term_pre = 2.0*PI/((double) L);
double _Complex cc2_pre = -1.0*I*((double)(L-1));
double _Complex cc_pre = -1.0*I*((double)(L));

#ifndef GGA_UNROLL

for(int w=0;w<W;w++)
{
LTFAT_COMPLEXH *cPtrTmp = (LTFAT_COMPLEXH*) cPtr+w*M;

for(int m=0;m<M;m++)
{
   double pik_term = pik_term_pre*indVecPtr[m];
   LTFAT_REAL cos_pik_term2 = (LTFAT_REAL) cos(pik_term)*2.0;
   LTFAT_COMPLEXH cc = (LTFAT_COMPLEXH) cexp(cc_pre*pik_term);
   LTFAT_COMPLEXH cc2 = (LTFAT_COMPLEXH) cexp(cc2_pre*pik_term);


      LTFAT_TYPE s0 =  0.0;
      LTFAT_TYPE s1 =  0.0;
      LTFAT_TYPE s2 =  0.0;
      LTFAT_TYPE *fPtrTmp = (LTFAT_TYPE*) fPtr+w*L;

      for(int ii=0;ii<L-1;ii++)
      {
         s0 = *fPtrTmp++ + cos_pik_term2*s1 - s2;
         s2=s1;
         s1=s0;
      }
      s0 = *fPtrTmp + cos_pik_term2*s1 - s2;

      *cPtrTmp++ = (s0*cc2 - s1*cc);
   }
}
#else
for(int w=0;w<W;w++)
{
   LTFAT_COMPLEXH *cPtrTmp = (LTFAT_COMPLEXH*) cPtr+w*M;
   int unrollRem = M%GGA_UNROLL;

//#pragma omp parallel for
   for(int m=0;m<M-unrollRem;m+=GGA_UNROLL)
   {
      double pik_term[GGA_UNROLL];
      LTFAT_REAL cos_pik_term2[GGA_UNROLL];
      LTFAT_COMPLEXH cc[GGA_UNROLL];
      LTFAT_COMPLEXH cc2[GGA_UNROLL];

      LTFAT_TYPE s0[GGA_UNROLL];
      LTFAT_TYPE s1[GGA_UNROLL];
      LTFAT_TYPE s2[GGA_UNROLL];

      for(int un=0;un<GGA_UNROLL;un++)
      {
         pik_term[un] = pik_term_pre*indVecPtr[m+un];
         cos_pik_term2[un] = (LTFAT_REAL) cos(pik_term[un])*2.0;
         cc[un] = (LTFAT_COMPLEXH) cexp(cc_pre*pik_term[un]);
         cc2[un] = (LTFAT_COMPLEXH) cexp(cc2_pre*pik_term[un]);
         s0[un] = 0.0;
         s1[un] = 0.0;
         s2[un] = 0.0;
      }

      LTFAT_TYPE *fPtrTmp = (LTFAT_TYPE*) fPtr+w*L;

      for(int ii=0;ii<L-1;ii++)
      {
         for(int un=0;un<GGA_UNROLL;un++)
         {
            s0[un] = *fPtrTmp + cos_pik_term2[un]*s1[un] - s2[un];
            s2[un]=s1[un];
            s1[un]=s0[un];
         }
         fPtrTmp++;
      }
      for(int un=0;un<GGA_UNROLL;un++)
      {
         s0[un] = *fPtrTmp + cos_pik_term2[un]*s1[un] - s2[un];
         cPtrTmp[m+un] = (s0[un]*cc2[un] - s1[un]*cc[un]);
      }
   }

   int m= M-unrollRem;

      double pik_term[GGA_UNROLL];
      LTFAT_REAL cos_pik_term2[GGA_UNROLL];
      LTFAT_COMPLEXH cc[GGA_UNROLL];
      LTFAT_COMPLEXH cc2[GGA_UNROLL];

      LTFAT_TYPE s0[GGA_UNROLL];
      LTFAT_TYPE s1[GGA_UNROLL];
      LTFAT_TYPE s2[GGA_UNROLL];

   for(int un=0;un<unrollRem;un++)
   {
       pik_term[un] = pik_term_pre*indVecPtr[m+un];
       cos_pik_term2[un] = (LTFAT_REAL) cos(pik_term[un])*2.0;
       cc[un] = (LTFAT_COMPLEXH) cexp(cc_pre*pik_term[un]);
       cc2[un] = (LTFAT_COMPLEXH) cexp(cc2_pre*pik_term[un]);
       s0[un] = 0.0;
       s1[un] = 0.0;
       s2[un] = 0.0;
   }

      LTFAT_TYPE *fPtrTmp = (LTFAT_TYPE*) fPtr+w*L;

      for(int ii=0;ii<L-1;ii++)
      {
         for(int un=0;un<unrollRem;un++)
         {
            s0[un] = *fPtrTmp + cos_pik_term2[un]*s1[un] - s2[un];
            s2[un]=s1[un];
            s1[un]=s0[un];
         }
         fPtrTmp++;
      }

      for(int un=0;un<unrollRem;un++)
      {
         s0[un] = *fPtrTmp + cos_pik_term2[un]*s1[un] - s2[un];
         cPtrTmp[m+un] = (s0[un]*cc2[un] - s1[un]*cc[un]);
      }

}
#endif


}

#endif // LTFAT_TYPE

