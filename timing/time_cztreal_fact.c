#define LTFAT_DOUBLE
#include <stdlib.h>
#include <stdio.h>
#include "ltfat.h"
#include "ltfat_types.h"
#include "ltfat_time.h"


#ifndef PI
#   define PI 3.141592653589793
#endif

#define CZT_WITH_PLAN
#define CZT_FLAG CZT_NEXTFASTFFT

int main( int argc, char *argv[] )
{
   if (argc<3)
   {
      printf("Correct parameters: L, K, nrep\n");
      return(1);
   }
   int nrep = 1;
   if (argc>3)
   {
     nrep = atoi(argv[3]);
   }

   const size_t L = atoi(argv[1]);
   const size_t K = atoi(argv[2]);

   LTFAT_REAL *f = (LTFAT_REAL*)ltfat_malloc(L*sizeof(LTFAT_REAL));
   LTFAT_COMPLEXH *c = (LTFAT_COMPLEXH*)ltfat_malloc(K*sizeof(LTFAT_COMPLEXH));

   LTFAT_NAME(fillRand)(f,L);


   double o = 0.1;
   double deltao = 2.0*PI/100.0;


   double st0,st1;
   #ifndef CZT_WITH_PLAN
   st0=ltfat_time();
   for (int jj=0;jj<nrep;jj++)
   {
      LTFAT_NAME(chzt_fact)(f,L,1,K,deltao,o,c);
   }
   st1=ltfat_time();
   #else
   LTFAT_NAME(chzt_plan) p = LTFAT_NAME(create_chzt_plan_fact)(K,L,deltao,o,FFTW_MEASURE,CZT_FLAG);
   st0=ltfat_time();
   for (int jj=0;jj<nrep;jj++)
   {
      LTFAT_NAME(chzt_with_plan_fact)(p,f,1,deltao,o,c);
   }
   st1=ltfat_time();
   LTFAT_NAME(destroy_chzt_plan)(p);
   #endif

   //printf("Length: %i, avr. %f ms \n",L,(st1-st0)/((double)nrep));
   printf("%i, %f\n",L,(st1-st0)/((double)nrep));

   ltfat_free(f);
   ltfat_free(c);

   return 0;
}
