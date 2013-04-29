#include "config.h"
#include <stdlib.h>
#include "ltfat.h"
#include "ltfat_time.h"

int main( int argc, char *argv[] )
{
   ltfat_complex *f, *g, *c;
   int a, M, M2, L, W, N, nrep, ii;
   double s0, s1;
   
   if (argc<5)
   {
      printf("Correct parameters: a, M, L, W, nrep\n");     
      return(1);
   }
   a = atoi(argv[1]);
   M = atoi(argv[2]);
   L = atoi(argv[3]);
   W = atoi(argv[4]);
   nrep = atoi(argv[5]);
   
   N=L/a;
   M2=M/2+1;
   
   f  = ltfat_malloc(L*W*sizeof(double));
   g  = ltfat_malloc(L*W*sizeof(double));
   c  = ltfat_malloc(M2*N*W*sizeof(ltfat_complex));
   
   LTFAT_NAME(dgtreal_long_plan) plan = LTFAT_NAME(dgtreal_long_init)((const double*)f, (const double*)g, 
					      L, W, a, M, c, FFTW_PATIENT);
   
   s0 = ltfat_time();
   for (ii=0;ii<nrep;ii++)
   {
      
      LTFAT_NAME(dgtreal_long_execute)(plan);
      
   }
   s1 = ltfat_time();
   
   LTFAT_NAME(dgtreal_long_done)(plan);
   
   printf("%i %i %i %i %f\n",a,M,L,W,(s1-s0)/nrep); 
   
   ltfat_free(f);
   ltfat_free(g);
   ltfat_free(c);
   
   return(0);
}
