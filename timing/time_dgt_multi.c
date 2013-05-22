#include "config.h"
#include <stdlib.h>
#include "complex.h"
#include "ltfat.h"
#include "ltfat_time.h"

int main( int argc, char *argv[] )
{
   double st0, st1;
   
   if (argc<9)
   {
      printf("Correct parameters: a, M, L, Lg, W, lt1, lt2, nrep\n");     
      return(1);
   }
   const int a = atoi(argv[1]);
   const int M = atoi(argv[2]);
   const int L = atoi(argv[3]);
   const int Lg = atoi(argv[4]);
   const int W = atoi(argv[5]);
   const int lt1 = atoi(argv[6]);
   const int lt2 = atoi(argv[7]);
   const int nrep = atoi(argv[8]);
   
   const int N=L/a;
   
   ltfat_complex *f = ltfat_malloc(L*W*sizeof(ltfat_complex));
   ltfat_complex *g = ltfat_malloc(L*W*sizeof(ltfat_complex));
   ltfat_complex *c = ltfat_malloc(M*N*W*sizeof(ltfat_complex));
   
   LTFAT_NAME(dgt_multi_plan) plan = LTFAT_NAME(dgt_multi_init)(f, g, L, Lg, W, a, M, lt1, lt2, c, FFTW_PATIENT);
   
   st0 = ltfat_time();
   for (int ii=0;ii<nrep;ii++)
   {
      LTFAT_NAME(dgt_multi_execute)(plan);
   }
   st1 = ltfat_time();
   
   LTFAT_NAME(dgt_multi_done)(plan);
   
   printf("%i %i %i %i %i %i %i %f\n",a,M,L,Lg,W,lt1,lt2,(st1-st0)/nrep); 
   
   ltfat_free(f);
   ltfat_free(g);
   ltfat_free(c);
   
   return(0);
}
