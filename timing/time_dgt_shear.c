#include <stdlib.h>
#include "ltfat.h"
#include "ltfat_time.h"

int main( int argc, char *argv[] )
{
   double st0, st1;
   
   if (argc<9)
   {
      printf("Correct parameters: a, M, L, W, s0, s1, br, nrep\n");     
      return(1);
   }
   const int a = atoi(argv[1]);
   const int M = atoi(argv[2]);
   const int L = atoi(argv[3]);
   const int W = atoi(argv[4]);
   const int s0 = atoi(argv[5]);
   const int s1 = atoi(argv[6]);
   const int br = atoi(argv[7]);
   const int nrep = atoi(argv[8]);
   
   const int N=L/a;
   
   const ltfat_complex *f = ltfat_malloc(L*W*sizeof(ltfat_complex));
   const ltfat_complex *g = ltfat_malloc(L*W*sizeof(ltfat_complex));
   ltfat_complex       *c = ltfat_malloc(M*N*W*sizeof(ltfat_complex));
   
   dgt_shear_plan plan = dgt_shear_init(f, g, L, W, a, M, s0, s1, br, c, FFTW_PATIENT);
   
   st0 = ltfat_time();
   for (int ii=0;ii<nrep;ii++)
   {
      dgt_shear_execute(plan);
   }
   st1 = ltfat_time();
   
   dgt_shear_done(plan);
   
   printf("%i %i %i %i %i %i %i %f\n",a,M,L,W,s0,s1,br,(st1-st0)/nrep); 
   
   ltfat_free(f);
   ltfat_free(g);
   ltfat_free(c);
   
   return(0);
}
