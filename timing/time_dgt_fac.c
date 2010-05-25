#include <stdlib.h>
#include "ltfat.h"
#include "ltfat_time.h"

int main( int argc, char *argv[] )
{
  ltfat_complex *f, *gf, *c;
  int a, M, L, W, N, nrep, ii;
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

  f  = ltfat_malloc(L*W*sizeof(ltfat_complex));
  gf = ltfat_malloc(L*W*sizeof(ltfat_complex));
  c  = ltfat_malloc(M*N*W*sizeof(ltfat_complex));
  
  dgt_plan plan = plan_dgt_fac(L, W, a, M, c, FFTW_PATIENT);

  s0 = ltfat_time();
  for (ii=0;ii<nrep;ii++)
  {

    execute_dgt_fac(plan,f,gf);
    
  }
  s1 = ltfat_time();

  destroy_dgt_fac(plan);

  printf("%i %i %i %i %f\n",a,M,L,W,(s1-s0)/nrep); 

  ltfat_free(f);
  ltfat_free(gf);
  ltfat_free(c);

  return(0);
}
