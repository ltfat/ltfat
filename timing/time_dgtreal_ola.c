#include "config.h"
#include <stdlib.h>
#include "ltfat.h"
#include "ltfat_time.h"

int main( int argc, char *argv[] )
{
  double *f, *g;
  ltfat_complex *c;
  int a, M, L, gl, W, N, bl, nrep, ii;
  double s0, s1;

  if (argc<8)
  {
     printf("Correct parameters: a, M, L, W, gl, bl, nrep\n");     
     return(1);
  }
  a = atoi(argv[1]);
  M = atoi(argv[2]);
  L = atoi(argv[3]);
  W = atoi(argv[4]);
  gl = atoi(argv[5]);
  bl = atoi(argv[6]);
  nrep = atoi(argv[7]);

  N=L/a;

  f  = ltfat_malloc(L*W*sizeof(double));
  g  = ltfat_malloc(L*W*sizeof(double));
  c  = ltfat_malloc(M*N*W*sizeof(ltfat_complex));

  LTFAT_NAME(dgtreal_ola_plan) plan = LTFAT_NAME(dgtreal_ola_init)((const double*)g, gl, W, a, M, bl, FFTW_PATIENT);
  
  s0 = ltfat_time();
  for (ii=0;ii<nrep;ii++)
  {
     LTFAT_NAME(dgtreal_ola_execute)(plan,(const double*)f,L,c);
  }
  s1 = ltfat_time();

  LTFAT_NAME(dgtreal_ola_done)(plan);

  printf("%i %i %i %i %i %i %f\n",a,M,L,W,gl,bl,(s1-s0)/nrep); 

  ltfat_free(f);
  ltfat_free(g);
  ltfat_free(c);

  return(0);
}
