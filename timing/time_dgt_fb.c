#include <stdlib.h>
#include "ltfat.h"
#include "ltfat_time.h"

int main( int argc, char *argv[] )
{
  ltfat_complex *f, *g, *c;
  int a, M, L, gl, W, N, nrep, ii;
  double s0, s1;

  if (argc<7)
  {
     printf("Correct parameters: a, M, L, W, gl, nrep\n");     
     return(1);
  }
  a = atoi(argv[1]);
  M = atoi(argv[2]);
  L = atoi(argv[3]);
  W = atoi(argv[4]);
  gl = atoi(argv[5]);
  nrep = atoi(argv[6]);

  N=L/a;

  f  = ltfat_malloc(L*W*sizeof(ltfat_complex));
  g  = ltfat_malloc(L*W*sizeof(ltfat_complex));
  c  = ltfat_malloc(M*N*W*sizeof(ltfat_complex));
  
  s0 = ltfat_time();
  for (ii=0;ii<nrep;ii++)
  {
     dgt_fb((const ltfat_complex*)f,(const ltfat_complex*)g,L,gl,W,a,M,c);
  }
  s1 = ltfat_time();

  printf("%i %i %i %i %i %f\n",a,M,L,W,gl,(s1-s0)/nrep); 

  ltfat_free(f);
  ltfat_free(g);
  ltfat_free(c);

  return(0);
}
