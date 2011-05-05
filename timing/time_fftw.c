#define LTFAT_DOUBLE
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include "fftw3.h"
#include "ltfat.h"
#include "ltfat_time.h"


int main( int argc, char *argv[] )
{
   if (argc<2)
   {
      printf("Correct parameters: L, nrep\n");     
      return(1);
   }
   
   const int L = atoi(argv[1]);
   const int nrep = atoi(argv[2]);
   
   LTFAT_FFTW(plan) p;
   
   ltfat_complex *f = (ltfat_complex*)ltfat_malloc(L*sizeof(ltfat_complex));
   
   p = LTFAT_FFTW(plan_dft_1d)(L, f, f, FFTW_FORWARD, FFTW_MEASURE);
      
   const double st0=ltfat_time();
   for (int jj=0;jj<nrep;jj++)
   {
      LTFAT_FFTW(execute)(p);
   }
   
   const double st1=ltfat_time();
   printf("%i %f\n",L,(st1-st0)/nrep);
   
   LTFAT_FFTW(destroy_plan)(p);
  
   ltfat_free(f);
   
}
