#define LTFAT_DOUBLE
#include <stdlib.h>
#include <stdio.h>
#include "fftw3.h"
#include "ltfat.h"
#include "config.h"
#include "ltfat_time.h"
#include "ltfat_types.h"

/*
FFTW_MEASURE 
FFTW_ESTIMATE
FFTW_PATIENT

You must create the plan before initializing the input, because FFTW_MEASURE 
overwrites the in/out arrays. (Technically, FFTW_ESTIMATE does not touch your arrays,
but you should always create plans first just to be sure.) 
*/

int main( int argc, char *argv[] )
{

   if (argc<2)
   {
      printf("Correct parameters: L, nrep\n");     
      return(1);
   }
   
   const size_t L = atoi(argv[1]);
   const int nrep = atoi(argv[2]);
   
      
   LTFAT_FFTW(plan) p;
   
   LTFAT_REAL *f = (LTFAT_REAL*)ltfat_malloc(L*sizeof(LTFAT_REAL));
   LTFAT_REAL (*c)[2] = (LTFAT_REAL (*)[2])ltfat_malloc((L/2+1)*2*sizeof(LTFAT_REAL));
   
   p = LTFAT_FFTW(plan_dft_r2c_1d)(L, f, c, FFTW_MEASURE);
   
   
   LTFAT_NAME_REAL(fillRand)(f,L);
      
   const double st0=ltfat_time();
   for (int jj=0;jj<nrep;jj++)
   {
      LTFAT_FFTW(execute)(p);
   }
   const double st1=ltfat_time();
   //printf("Length: %i, avr. %f ms \n",L,(st1-st0)/((double)nrep));
   printf("%i, %f\n",L,(st1-st0)/((double)nrep));
   LTFAT_FFTW(destroy_plan)(p);
  
   ltfat_free(f);
   
}
