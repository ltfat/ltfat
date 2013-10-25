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

#define FFTW_OPTIM FFTW_ESTIMATE

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
   
   LTFAT_REAL (*f)[2] = (LTFAT_REAL (*)[2])ltfat_malloc(L*2*sizeof(LTFAT_REAL));
   
   //const double pt0=ltfat_time();
   p = LTFAT_FFTW(plan_dft_1d)(L, f, f, FFTW_FORWARD, FFTW_OPTIM);
   //const double pt1=ltfat_time();
   //printf("Plan: %f\n",(pt1-pt0)/((double)nrep));
   
   LTFAT_NAME_COMPLEX(fillRand)((double _Complex*)f,L);
      
   const double st0=ltfat_time();
   for (int jj=0;jj<nrep;jj++)
   {
      LTFAT_FFTW(execute)(p);
   }
   const double st1=ltfat_time();
   printf("%i, %f\n",L,(st1-st0)/((double)nrep));
   
   LTFAT_FFTW(destroy_plan)(p);
   /*
   LTFAT_FFTW(plan) p2;
   const double p2t0=ltfat_time();
   p2 = LTFAT_FFTW(plan_dft_1d)(L/2, f, f, FFTW_FORWARD, FFTW_OPTIM);
   const double p2t1=ltfat_time();
   printf("Plan2: %f\n",(p2t1-p2t0)/((double)nrep));
  
   LTFAT_FFTW(destroy_plan)(p2);
   
   LTFAT_FFTW(plan) p3;
   const double p3t0=ltfat_time();
   p3 = LTFAT_FFTW(plan_dft_1d)(L, f, f, FFTW_FORWARD, FFTW_OPTIM);
   const double p3t1=ltfat_time();
   printf("Plan2: %f\n",(p3t1-p3t0)/((double)nrep));
  
   LTFAT_FFTW(destroy_plan)(p3);
   */
 
   ltfat_free(f);
   
}
