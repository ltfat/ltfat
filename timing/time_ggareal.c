#define LTFAT_DOUBLE
#include <stdlib.h>
#include <stdio.h>
#include "ltfat.h"
#include "ltfat_types.h"
#include "ltfat_time.h"


int main( int argc, char *argv[] )
{
   if (argc<3)
   {
      printf("Correct parameters: L, M, nrep\n");     
      return(1);
   }
   int nrep = 1;
   if (argc>3)
   {
     nrep = atoi(argv[3]);
   }
   
   const int L = atoi(argv[1]);
   const int M = atoi(argv[2]);
   
   
   LTFAT_REAL *f = (LTFAT_REAL *)ltfat_malloc(L*sizeof(LTFAT_REAL));
   LTFAT_COMPLEXH *c = (LTFAT_COMPLEXH*)ltfat_malloc(M*sizeof(LTFAT_COMPLEXH));
   double *indVec = (double*)ltfat_malloc(M*sizeof(double));
   
   LTFAT_NAME_REAL(fillRand)(f,L);

   for(int m=0;m<M;m++)
   {
      indVec[m]=m;
   }
      
   const double st0=ltfat_time();

   for (int jj=0;jj<nrep;jj++)
   {
      LTFAT_NAME_REAL(gga)(f,indVec,L,1,M,c);
   }
   
   const double st1=ltfat_time();
   //printf("L=%d, M=%d, nrep=%d\n",L,M,nrep);
   //printf("Length: %i, avr. %f ms \n",L,(st1-st0)/((double)nrep));
    printf("%i, %f\n",L,(st1-st0)/((double)nrep));
   
  
   ltfat_free(f);
   ltfat_free(c);
   ltfat_free(indVec);
   
}
