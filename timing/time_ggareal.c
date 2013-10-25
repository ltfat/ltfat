#define LTFAT_DOUBLE
#include <stdlib.h>
#include <stdio.h>
#include "ltfat.h"
#include "ltfat_types.h"
#include "ltfat_time.h"

#define GGA_WITH_PLAN

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
   
   const size_t L = atoi(argv[1]);
   const size_t M = atoi(argv[2]);
   
   //printf("1) Before allocation\n");
   LTFAT_REAL *f = (LTFAT_REAL *)ltfat_malloc(L*sizeof(LTFAT_REAL));
   LTFAT_COMPLEXH *c = (LTFAT_COMPLEXH*)ltfat_malloc(M*sizeof(LTFAT_COMPLEXH));
   double *indVec = (double*)ltfat_malloc(M*sizeof(double));
   //printf("2) Before filling\n");
   LTFAT_NAME_REAL(fillRand)(f,L);

   for(size_t m=0;m<M;m++)
   {
      indVec[m]=0.1+((double)m)/100.0;
   }
      
	  
   //printf("3) Before call\n");
   double st0,st1;
   #ifndef GGA_WITH_PLAN
      st0=ltfat_time();
      for (int jj=0;jj<nrep;jj++)
      {
         LTFAT_NAME_REAL(gga)(f,indVec,L,1,M,c);
      }
      st1=ltfat_time();
   #else
      //printf("4) Before plan\n");
      LTFAT_NAME_REAL(gga_plan) p = LTFAT_NAME_REAL(create_gga_plan)(indVec,M,L);
      st0=ltfat_time();  

	  //printf("5) Before loop\n");
      for (int jj=0;jj<nrep;jj++)
      {
         LTFAT_NAME_REAL(gga_with_plan)(p,f,c,1);
      }
	  st1=ltfat_time();
	  //printf("6) Before plan delete\n");
      LTFAT_NAME_REAL(destroy_gga_plan)(p);

   #endif
   
   
   //printf("L=%d, M=%d, nrep=%d\n",L,M,nrep);
   //printf("Length: %i, avr. %f ms \n",L,(st1-st0)/((double)nrep));
    printf("%i, %f\n",L,(st1-st0)/((double)nrep));
   
  
   ltfat_free(f);
   ltfat_free(c);
   ltfat_free(indVec);
   return 0;
}
