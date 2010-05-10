#include "config.h"
#include <math.h>

const double pi = 3.141592653589793;

void LTFAT_NAME(pgauss)(const int L, const double w, const double center,
	    LTFAT_REAL *g)
{
   
   int lr,k,nk;
   double tmp,sqrtl, safe, gnorm;
   
   sqrtl=sqrt((double)L);
   safe=4;
   gnorm=0;
   
   /* Outside the interval [-safe,safe] then exp(-pi*x.^2) is numerically zero. */
   nk=(int)ceil(safe/sqrt((double)L/sqrt(w)));
   
   for ( lr=0; lr<L; lr++)
   {
      g[lr]=0.0;
      for (k=-nk;k<=nk; k++)
      {  
	 /* Use a tmp variable to calculate squaring */
	 tmp = ((double)lr+center)/sqrtl-(double)k*sqrtl;
	 g[lr]+=exp(-pi*tmp*tmp/w);
      }
      gnorm +=g[lr]*g[lr];
   }
   
   /* Normalize it exactly. */
   gnorm=sqrt(gnorm);

   for ( lr=0; lr<L; lr++)
   {   
      g[lr] /= gnorm;
   }   
}
