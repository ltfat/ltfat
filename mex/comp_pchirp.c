#include "math.h"
#include "complex.h"
#include "mex.h"
#include "config.h"
#include "ltfat.h"

#define PI 3.1415926535897932384626433832795

static inline long positiverem_long(long a,long b)
{
  const long c = a%b;
  return(c<0 ? c+b : c);
}

/* Calling convention:
 *  pchirp(L,n);
 */
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
{ 
   const long L=(long)mxGetScalar(prhs[0]);
   const long n=(long)mxGetScalar(prhs[1]);

   plhs[0] = mxCreateDoubleMatrix(L, 1, mxCOMPLEX);
   double *gr = mxGetPr(plhs[0]);
   double *gi = mxGetPi(plhs[0]);


   const long LL=2*L;
   const long Lponen=positiverem_long((L+1)*n,LL);
   
   for (long m=0;m<L;m++)
   {
      const long idx = positiverem_long(
   	 positiverem_long(Lponen*m,LL)*m,LL);

      gr[m] = cos(PI*idx/L);
      gi[m] = sin(PI*idx/L);
   }

   return;
  
}


