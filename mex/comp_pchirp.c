#include "math.h"
#include "complex.h"
#include "mex.h"
#include "config.h"
#include "ltfat.h"

#define PI 3.1415926535897932384626433832795

static inline int positiverem_long(long int a,int b)
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
   const int L=(int)mxGetScalar(prhs[0]);
   const int n=(int)mxGetScalar(prhs[1]);

   plhs[0] = mxCreateDoubleMatrix(L, 1, mxCOMPLEX);
   double *gr = mxGetPr(plhs[0]);
   double *gi = mxGetPi(plhs[0]);

   const double LL=2.0*L;
   const double Lpone=L+1;
   
   for (int m=0;m<L;m++)
   {
      const double work = PI*fmod(fmod(fmod(Lpone*n,LL)*m,LL)*m,LL)/L;
      gr[m] = cos(work);
      gi[m] = sin(work);
   }

   return;
  
}


