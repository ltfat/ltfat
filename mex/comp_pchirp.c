#include "math.h"
#include "mex.h"
#include "config.h"
#include "ltfat.h"
#include "ltfat-mex-helper.h"

#define PI 3.1415926535897932384626433832795

/* Calling convention:
 *  pchirp(L,n);
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
  int l, L, n;
  /*long long workl;*/
   double work, LL, Lpone;
   double *gr,*gi;

   L=(int)mxGetScalar(prhs[0]);
   n=(int)mxGetScalar(prhs[1]);

   plhs[0] = mxCreateDoubleMatrix(L, 1, mxCOMPLEX);
   gr = mxGetPr(plhs[0]);
   gi = mxGetPi(plhs[0]);

   LL=2.0*L;
   Lpone=L+1;

   for (l=0;l<L;l++)
   {
      /* mod(n*m.^2*(L+1),LL); */
      /* workl = (n*m*(L+1))%LL; */
      work = PI*fmod(Lpone*n*l*l,LL)/L;
      gr[l] = cos(work);
      gi[l] = sin(work);

   }

   return;
  
}


