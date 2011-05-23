#include "config.h"
#include <math.h>
#include "ltfat.h"


LTFAT_EXTERN void
LTFAT_NAME(gabreassign)(const LTFAT_REAL *s, const LTFAT_REAL *tgrad,
			     const LTFAT_REAL *fgrad, const int L, const int W, 
			     const int a, const int M, LTFAT_REAL *sr)
{

   int ii, posi, posj;


   const int N=L/a;
   const int b=L/M;

   int *timepos = (int*)ltfat_malloc(N*sizeof(int));
   int *freqpos = (int*)ltfat_malloc(M*sizeof(int));
   
   fftindex(N,timepos);
   fftindex(M,freqpos);

   /* Zero the output array. */
   for (ii=0;ii<M*N*W;ii++)
   {
      sr[ii]=0.0;
   }      

   for (int w=0;w<W;w++)
   {
      for (ii=0;ii<M;ii++)
      {
	 for (int jj=0;jj<N;jj++)
	 {
	    /* Do a 'round' followed by a 'mod'. 'round' is not
	     * present in all libraries, so use trunc(x+.5) instead */
	    /*posi=positiverem((int)trunc(tgrad[ii+jj*M]/b+freqpos[ii]+.5),M);
	      posj=positiverem((int)trunc(fgrad[ii+jj*M]/a+timepos[jj]+.5),N);*/
	   posi=positiverem(ltfat_round(tgrad[ii+jj*M]/b+freqpos[ii]),M);
	   posj=positiverem(ltfat_round(fgrad[ii+jj*M]/a+timepos[jj]),N);


	    
	    sr[posi+posj*M]+=s[ii+jj*M];
	 }
      }
   }

   ltfat_free(freqpos);
   ltfat_free(timepos);

}
