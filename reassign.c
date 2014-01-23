#include "ltfat.h"
#include "ltfat_types.h"


LTFAT_EXTERN void
LTFAT_NAME(gabreassign)(const LTFAT_REAL *s, const LTFAT_REAL *tgrad,
                        const LTFAT_REAL *fgrad, const ltfatInt L, const ltfatInt W,
                        const ltfatInt a, const ltfatInt M, LTFAT_REAL *sr)
{

    ltfatInt ii, posi, posj;


    const ltfatInt N=L/a;
    const ltfatInt b=L/M;

    ltfatInt *timepos = ltfat_malloc(N*sizeof*timepos);
    ltfatInt *freqpos = ltfat_malloc(M*sizeof*freqpos);

    fftindex(N,timepos);
    fftindex(M,freqpos);

    /* Zero the output array. */
    memset(sr,0,M*N*W*sizeof*sr);

    for (ltfatInt w=0; w<W; w++)
    {
        for (ii=0; ii<M; ii++)
        {
            for (ltfatInt jj=0; jj<N; jj++)
            {
                /* Do a 'round' followed by a 'mod'. 'round' is not
                 * present in all libraries, so use trunc(x+.5) instead */
                /*posi=positiverem((ltfatInt)trunc(tgrad[ii+jj*M]/b+freqpos[ii]+.5),M);
                  posj=positiverem((ltfatInt)trunc(fgrad[ii+jj*M]/a+timepos[jj]+.5),N);*/
                posi=positiverem(ltfat_round(tgrad[ii+jj*M]/b+freqpos[ii]),M);
                posj=positiverem(ltfat_round(fgrad[ii+jj*M]/a+timepos[jj]),N);

                sr[posi+posj*M]+=s[ii+jj*M];
            }
        }
    }

    LTFAT_SAFEFREEALL(freqpos,timepos);
}
