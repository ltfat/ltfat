/* NOT PROCESSED DIRECTLY, dst_ci.c */
#ifdef LTFAT_TYPE

#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN LTFAT_FFTW(plan)
LTFAT_NAME(dst_init)( const ltfatInt L, const ltfatInt W, LTFAT_TYPE *cout,
                      const dst_kind kind)
{
    LTFAT_FFTW(iodim) dims, howmanydims;
    LTFAT_FFTW(plan) p;

#ifdef LTFAT_COMPLEXTYPE
    dims.n = L;
    dims.is = 2;
    dims.os = 2;

    howmanydims.n = W;
    howmanydims.is = 2*L;
    howmanydims.os = 2*L;

    unsigned flag = FFTW_ESTIMATE | FFTW_UNALIGNED;
#else
    dims.n = L;
    dims.is = 1;
    dims.os = 1;

    howmanydims.n = W;
    howmanydims.is = L;
    howmanydims.os = L;

    unsigned flag = FFTW_ESTIMATE;
#endif

    LTFAT_FFTW(r2r_kind) kindFftw = (LTFAT_FFTW(r2r_kind)) kind;
    p = LTFAT_FFTW(plan_guru_r2r)(1, &dims,
                                  1, &howmanydims,
                                  (LTFAT_REAL*)cout, (LTFAT_REAL*)cout,
                                  &kindFftw, flag);

    return p;
}


// f and cout cannot be equal, because creating plan can tamper with the array
LTFAT_EXTERN void
LTFAT_NAME(dst)(const LTFAT_TYPE *f, const ltfatInt L, const ltfatInt W,
                LTFAT_TYPE *cout, const dst_kind kind)
{
    LTFAT_FFTW(plan) p = LTFAT_NAME(dst_init)( L, W, cout, kind);

    LTFAT_NAME(dst_execute)(p, f,  L,  W, cout, kind);

    LTFAT_FFTW(destroy_plan)(p);
}

// f and cout can be equal, provided plan was already created
LTFAT_EXTERN void
LTFAT_NAME(dst_execute)(LTFAT_FFTW(plan) p, const LTFAT_TYPE *f,
                        const ltfatInt L, const ltfatInt W, LTFAT_TYPE *cout,
                        const dst_kind kind)
{
    // Copy input to the output
    if(cout!=f)
        memcpy(cout,f,L*W*sizeof*f);

    if(L==1)
        return;

    ltfatInt N = 2*L;
    LTFAT_REAL sqrt2 = (LTFAT_REAL) sqrt(2.0);
    LTFAT_REAL postScale = (LTFAT_REAL) 1.0/sqrt2;
    LTFAT_REAL scale = (LTFAT_REAL) sqrt2*(1.0/(double)N)*sqrt((double)L);

    if(kind==DSTIII)
    {
        for(ltfatInt ii=0; ii<W; ii++)
        {
            cout[(ii+1)*L-1] *= sqrt2;
        }
    }

    if(kind==DSTI)
    {
        N += 2;
        scale = (LTFAT_REAL) sqrt2*(1.0/((double)N))*sqrt((double)L+1);
    }

    LTFAT_REAL* c_r = (LTFAT_REAL*)cout;

    LTFAT_FFTW(execute_r2r)(p,c_r,c_r);
#ifdef LTFAT_COMPLEXTYPE
    LTFAT_REAL* c_i = c_r+1;
    LTFAT_FFTW(execute_r2r)(p,c_i,c_i);
#endif

    // Post-scaling
    for(ltfatInt ii=0; ii<L*W; ii++)
    {
        cout[ii] *= scale;
    }

    if(kind==DSTII)
    {
        // Scale AC component
        for(ltfatInt ii=0; ii<W; ii++)
        {
            cout[(ii+1)*L-1] *= postScale;
        }
    }
}

#endif
