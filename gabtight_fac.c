#include "config.h"
#include "ltfat.h"
#include "ltfat_types.h"

LTFAT_EXTERN void
LTFAT_NAME(gabtight_fac)(const LTFAT_COMPLEX *gf, const ltfatInt L,const ltfatInt R,
                         const ltfatInt a, const ltfatInt M,
                         LTFAT_COMPLEX *gtightf)
{

    ltfatInt h_a, h_m;

    LTFAT_COMPLEX *Sf, *U, *VT, *gfwork;
    LTFAT_REAL *S;

    const LTFAT_COMPLEX zzero = (LTFAT_COMPLEX) 0.0;//{0.0, 0.0 };
    const LTFAT_COMPLEX alpha = (LTFAT_COMPLEX) 1.0+0.0*I;//{1.0, 0.0 };

    const ltfatInt N=L/a;

    const ltfatInt c=gcd(a, M,&h_a, &h_m);
    const ltfatInt p=a/c;
    const ltfatInt q=M/c;
    const ltfatInt d=N/q;

    S  = (LTFAT_REAL*)ltfat_malloc(p*sizeof(LTFAT_REAL));
    Sf = (LTFAT_COMPLEX*)ltfat_malloc(p*p*sizeof(LTFAT_COMPLEX));
    U  = (LTFAT_COMPLEX*)ltfat_malloc(p*p*sizeof(LTFAT_COMPLEX));
    VT = (LTFAT_COMPLEX*)ltfat_malloc(p*q*R*sizeof(LTFAT_COMPLEX));
    gfwork = (LTFAT_COMPLEX*)ltfat_malloc(L*R*sizeof(LTFAT_COMPLEX));

    /* Copy the contents of gf to gfwork because LAPACK overwrites
     * the input.
     */
    memcpy(gfwork,gf,sizeof(LTFAT_COMPLEX)*L*R);

    for (ltfatInt rs=0; rs<c*d; rs++)
    {
        /* Compute the thin SVD */
        LTFAT_NAME(ltfat_gesvd)(p, q*R, gfwork+rs*p*q*R, p,
                                S, U, p, VT, p);

        /* Combine U and V. */
        LTFAT_NAME(ltfat_gemm)(CblasNoTrans,CblasNoTrans,p,q*R,p,
                               &alpha,(const LTFAT_COMPLEX*)U,p,
                               (const LTFAT_COMPLEX*)VT,p,
                               &zzero,gtightf+rs*p*q*R, p);


    }

    LTFAT_SAFEFREEALL(gfwork,Sf,S,U,VT);

}


LTFAT_EXTERN void
LTFAT_NAME(gabtightreal_fac)(const LTFAT_COMPLEX *gf, const ltfatInt L, const ltfatInt R,
                             const ltfatInt a, const ltfatInt M,
                             LTFAT_COMPLEX *gtightf)
{

    ltfatInt h_a, h_m;

    LTFAT_COMPLEX *Sf, *U, *VT, *gfwork;
    LTFAT_REAL *S;

    const LTFAT_COMPLEX zzero = (LTFAT_COMPLEX) 0.0;
    const LTFAT_COMPLEX alpha = (LTFAT_COMPLEX) 1.0+0.0*I;//{1.0, 0.0 };

    const ltfatInt N=L/a;

    const ltfatInt c=gcd(a, M,&h_a, &h_m);
    const ltfatInt p=a/c;
    const ltfatInt q=M/c;
    const ltfatInt d=N/q;

    /* This is a floor operation. */
    const ltfatInt d2= d/2+1;

    S  = (LTFAT_REAL*)ltfat_malloc(p*sizeof(LTFAT_REAL));
    Sf = (LTFAT_COMPLEX*)ltfat_malloc(p*p*sizeof(LTFAT_COMPLEX));
    U  = (LTFAT_COMPLEX*)ltfat_malloc(p*p*sizeof(LTFAT_COMPLEX));
    VT = (LTFAT_COMPLEX*)ltfat_malloc(p*q*R*sizeof(LTFAT_COMPLEX));
    gfwork = (LTFAT_COMPLEX*)ltfat_malloc(L*R*sizeof(LTFAT_COMPLEX));

    /* Copy the contents of gf to gfwork because LAPACK overwrites
     * the input.
     */
    memcpy(gfwork,gf,sizeof(LTFAT_COMPLEX)*L*R);

    for (ltfatInt rs=0; rs<c*d2; rs++)
    {
        /* Compute the thin SVD */
        LTFAT_NAME(ltfat_gesvd)(p, q*R, gfwork+rs*p*q*R, p,
                                S, U, p, VT, p);

        /* Combine U and V. */
        LTFAT_NAME(ltfat_gemm)(CblasNoTrans,CblasNoTrans,p,q*R,p,
                               &alpha,(const LTFAT_COMPLEX*)U,p,
                               (const LTFAT_COMPLEX*)VT,p,
                               &zzero,gtightf+rs*p*q*R, p);
    }

    LTFAT_SAFEFREEALL(gfwork,Sf,S,U,VT);
}

