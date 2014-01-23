#include "ltfat.h"
#include "ltfat_types.h"


LTFAT_EXTERN LTFAT_NAME(dgt_shearola_plan)
LTFAT_NAME(dgt_shearola_init)(const LTFAT_COMPLEX *g, const ltfatInt gl,
                              const ltfatInt W, const ltfatInt a, const ltfatInt M,
                              const ltfatInt s0, const ltfatInt s1, const ltfatInt br,
                              const ltfatInt bl,
                              unsigned flags)
{

    LTFAT_NAME(dgt_shearola_plan) plan;

    plan.bl = bl;
    plan.gl = gl;
    plan.W  = W;

    const ltfatInt Lext    = bl+gl;
    const ltfatInt Nblocke = Lext/a;

    plan.buf  = ltfat_malloc(Lext*W*sizeof(LTFAT_COMPLEX));
    plan.gext = ltfat_malloc(Lext*sizeof(LTFAT_COMPLEX));
    plan.cbuf = ltfat_malloc(M*Nblocke*W*sizeof(LTFAT_COMPLEX));

    LTFAT_NAME(fir2long_c)(g, gl, Lext, plan.gext);

    /* Zero the last part of the buffer, it will always be zero. */
    for (ltfatInt w=0; w<W; w++)
    {
        for (ltfatInt jj=bl; jj<Lext; jj++)
        {
            plan.buf[jj+w*Lext] = (LTFAT_COMPLEX) 0.0;
        }
    }

    plan.plan =
        LTFAT_NAME(dgt_shear_init)((const LTFAT_COMPLEX*)plan.buf,
                                   (const LTFAT_COMPLEX*)plan.gext,
                                   Lext, W, a, M,
                                   s0, s1, br,
                                   plan.cbuf, flags);

    return (plan);

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola_execute)(const LTFAT_NAME(dgt_shearola_plan) plan,
                                 const LTFAT_COMPLEX *f, const ltfatInt L,
                                 LTFAT_COMPLEX *cout)

{
    const ltfatInt bl      = plan.bl;
    const ltfatInt gl      = plan.gl;
    const ltfatInt a       = plan.plan.a;
    const ltfatInt M       = plan.plan.M;
    const ltfatInt N       = L/a;
    const ltfatInt Lext    = bl+gl;
    const ltfatInt Nb      = L/bl;
    const ltfatInt b2      = gl/a/2;
    const ltfatInt Nblock  = bl/a;
    const ltfatInt Nblocke = Lext/a;
    const ltfatInt W       = plan.W;


    /* Zero the output array, as we will be adding to it */
    for (ltfatInt ii=0; ii<M*N*W; ii++)
    {
        cout[ii] = (LTFAT_COMPLEX) 0.0;
    }

    for (ltfatInt ii=0; ii<Nb; ii++)
    {
        ltfatInt s_ii;

        /* Copy to working buffer. */
        for (ltfatInt w=0; w<W; w++)
        {
            memcpy(plan.buf+Lext*w,f+ii*bl+w*L,sizeof(LTFAT_COMPLEX)*bl);
        }

        /* Execute the short DGT */
        LTFAT_NAME(dgt_shear_execute)(plan.plan);

        /* Place the results */
        for (ltfatInt w=0; w<W; w++)
        {
            /* Place large block */
            LTFAT_COMPLEX *cout_p = cout +      ii*M*Nblock+w*M*N ;
            LTFAT_COMPLEX *cbuf_p = plan.cbuf +  w*M*Nblocke;
            for (ltfatInt m=0; m<M; m++)
            {
                for (ltfatInt n=0; n<Nblock; n++)
                {
                    cout_p[m+n*M] += cbuf_p[m+n*M];
                }
            }

            /* Small block + */
            s_ii=positiverem(ii+1,Nb);
            cout_p = cout + s_ii*M*Nblock+w*M*N ;
            cbuf_p = plan.cbuf +      M*Nblock+w*M*Nblocke;
            for (ltfatInt m=0; m<M; m++)
            {
                for (ltfatInt n=0; n<b2; n++)
                {
                    cout_p[m+n*M] += cbuf_p[m+n*M];
                }
            }


            /* Small block - */
            s_ii=positiverem(ii-1,Nb)+1;
            cout_p = cout + M*(s_ii*Nblock-b2)+w*M*N ;
            cbuf_p = plan.cbuf + M*(Nblock+b2)     +w*M*Nblocke;
            for (ltfatInt m=0; m<M; m++)
            {
                for (ltfatInt n=0; n<b2; n++)
                {
                    cout_p[m+n*M] += cbuf_p[m+n*M];
                }
            }

        }

    }


}

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola_done)(LTFAT_NAME(dgt_shearola_plan) plan)
{
    LTFAT_NAME(dgt_shear_done)(plan.plan);

    /* ltfat_free(plan.cbuf); */

    LTFAT_SAFEFREEALL(plan.gext,plan.buf);

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_shearola)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                         const ltfatInt L, const ltfatInt gl, const ltfatInt W, const ltfatInt a, const ltfatInt M,
                         const ltfatInt s0, const ltfatInt s1, const ltfatInt br, const ltfatInt bl,
                         LTFAT_COMPLEX *cout)
{

    LTFAT_NAME(dgt_shearola_plan) plan = LTFAT_NAME(dgt_shearola_init)(
            g,gl,W,a,M,s0,s1,br,bl,FFTW_ESTIMATE);

    LTFAT_NAME(dgt_shearola_execute)(plan,f,L,cout);

    LTFAT_NAME(dgt_shearola_done)(plan);

}
