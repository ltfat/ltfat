#include "ltfat.h"
#include "ltfat_types.h"



// long is only "at least 32 bit"
static inline long long positiverem_long(long long a, long long b)
{
    const long long c = a % b;
    return (c < 0 ? c + b : c);
}


LTFAT_EXTERN void
LTFAT_NAME(pchirp)(const long long L, const long long n, LTFAT_COMPLEX *g)
{

    const long long LL = 2 * L;
    const long long Lponen = positiverem_long((L + 1) * n, LL);

    for (long long m = 0; m < L; m++)
    {
        const long long idx = positiverem_long(
                                  positiverem_long(Lponen * m, LL) * m, LL);

        g[m] = cexp(1.0 * I * PI * idx / L);
    }


    /* const LTFAT_REAL LL=2.0*L; */
    /* const LTFAT_REAL Lpone=L+1; */

    /* for (ltfatInt m=0;m<L;m++) */
    /* { */
    /*    //g[m] = cexp(I*PI*fmod(Lpone*n*m*m,LL)/L); */
    /*    g[m] = cexp(I*PI*fmod(fmod(fmod(Lpone*n,LL)*m,LL)*m,LL)/L); */
    /* } */

}


LTFAT_EXTERN LTFAT_NAME(dgt_shear_plan)
LTFAT_NAME(dgt_shear_init)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                           const ltfatInt L, const ltfatInt W, const ltfatInt a, const ltfatInt M,
                           const ltfatInt s0, const ltfatInt s1, const ltfatInt br,
                           LTFAT_COMPLEX *cout,
                           unsigned flags)
{
    LTFAT_NAME(dgt_shear_plan) plan;

    plan.a = a;
    plan.M = M;
    plan.L = L;
    plan.W = W;

    plan.s0 = s0;
    plan.s1 = s1;
    plan.br = br;

    const ltfatInt b = L / M;
    const ltfatInt N = L / a;

    const ltfatInt ar = a * b / br;
    const ltfatInt Mr = L / br;
    const ltfatInt Nr = L / ar;

    plan.f     = (LTFAT_COMPLEX *)f;
    plan.fwork = (LTFAT_COMPLEX *)f;
    plan.gwork = (LTFAT_COMPLEX *)g;
    plan.cout  = cout;

    plan.c_rect = ltfat_malloc(M * N * W * sizeof(LTFAT_COMPLEX));

    LTFAT_COMPLEX *f_before_fft = (LTFAT_COMPLEX *)f;
    LTFAT_COMPLEX *g_before_fft = (LTFAT_COMPLEX *)g;

    if ((s0 != 0) || (s1 != 0))
    {
        plan.fwork = ltfat_malloc(L * W * sizeof(LTFAT_COMPLEX));
        plan.gwork = ltfat_malloc(L * sizeof(LTFAT_COMPLEX));
    }


    if (s1)
    {
        plan.p1 = ltfat_malloc(L * sizeof(LTFAT_COMPLEX));

        LTFAT_NAME(pchirp)(L, s1, plan.p1);

        for (ltfatInt l = 0; l < L; l++)
        {
            plan.gwork[l] = g[l] * plan.p1[l];
        }

        f_before_fft = plan.fwork;
        g_before_fft = plan.gwork;

    }

    if (s0 == 0)
    {

        /* Call the rectangular computation in the time domain */
        /* LTFAT_NAME(dgt_long)(plan.fwork,plan.gwork,L,W,ar,Mr,plan.c_rect); */

        plan.rect_plan = LTFAT_NAME(dgt_long_init)(plan.fwork, plan.gwork,
                         L, W, ar, Mr, plan.c_rect, 0, flags);
    }
    else
    {

        /* Allocate memory and compute the pchirp */
        plan.p0 = ltfat_malloc(L * sizeof(LTFAT_COMPLEX));
        LTFAT_NAME(pchirp)(L, -s0, plan.p0);

        /* if data has already been copied to the working arrays, use
         * inline FFTs. Otherwise, if this is the first time they are
         * being used, do the copying using the fft. */

        // Downcasting to int
        int Lint = (int) L;
        plan.f_plan = LTFAT_FFTW(plan_many_dft)(1, &Lint, W,
                                                f_before_fft, NULL, 1, L,
                                                plan.fwork, NULL, 1, L,
                                                FFTW_FORWARD, flags);

        plan.g_plan = LTFAT_FFTW(plan_dft_1d)(L, g_before_fft, plan.gwork, FFTW_FORWARD, flags);

        /* Execute the FFTs */
        LTFAT_FFTW(execute)(plan.g_plan);

        /* Multiply g by the chirp and scale by 1/L */
        for (ltfatInt l = 0; l < L; l++)
        {
            plan.gwork[l] = plan.gwork[l] * plan.p0[l] / L;
        }

        /* Call the rectangular computation in the frequency domain*/
        /* LTFAT_NAME(dgt_long)(plan.fwork,plan.gwork,L,W,br,Nr,plan.c_rect); */
        /* Call the rectangular computation in the frequency domain*/
        plan.rect_plan = LTFAT_NAME(dgt_long_init)(plan.fwork, plan.gwork, L, W,
                         br, Nr, plan.c_rect, 0, flags);

    }

    plan.finalmod = ltfat_malloc(2 * N * sizeof(LTFAT_COMPLEX));

    for (ltfatInt n = 0; n < 2 * N; n++)
    {
        plan.finalmod[n] = cexp(PI * I * n / N);
    }

    return plan;

}

LTFAT_EXTERN void
LTFAT_NAME(dgt_shear_execute)(const LTFAT_NAME(dgt_shear_plan) plan)
{

    const ltfatInt a = plan.a;
    const ltfatInt M = plan.M;
    const ltfatInt L = plan.L;

    const ltfatInt b = plan.L / plan.M;
    const ltfatInt N = plan.L / plan.a;
    const long s0 = plan.s0;
    const long s1 = plan.s1;

    const ltfatInt ar = plan.a * b / plan.br;
    const ltfatInt Mr = plan.L / plan.br;
    const ltfatInt Nr = plan.L / ar;


    if (s1)
    {
        for (ltfatInt w = 0; w < plan.W; w++)
        {
            for (ltfatInt l = 0; l < plan.L; l++)
            {
                plan.fwork[l + w * plan.L] = plan.f[l + w * plan.L] * plan.p1[l];
            }
        }

    }


    if (s0 == 0)
    {

        const ltfatInt twoN = 2 * N;

        /* In this case, cc1=1 */

        const long cc3 = positiverem_long(s1 * (L + 1), twoN);

        const long tmp1 = positiverem_long(cc3 * a, twoN);

        LTFAT_NAME(dgt_long_execute)(plan.rect_plan);

        for (ltfatInt k = 0; k < N; k++)
        {
            const ltfatInt phsidx = positiverem_long((tmp1 * k) % twoN * k, twoN);
            const long part1 = positiverem_long(-s1 * k * a, L);
            for (ltfatInt m = 0; m < M; m++)
            {
                /* The line below has a hidden floor operation when dividing with the last b */
                const ltfatInt idx2 = ((part1 + b * m) % L) / b;

                const ltfatInt inidx  =    m + k * M;
                const ltfatInt outidx = idx2 + k * M;
                for (ltfatInt w = 0; w < plan.W; w++)
                {
                    plan.cout[outidx + w * M * N] = plan.c_rect[inidx + w * M * N] * plan.finalmod[phsidx];
                }
            }
        }


    }
    else
    {

        const ltfatInt twoN = 2 * N;
        const long cc1 = ar / a;
        const long cc2 = positiverem_long(-s0 * plan.br / a, twoN);
        const long cc3 = positiverem_long(a * s1 * (L + 1), twoN);
        const long cc4 = positiverem_long(cc2 * plan.br * (L + 1), twoN);
        const long cc5 = positiverem_long(2 * cc1 * plan.br, twoN);
        const long cc6 = positiverem_long((s0 * s1 + 1) * plan.br, L);

        LTFAT_FFTW(execute)(plan.f_plan);

        for (ltfatInt w = 0; w < plan.W; w++)
        {
            for (ltfatInt l = 0; l < plan.L; l++)
            {

                plan.fwork[l + w * plan.L] = plan.fwork[l + w * plan.L] * plan.p0[l];
            }
        }

        LTFAT_NAME(dgt_long_execute)(plan.rect_plan);

        for (ltfatInt k = 0; k < Nr; k++)
        {
            const long part1 = positiverem_long(-s1 * k * ar, L);
            for (ltfatInt m = 0; m < Mr; m++)
            {
                const long sq1 = k * cc1 + cc2 * m;

                const ltfatInt phsidx = positiverem_long(
                                            (cc3 * sq1 * sq1) % twoN - (m * (cc4 * m + k * cc5)) % twoN, twoN);

                /* The line below has a hidden floor operation when dividing with the last b */
                const ltfatInt idx2 = ((part1 + cc6 * m) % L) / b;

                const ltfatInt inidx  = positiverem(-k, Nr) + m * Nr;
                const ltfatInt outidx = idx2 + (sq1 % N) * M;
                for (ltfatInt w = 0; w < plan.W; w++)
                {
                    plan.cout[outidx + w * M * N] = plan.c_rect[inidx + w * M * N] * plan.finalmod[phsidx];

                }
            }
        }

    }
}


LTFAT_EXTERN void
LTFAT_NAME(dgt_shear_done)(LTFAT_NAME(dgt_shear_plan) plan)
{
    LTFAT_NAME(dgt_long_done)(plan.rect_plan);
    LTFAT_SAFEFREEALL(plan.finalmod, plan.c_rect, plan.fwork, plan.gwork, plan.p0, plan.p1);
}


LTFAT_EXTERN void
LTFAT_NAME(dgt_shear)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                      const ltfatInt L, const ltfatInt W, const ltfatInt a, const ltfatInt M,
                      const ltfatInt s0, const ltfatInt s1, const ltfatInt br,
                      LTFAT_COMPLEX *cout)
{

    LTFAT_NAME(dgt_shear_plan) plan = LTFAT_NAME(dgt_shear_init)(
                                          f, g, L, W, a, M, s0, s1, br, cout, FFTW_ESTIMATE);

    LTFAT_NAME(dgt_shear_execute)(plan);

    LTFAT_NAME(dgt_shear_done)(plan);

}




