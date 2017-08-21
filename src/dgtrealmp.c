#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"
#include "dgtrealmp_private.h"


LTFAT_API LTFAT_NAME(dgtreal_plan)**
LTFAT_NAME(dgtrealmp_getdgtrealplan)(LTFAT_NAME(dgtrealmp_plan)* p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);

    return p->dgtplans;
error:
    return NULL;

}

LTFAT_API LTFAT_COMPLEX**
LTFAT_NAME(dgtrealmp_getresidualcoef)(LTFAT_NAME(dgtrealmp_plan)* p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);

    return p->iterstate->c;
error:
    return NULL;

}

LTFAT_API ltfat_dgtrealmp_params*
LTFAT_NAME(dgtrealmp_getparams)(LTFAT_NAME(dgtrealmp_plan)* p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);

    return p->params;
error:
    return NULL;
}

LTFAT_API int
LTFAT_NAME(dgtrealmp_init)(
    const LTFAT_REAL* g[], ltfat_int gl[], ltfat_int L, ltfat_int P, ltfat_int a[],
    ltfat_int M[], ltfat_dgtrealmp_params* params,
    LTFAT_NAME(dgtrealmp_plan)** pout)
{
    int status = LTFATERR_FAILED;
    const LTFAT_REAL* gtmp[2]; ltfat_int gltmp[2]; ltfat_int atmp[2];
    ltfat_int Mtmp[2];
    LTFAT_NAME(dgtrealmp_plan)* p = NULL;

    ltfat_int nextL = ltfat_dgtlengthmulti(L, P, a, M);

    DEBUG("L=%td, nextL=%td, P=%td, gl[0]=%td, a[0]=%td, M[0]=%td",
          L, nextL, P, gl[0], a[0], M[0]);

    CHECKNULL( pout );
    CHECK(LTFATERR_BADTRALEN, L == nextL,
          "Next compatible transform length is %d (passed %d).", nextL, L);

    CHECKMEM( p = LTFAT_NEW( LTFAT_NAME(dgtrealmp_plan)) );
    CHECKMEM( p->params = ltfat_dgtrealmp_params_allocdef() );

    if (params)
        memcpy( p->params, params, sizeof * p->params);

    if ( p->params->maxatoms == 0 )
        p->params->maxatoms =  (size_t) ( 0.1 * L);

    if ( p->params->maxit == 0 )
        p->params->maxit = 2 * p->params->maxatoms;

    p->params->initwasrun = 1;

    DEBUG("maxit=%zu, maxatoms=%zu, iterstep=%zu, errtoldb=%f",
          p->params->maxit, p->params->maxatoms, p->params->iterstep,
          p->params->errtoldb);

    CHECKMEM( p->dgtplans  = LTFAT_NEWARRAY( LTFAT_NAME(dgtreal_plan)*, P) );
    CHECKMEM( p->gramkerns = LTFAT_NEWARRAY( LTFAT_NAME(kerns)*, P * P) );
    CHECKMEM( p->a  = LTFAT_NEWARRAY( ltfat_int, P));
    CHECKMEM( p->M  = LTFAT_NEWARRAY( ltfat_int, P));
    CHECKMEM( p->M2 = LTFAT_NEWARRAY( ltfat_int, P));
    CHECKMEM( p->N  = LTFAT_NEWARRAY( ltfat_int, P));

    for (ltfat_int k = 0; k < P; k++)
    {
        p->a[k] = a[k]; p->M[k] = M[k];
        p->M2[k] = M[k] / 2 + 1; p->N[k] = L / a[k];
    }

    p->P = P; p->L = L;

    for (ltfat_int k = 0; k < P; k++)
    {
        CHECKSTATUS(
            LTFAT_NAME(dgtreal_init_gen)(g[k], gl[k], g[k], gl[k], L, 1, a[k], M[k],
                                         NULL, NULL, NULL, &p->dgtplans[k] ),
            "dgtreal_init failed" );
    }

    for (ltfat_int k1 = 0; k1 < P; k1++)
    {
        for (ltfat_int k2 = 0; k2 < P; k2++)
        {
            gtmp[0] = g[k1]; gtmp[1] = g[k2]; atmp[0] = a[k1]; atmp[1] = a[k2];
            Mtmp[0] = M[k1]; Mtmp[1] = M[k2]; gltmp[0] = gl[k1]; gltmp[1] = gl[k2];

            CHECKSTATUS( LTFAT_NAME(dgtrealmp_kernel_init)( gtmp, gltmp,
                         atmp, Mtmp, L, p->params->kernrelthr, &p->gramkerns[k1 + k2 * P]),
                         "dgtrealmp_kernel_init failed");
        }
    }

    CHECKSTATUS( LTFAT_NAME(dgtrealmpiter_init)(a, M, P, L, &p->iterstate),
                 "dgtrealmpiter_init failed" );

    *pout = p;
    return LTFATERR_SUCCESS;
error:
    if (p) LTFAT_NAME(dgtrealmp_done)(&p);
    *pout = NULL;
    return status;
}


LTFAT_API int
LTFAT_NAME(dgtrealmp_reset)(LTFAT_NAME(dgtrealmp_plan)* p, const LTFAT_REAL* f)
{
    int status = LTFATERR_SUCCESS;

    CHECKNULL(p); CHECKNULL(f);

    p->iterstate->currit = 0;
    p->iterstate->curratoms = 0;
    p->iterstate->err = 0.0;

    /* DEBUG("L=%td, M=%td, M2=%td, N=%td,P=%td", p->L, p->M[0], p->M2[0], p->N[0], */
    /*       p->P); */

    for (ltfat_int l = 0; l < p->L; l++)
        p->iterstate->err += f[l] * f[l];

    p->iterstate->fnorm2 = p->iterstate->err;
    p->params->errtoladj = pow(10.0,
                               p->params->errtoldb / 10.0) * p->iterstate->fnorm2;

    for (ltfat_int k = 0; k < p->P; k++)
    {
        LTFAT_REAL* sEl = p->iterstate->s[k];
        LTFAT_COMPLEX* cEl = p->iterstate->c[k];

        CHECKSTATUS(
            LTFAT_NAME(dgtreal_execute_ana_newarray)(p->dgtplans[k], f, cEl),
            "dgtreal_execute failed");

        LTFAT_NAME_COMPLEX(dgtreal2dgt)(cEl, p->M[k], p->N[k], cEl);

        for (size_t l = 0; l < (size_t) ( p->M[k] * p->N[k]); l++ )
            sEl[l] = ltfat_abs(cEl[l]);

        LTFAT_NAME_REAL(findmaxincols)(sEl, p->M[k], p->M2[k], p->N[k],
                                       p->iterstate->maxcols[k],
                                       p->iterstate->maxcolspos[k]);

        memset( p->iterstate->suppindCount[k], 0 ,
                p->M2[k] * p->N[k] * sizeof * p->iterstate->suppindCount[k] );
    }

error:
    return status;
}

#define LTFAT_DGTREALMP_SUBSTRACTKERNEL(ctmp)\
for (ltfat_int nidx = nstart, knidx = knstart; knidx < currkern->size.width;\
     nidx = ltfat_positiverem(nidx + 1, N), knidx+=astep )\
{\
\
    LTFAT_COMPLEX* currcCol = currc + nidx * M;\
    LTFAT_REAL*    currsCol = currs + nidx * M;\
    LTFAT_COMPLEX* kcurrCol = currkernvals + knidx * currkern->size.height;\
\
    for (ltfat_int midx = mstart, kmidx = kmstart ; kmidx < currkern->size.height;\
         midx = ltfat_positiverem(midx + 1, M), kmidx+=Mstep )\
    {\
        currcCol[midx] -= ctmp * kcurrCol[kmidx];\
        currsCol[midx] = ltfat_abs(currcCol[midx]);\
    }\
}

#define LTFAT_DGTREALMP_COPYTOCONJ \
for (ltfat_int nidx = nstart, knidx = knstart; knidx < currkern->size.width;\
     nidx = ltfat_positiverem(nidx + 1, N), knidx+=astep )\
{\
\
    LTFAT_COMPLEX*     currcCol = currc + nidx * M + mstart;\
    LTFAT_COMPLEX* currcconjCol = currc + nidx * M + mconjstop;\
\
    for (ltfat_int kmidx = kmstart; kmidx < currkern->size.height;\
            kmidx+=Mstep )\
          *currcconjCol-- = conj(*currcCol++);\
}

#define LTFAT_DGTREALMP_UPDATEMAX(wtmp)\
for (ltfat_int nidx = nstart, knidx = knstart; knidx < currkern->size.width;\
     nidx = ltfat_positiverem(nidx + 1, N), knidx+=astep )\
{\
    LTFAT_NAME_REAL(findmaxinarray)(currs + nidx * M, M2,\
                                    s->maxcols[wtmp] + nidx,\
                                    s->maxcolspos[wtmp] + nidx);\
}


LTFAT_API int
LTFAT_NAME(dgtrealmp_execute_niters)(LTFAT_NAME(dgtrealmp_plan)* p,
                                     ltfat_int itno, LTFAT_COMPLEX** cout)
{
    int status = LTFAT_DGTREALMP_STATUS_WILLCONTINUE;

    LTFAT_NAME(dgtrealmpiter_state)* s = p->iterstate;

    for (ltfat_int iter = 0; iter < itno; iter++)
    {
        s->currit++;

        // 1) Find max m and n
        ltfat_int n = 0, m = 0, mconj = 0, wIdx = 0, M, M2, N;
        LTFAT_COMPLEX cval = 0.0, cval2 = 0.0;
        LTFAT_REAL cvalabs = 0.0, val = 0.0;

        for (ltfat_int k = 0; k < s->P; k++)
        {
            LTFAT_REAL valTmp;
            ltfat_int nTmp;
            LTFAT_NAME_REAL(findmaxinarray)(s->maxcols[k], p->N[k], &valTmp, &nTmp);

            if ( valTmp > val )
            {
                val = valTmp;
                m = s->maxcolspos[k][nTmp];
                n = nTmp;
                wIdx = k;
            }
        }

        LTFAT_NAME(kerns)* currkern = p->gramkerns[wIdx + s->P * wIdx];

        /* DEBUG("m=%td, n=%td, wIdx=%td, N=%td, ", m, n, wIdx, p->N[wIdx]); */

        M  = p->M[wIdx]; M2 = p->M2[wIdx]; N  = p->N[wIdx]; mconj = M - m;
        int uniquenyquest = M % 2 == 0;
        int do_conj = !( m == 0 || (m == M2 - 1 && uniquenyquest));
        // true if the support of the kernel is completely contained in the
        // positive frequencies
        int do_copy = (m - currkern->mid.hmid) > 0 &&
                      (m + currkern->size.height - currkern->mid.hmid - 1) < M2 - 1 - uniquenyquest;

        if ( s->suppindCount[wIdx][m + M2 * n] == 0 )
            s->curratoms++;

        s->suppindCount[wIdx][m + M2 * n] += 1;

        // 2) The maximum coeffcicient
        cval = s->c[wIdx][m + M * n];

        // Update the solution
        cout[wIdx][m + M2 * n] += cval;

        // Update error
        cvalabs = ltfat_abs(cval);

        /* DEBUG("cvalreal=%.4f, cvalimag=%.4f, abs=%.4f", ltfat_real(cval), */
        /*       ltfat_imag(cval), cvalabs); */

        s->err -= cvalabs * cvalabs;
        if (do_conj)
            s->err -= cvalabs * cvalabs;

        if (s->err < 0)
        {
            status = LTFAT_DGTREALMP_STATUS_STALLED;
            break;
        }

        LTFAT_COMPLEX* currkernvals = currkern->kval[n % currkern->kNo];

        ltfat_int Mstep   = 1;
        ltfat_int astep   = 1;
        ltfat_int mstart  = ltfat_positiverem(m - currkern->mid.hmid, M);
        ltfat_int nstart  = ltfat_positiverem(n - currkern->mid.wmid, N);
        ltfat_int kmstart = 0;
        ltfat_int knstart = 0;

        LTFAT_COMPLEX* currc = s->c[wIdx];
        LTFAT_REAL* currs    = s->s[wIdx];

        LTFAT_DGTREALMP_SUBSTRACTKERNEL(cval)

        if (do_conj)
        {
            if (!do_copy)
            {
                cval2 = currc[mconj + M * n];
                mstart = ltfat_positiverem(mconj - currkern->mid.hmid, M);

                LTFAT_DGTREALMP_SUBSTRACTKERNEL(cval2)
            }
            else
            {
                ltfat_int mconjstop = mconj +
                                      (currkern->size.height - currkern->mid.hmid - 1);

                LTFAT_DGTREALMP_COPYTOCONJ
            }
        }

        LTFAT_DGTREALMP_UPDATEMAX(wIdx)

        for (ltfat_int secondwIdx = 0; secondwIdx < s->P; secondwIdx++)
        {
            if (secondwIdx == wIdx)
                continue;

            M  = p->M[secondwIdx]; M2 = p->M2[secondwIdx]; N = p->N[secondwIdx];


            double Mrat = p->M[wIdx] / M;
            double arat = p->a[secondwIdx] / p->a[wIdx];

            Mstep = Mrat > 1 ? (ltfat_int) Mrat : 1;
            astep = arat > 1 ? (ltfat_int) arat : 1;

            ltfat_int nsecond    = (ltfat_int) ltfat_round(n / arat);
            ltfat_int nsecondoff = n - (ltfat_int)(nsecond * arat);
            ltfat_int msecond    = (ltfat_int) ltfat_round( m / Mrat);
            ltfat_int msecondoff = m - (ltfat_int)(msecond * Mrat);

            currkern = p->gramkerns[wIdx + s->P * secondwIdx];
            currkernvals = currkern->kval[n % currkern->kNo];

            ltfat_int Mremprev = (currkern->mid.hmid - msecondoff) / Mstep;
            ltfat_int aremprev = (currkern->mid.wmid - nsecondoff) / astep;

            mstart = ltfat_positiverem(msecond - Mremprev, M);
            nstart = ltfat_positiverem(nsecond - aremprev, N);

            knstart = currkern->mid.wmid - nsecondoff - aremprev * astep;
            kmstart = currkern->mid.hmid - msecondoff - Mremprev * Mstep;

            currc = s->c[secondwIdx];
            currs = s->s[secondwIdx];

            LTFAT_DGTREALMP_SUBSTRACTKERNEL(cval)

            if (do_conj)
            {
                ltfat_int msecondconj = (ltfat_int) ltfat_round( mconj / Mrat);
                ltfat_int msecondoffconj = mconj - (ltfat_int)(msecondconj * Mrat);

                Mremprev = (currkern->mid.hmid - msecondoffconj) / Mstep;
                mstart = ltfat_positiverem(msecond - Mremprev, M);
                kmstart = currkern->mid.hmid - msecondoffconj - Mremprev * Mstep;

                LTFAT_DGTREALMP_SUBSTRACTKERNEL(cval2)
            }

            LTFAT_DGTREALMP_UPDATEMAX(secondwIdx)
        }

        if (s->err <= p->params->errtoladj)
        {
            status = LTFAT_DGTREALMP_STATUS_TOLREACHED;
            break;
        }

        if (s->curratoms >= p->params->maxatoms)
        {
            status = LTFAT_DGTREALMP_STATUS_MAXATOMS;
            break;
        }

        if (s->currit >= p->params->maxit)
        {
            status = LTFAT_DGTREALMP_STATUS_MAXITER;
            break;
        }
    }

    return status;
}

#undef LTFAT_DGTREALMP_SUBSTRACTKERNEL
#undef LTFAT_DGTREALMP_COPYTOCONJ
#undef LTFAT_DGTREALMP_UPDATEMAX

LTFAT_API int
LTFAT_NAME(dgtrealmp_execute)(LTFAT_NAME(dgtrealmp_plan)* p,
                              const LTFAT_REAL* f,
                              LTFAT_COMPLEX** cout,
                              LTFAT_REAL* fout)
{
    int status = LTFATERR_SUCCESS;
    int status2 = LTFATERR_SUCCESS;
    size_t iterstep = p->params->iterstep;
    LTFAT_REAL* ftmp = NULL;

    CHECKNULL(p); CHECKNULL(f); CHECKNULL(cout); CHECKNULL(fout);

    CHECKSTATUS( LTFAT_NAME(dgtrealmp_reset)( p, f),
                 "dgtrealmp_reset failed" );

    for (ltfat_int k = 0; k < p->P; k++)
        memset(cout[k], 0, p->M2[k] * p->N[k] * sizeof * cout[k]);

    while ( LTFAT_DGTREALMP_STATUS_WILLCONTINUE ==
            ( status2 =
                  LTFAT_NAME(dgtrealmp_execute_niters)( p, iterstep, cout)))
    {
        if (p->params->verbose)
        {
            printf("Error %.2f \n", 10.0 * log10( p->iterstate->err /
                                                  p->iterstate->fnorm2 ) );
        }
    }

    if (status2 < 0)
        CHECKSTATUS(status2, "Iterations exited with an error");
    else
    {
        // TO DO: Finished sucessfully with some exit code
    }

    CHECKSTATUS( LTFAT_NAME(dgtreal_execute_syn_newarray)(
                     p->dgtplans[0], cout[0], fout),
                 "dgtreal_execute_syn failed");

    if (p->P > 1)
    {
        CHECKMEM( ftmp = LTFAT_NAME_REAL(malloc)(p->L) );

        for (ltfat_int k = 1; k < p->P; k++)
        {
            LTFAT_NAME(dgtreal_execute_syn_newarray)(p->dgtplans[k], cout[k], ftmp);

            for (ltfat_int l = 0; l < p->L; l++)
                fout[l] += ftmp[l];
        }

        ltfat_free(ftmp);
    }

    return status2;
error:
    ltfat_safefree(ftmp);
    return status;
}


LTFAT_API int
LTFAT_NAME(dgtrealmp_done)( LTFAT_NAME(dgtrealmp_plan)** p)
{
    LTFAT_NAME(dgtrealmp_plan)* pp;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(*p);
    pp = *p;

    ltfat_safefree(pp->a);
    ltfat_safefree(pp->M);
    ltfat_safefree(pp->M2);
    ltfat_safefree(pp->N);

    if (pp->params)
        ltfat_dgtrealmp_params_free(pp->params);

    if (pp->dgtplans)
    {
        for (ltfat_int k = 0; k < pp->P; k++)
            LTFAT_NAME(dgtreal_done)(&pp->dgtplans[k]);

        ltfat_free(pp->dgtplans);
        pp->dgtplans = NULL;
    }

    if (pp->gramkerns)
    {
        for (ltfat_int k = 0; k < pp->P * pp->P; k++)
            LTFAT_NAME(dgtrealmp_kernel_done)( &pp->gramkerns[k]);

        ltfat_free(pp->gramkerns);
        pp->gramkerns = NULL;
    }

    if (pp->iterstate)
        LTFAT_NAME(dgtrealmpiter_done)(&pp->iterstate);

    ltfat_free(pp);
    *p = NULL;
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(dgtrealmpiter_init)(ltfat_int a[], ltfat_int M[], ltfat_int P,
                               ltfat_int L, LTFAT_NAME(dgtrealmpiter_state)** state)
{
    LTFAT_NAME(dgtrealmpiter_state)* s = NULL;
    int status = LTFATERR_FAILED;

    CHECKNULL( state );
    CHECKMEM( s =    LTFAT_NEW( LTFAT_NAME(dgtrealmpiter_state)) );
    CHECKMEM( s->s = LTFAT_NEWARRAY(LTFAT_REAL*,    P));
    CHECKMEM( s->c = LTFAT_NEWARRAY(LTFAT_COMPLEX*, P));
    CHECKMEM( s->suppindCount = LTFAT_NEWARRAY(int*, P));
    CHECKMEM( s->maxcols =      LTFAT_NEWARRAY(LTFAT_REAL*, P));
    CHECKMEM( s->maxcolspos =   LTFAT_NEWARRAY(ltfat_int*, P));
    s->P = P;

    for (ltfat_int p = 0; p < P; p++)
    {
        ltfat_int N = L / a[p];
        CHECKMEM( s->s[p] = LTFAT_NAME_REAL(malloc)(N * M[p]) );
        CHECKMEM( s->c[p] = LTFAT_NAME_COMPLEX(malloc)(N * M[p]) );
        CHECKMEM( s->suppindCount[p] = LTFAT_NEWARRAY(int, N * (M[p] / 2 + 1)) );
        CHECKMEM( s->maxcols[p] = LTFAT_NAME_REAL(malloc)(N) );
        CHECKMEM( s->maxcolspos[p] = LTFAT_NEWARRAY(ltfat_int, N) );
    }

    *state = s;
    return LTFATERR_SUCCESS;
error:
    if (s) LTFAT_NAME(dgtrealmpiter_done)(&s);
    *state = NULL;
    return status;
}

LTFAT_API int
LTFAT_NAME(dgtrealmpiter_done)(LTFAT_NAME(dgtrealmpiter_state)** state)
{
    LTFAT_NAME(dgtrealmpiter_state)* s = NULL;
    int status = LTFATERR_SUCCESS;
    CHECKNULL(state); CHECKNULL(*state);
    s = *state;

    if (s->s)
    {
        for (ltfat_int p = 0; p < s->P; p++)
            ltfat_safefree(s->s[p]);

        ltfat_free(s->s);
    }

    if (s->c)
    {
        for (ltfat_int p = 0; p < s->P; p++)
            ltfat_safefree(s->c[p]);

        ltfat_free(s->c);
    }

    if (s->suppindCount)
    {
        for (ltfat_int p = 0; p < s->P; p++)
            ltfat_safefree(s->suppindCount[p]);

        ltfat_free(s->suppindCount);
    }

    if (s->maxcols)
    {
        for (ltfat_int p = 0; p < s->P; p++)
            ltfat_safefree(s->maxcols[p]);

        ltfat_free(s->maxcols);
    }

    if (s->maxcolspos)
    {
        for (ltfat_int p = 0; p < s->P; p++)
            ltfat_safefree(s->maxcolspos[p]);

        ltfat_free(s->maxcolspos);
    }

    ltfat_free(s);
    *state = NULL;
error:
    return status;
}

int
LTFAT_NAME(dgtrealmp_kernel_init)( const LTFAT_REAL* g[], ltfat_int gl[],
                                   ltfat_int a[], ltfat_int M[],
                                   ltfat_int L, LTFAT_REAL reltol,
                                   LTFAT_NAME(kerns)** pout)
{
    double arat;
    ltfat_int modNo, amin, Mmax, lefttail0, righttail0, lefttail1, righttail1,
              Lshort, Nshort;
    LTFAT_REAL* g0tmp = NULL, *g1tmp = NULL;
    LTFAT_COMPLEX* kernlarge = NULL;
    LTFAT_NAME(kerns)* ktmp = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKMEM( ktmp = LTFAT_NEW(LTFAT_NAME(kerns)) );

    amin = a[0] > a[1] ? a[0] : a[1];
    Mmax = M[0] > M[1] ? M[0] : M[1];

    LTFAT_NAME(dgtrealmp_essentialsupport)(g[0], gl[0], 1e-6, &lefttail0,
                                           &righttail0);
    LTFAT_NAME(dgtrealmp_essentialsupport)(g[1], gl[1], 1e-6, &lefttail1,
                                           &righttail1);

    DEBUG("lefttail=%td, righttail=%td", lefttail0, righttail0);

    Lshort =
        (lefttail0 > righttail0 ? 2 * lefttail0 : 2 * righttail0) +
        (lefttail1 > righttail1 ? 2 * lefttail1 : 2 * righttail1);

    Lshort = ltfat_dgtlength(Lshort > L ? L : Lshort, amin, Mmax);

    DEBUG("Lshort=%td", Lshort);

    Nshort = Lshort / amin;

    CHECKMEM(g0tmp     = LTFAT_NAME_REAL(malloc)(Lshort));
    CHECKMEM(g1tmp     = LTFAT_NAME_REAL(malloc)(Lshort));
    CHECKMEM(kernlarge = LTFAT_NAME_COMPLEX(malloc)(Mmax * Nshort));

    LTFAT_NAME(middlepad)(g[0], gl[0], LTFAT_WHOLEPOINT, Lshort, g0tmp);
    LTFAT_NAME(middlepad)(g[1], gl[1], LTFAT_WHOLEPOINT, Lshort, g1tmp);

    LTFAT_NAME_REAL(dgt_long)(g0tmp, g1tmp, Lshort, 1, amin, Mmax,
                              LTFAT_FREQINV, kernlarge);

    LTFAT_NAME(dgtrealmp_kernel_findsmallsize)(kernlarge, Mmax, Nshort,
            reltol, &ktmp->size, &ktmp->mid);

    /* for (ltfat_int m = 0; m < Mmax; m++) */
    /* { */
    /*     printf("\n"); */
    /*     for (ltfat_int n = 0; n < Nshort; n++) */
    /*         printf("%2.2e ", ltfat_abs(kernlarge[m + Mmax * n]) ); */
    /* } */
    /*  */
    /* printf("\n"); */

    LTFAT_NAME_COMPLEX(circshift2)(kernlarge, Mmax, Nshort,
                                   ktmp->mid.hmid, ktmp->mid.wmid, kernlarge);

    /* for (ltfat_int m = 0; m < Mmax; m++) */
    /* { */
    /*     printf("\n"); */
    /*     for (ltfat_int n = 0; n < Nshort; n++) */
    /*         printf("%2.2e ", ltfat_abs(kernlarge[m + Mmax * n]) ); */
    /* } */
    /*  */
    /* printf("\n"); */

    if (a[0] > a[1])
    {
        modNo = ltfat_lcm(amin, Mmax) / a[0];
        arat = a[0] / a[1];
    }
    else
    {
        modNo = ltfat_lcm(amin, Mmax) / amin;
        arat = 1;
    }

    ktmp->kNo = modNo;


    DEBUG("h=%td,w=%td,hmid=%td,wmid=%td,kno=%td",
          ktmp->size.height, ktmp->size.width, ktmp->mid.hmid, ktmp->mid.wmid, modNo );

    CHECKMEM( ktmp->kval = LTFAT_NEWARRAY(LTFAT_COMPLEX*, ktmp->kNo));

    for (ltfat_int k = 0; k < ktmp->kNo; k++)
        CHECKMEM( ktmp->kval[k] =
                      LTFAT_NAME_COMPLEX(malloc)( ktmp->size.height * ktmp->size.width));

    // Copy the zero-th kernel
    for (ltfat_int n = 0; n < ktmp->size.width; n++)
        memcpy(ktmp->kval[0] + n * ktmp->size.height,
               kernlarge + n * Mmax,
               ktmp->size.height * sizeof * kernlarge);


    /* for (ltfat_int m = 0; m < ktmp->size.height; m++) */
    /* { */
    /*     printf("\n"); */
    /*     for (ltfat_int n = 0; n < ktmp->size.width; n++) */
    /*         printf("%2.2e ", ltfat_abs(ktmp->kval[0][m + ktmp->size.height * n]) ); */
    /* } */

    printf("\n");

    // Compute kernel modulations
    for (ltfat_int n = 1; n < ktmp->kNo; n++)
        LTFAT_NAME(dgtrealmp_kernel_modfi)(ktmp->kval[0], ktmp->size, ktmp->mid,
                                           arat * n, amin, Mmax, ktmp->kval[n]);

    *pout = ktmp;
    LTFAT_SAFEFREEALL(g0tmp, g1tmp, kernlarge);
    return status;
error:
    LTFAT_SAFEFREEALL(g0tmp, g1tmp, kernlarge);
    if (ktmp) LTFAT_NAME(dgtrealmp_kernel_done)(&ktmp);
    return status;
}

int
LTFAT_NAME(dgtrealmp_kernel_done)(LTFAT_NAME(kerns)** k)
{
    LTFAT_NAME(kerns)* kk;
    int status = LTFATERR_SUCCESS;

    CHECKNULL(k); CHECKNULL(*k);
    kk = *k;

    for (ltfat_int kIdx = 0; kIdx < kk->kNo; kIdx++)
        ltfat_safefree( kk->kval[kIdx] );

    ltfat_safefree( kk->kval );
    ltfat_free(kk);
    *k = NULL;

error:
    return status;
}

int
LTFAT_NAME(dgtrealmp_kernel_modfi)(LTFAT_COMPLEX* kfirst, ksize size,
                                   kanchor mid, ltfat_int n, ltfat_int a, ltfat_int M,
                                   LTFAT_COMPLEX* kmod)
{
    for (ltfat_int nn = 0; nn < size.width; nn++)
    {
        const LTFAT_COMPLEX* kfirstCol = kfirst + nn * size.height;
        LTFAT_COMPLEX* kmodCol   = kmod + nn * size.height;

        for (ltfat_int m = 0; m < size.height; m++ )
        {
            ltfat_int xval = m - mid.hmid;
            LTFAT_REAL exparg = -2.0 * M_PI * n * xval * a / ((double) M);
            kmodCol[m] = exp( I * exparg) * kfirstCol[m];
        }
    }

    return 0;
}


int
LTFAT_NAME(dgtrealmp_essentialsupport)(const LTFAT_REAL g[], ltfat_int gl,
                                       LTFAT_REAL reltol,
                                       ltfat_int* lefttail, ltfat_int* righttail)
{
    ltfat_int gl2 = gl / 2 + 1;
    LTFAT_REAL gthr;
    LTFAT_TYPE gmax;
    ltfat_int gmaxPos;

    LTFAT_NAME_REAL(findmaxinarray)(g, gl, &gmax, &gmaxPos);
    gthr = gmax * reltol;

    *righttail = 0;
    *lefttail  = 0;

    for (ltfat_int l = 0; l < gl2; l++)
        if (g[l] > gthr)
            *righttail = l + 1;

    for (ltfat_int l = gl - 1; l >= gl2; l--)
        if (g[l] > gthr)
            *lefttail = gl - l;

    return 0;
}


int
LTFAT_NAME(dgtrealmp_kernel_findsmallsize)(const LTFAT_COMPLEX kernlarge[],
        ltfat_int M, ltfat_int N, LTFAT_REAL reltol,
        ksize* size, kanchor* anchor)
{
    LTFAT_REAL kthr;
    LTFAT_COMPLEX maxcoef;
    ltfat_int maxcoefIdx, lastrow = 0, lastcol1 = 0, lastcol2 = 0, M2;
    size->width = 0; size->height = 0;
    anchor->hmid = 0; anchor->wmid = 0;
    M2 = M / 2 + 1;

    LTFAT_NAME_COMPLEX(findmaxinarray)(kernlarge, M * N, &maxcoef, &maxcoefIdx);

    kthr = reltol * ltfat_abs(maxcoef);

    for (ltfat_int n = 0; n < N; n++)
    {
        const LTFAT_COMPLEX* kernlargeCol = kernlarge + n * M;

        for (ltfat_int m = 0; m < M2; m++)
            if ( ltfat_abs(kernlargeCol[m]) > kthr && m > lastrow )
                lastrow = m;
    }

    for (ltfat_int m = 0; m < M2; m++)
    {
        const LTFAT_COMPLEX* kernlargeRow1 = kernlarge + m;
        const LTFAT_COMPLEX* kernlargeRow2 = kernlarge + (N - 1) * M + m;

        for (ltfat_int n = 0; n < N / 2; n++)
        {
            if ( ltfat_abs(kernlargeRow1[n * M]) > kthr && n > lastcol1 )
                lastcol1 = n;

            if ( ltfat_abs(kernlargeRow2[-n * M]) > kthr && n > lastcol2 )
                lastcol2 = n;
        }

    }

    DEBUG("lastrow=%td,lastcol1=%td,lasrcol2=%td", lastrow, lastcol1, lastcol2);

    size->height = 2 * (lastrow + 1) - 1;
    size->width = lastcol1 + lastcol2 + 2;

    anchor->hmid = lastrow ;
    anchor->wmid = lastcol1 ;

    return 0;
}


LTFAT_API int
LTFAT_NAME(dgtrealmp_set_iterstep)(LTFAT_NAME(dgtrealmp_plan)* p,
                                   size_t iterstep)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    CHECK(LTFATERR_NOTPOSARG, iterstep > 0, "iterstep must be greater than 0");
    p->params->iterstep = iterstep;
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(dgtrealmp_set_maxatoms)(LTFAT_NAME(dgtrealmp_plan)* p,
                                   size_t maxatoms)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    CHECK(LTFATERR_NOTPOSARG, maxatoms > 0, "iterstep must be greater than 0");
    p->params->maxatoms = maxatoms;
    p->params->maxit = 2 * maxatoms;
error:
    return status;
}
