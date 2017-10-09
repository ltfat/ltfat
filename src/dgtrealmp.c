#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"
#include "dgtrealmp_private.h"

static inline LTFAT_REAL
ltfat_norm(LTFAT_COMPLEX c)
{
    return ltfat_real(c) * ltfat_real(c) + ltfat_imag(c) * ltfat_imag(c);
}

LTFAT_API int
LTFAT_NAME(dgtrealmp_init_compact)(
    const LTFAT_REAL g[], ltfat_int gl[], ltfat_int L, ltfat_int P, ltfat_int a[],
    ltfat_int M[], ltfat_dgtrealmp_params* params,
    LTFAT_NAME(dgtrealmp_plan)** pout)
{
    int status = LTFATERR_SUCCESS;
    const LTFAT_REAL** multig = NULL;
    CHECK( LTFATERR_NOTPOSARG, P > 0, "P must be positive (passed %td)", P);
    CHECKMEM( multig = LTFAT_NEWARRAY(const LTFAT_REAL*, P));

    for (ltfat_int p = 0, glaccum = 0; p < P;  glaccum += gl[p], p++)
    {
        CHECK(LTFATERR_NOTPOSARG, gl[p] > 0,
              "gl[%td] must be positive (passed %td)", p, gl[p]);
        multig[p] = g + glaccum;
    }

    CHECKSTATUS(
        LTFAT_NAME(dgtrealmp_init)(multig, gl, L, P, a, M, params, pout),
        "dgtrealmp_init failed");

error:
    ltfat_safefree(multig);
    return status;

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
    ltfat_int nextL;
    LTFAT_NAME(dgtrealmp_plan)* p = NULL;
    ltfat_dgtreal_params* dgtparams = NULL;

    CHECK(LTFATERR_NOTPOSARG, P > 0, "P must be positive (passed %td)", P);
    CHECK(LTFATERR_NOTPOSARG, L > 0, "L must be positive (passed %td)", L);
    CHECKNULL(gl); CHECKNULL(g); CHECKNULL(a); CHECKNULL(M); CHECKNULL(pout);

    for (ltfat_int p = 0; p < P; p++)
    {
        CHECKNULL(g[p]);
        CHECK(LTFATERR_NOTPOSARG, gl[p] > 0,
              "gl[%td] must be positive (passed %td)", p, gl[p]);
        CHECK(LTFATERR_NOTPOSARG, a[p] > 0,
              "a[%td] must be positive (passed %td)", p, a[p]);
        CHECK(LTFATERR_NOTPOSARG, M[p] > 0,
              "m[%td] must be positive (passed %td)", p, M[p]);
    }

    nextL = ltfat_dgtlengthmulti(L, P, a, M);

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

#ifdef NOBLASLAPACK
    CHECK( LTFATERR_NOBLASLAPACK,
           p->params->alg == ltfat_dgtrealmp_alg_locomp,
           "LocOMP requires LAPACK, but libltfat was compiled without it.");
#endif

    p->params->initwasrun = 1;

    CHECKMEM( p->dgtplans  = LTFAT_NEWARRAY( LTFAT_NAME(dgtreal_plan)*, P) );
    CHECKMEM( p->gramkerns = LTFAT_NEWARRAY( LTFAT_NAME(kerns)*, P * P) );
    CHECKMEM( p->a  = LTFAT_NEWARRAY( ltfat_int, P));
    CHECKMEM( p->M  = LTFAT_NEWARRAY( ltfat_int, P));
    CHECKMEM( p->M2 = LTFAT_NEWARRAY( ltfat_int, P));
    CHECKMEM( p->N  = LTFAT_NEWARRAY( ltfat_int, P));
    CHECKMEM( p->couttmp = LTFAT_NEWARRAY( LTFAT_COMPLEX*, P));

    for (ltfat_int k = 0; k < P; k++)
    {
        p->a[k] = a[k]; p->M[k] = M[k];
        p->M2[k] = M[k] / 2 + 1; p->N[k] = L / a[k];
    }

    p->P = P; p->L = L;

    CHECKMEM( dgtparams = ltfat_dgtreal_params_allocdef());
    ltfat_dgtreal_setpar_phaseconv(dgtparams, p->params->ptype);

    for (ltfat_int k = 0; k < P; k++)
    {
        CHECKSTATUS(
            LTFAT_NAME(dgtreal_init_gen)(g[k], gl[k], g[k], gl[k], L, 1, a[k], M[k],
                                         NULL, NULL, dgtparams, &p->dgtplans[k]),
            "dgtreal_init failed" );

    }
    ltfat_dgtreal_params_free(dgtparams); dgtparams = NULL;

    for (ltfat_int k1 = 0; k1 < P; k1++)
    {
        for (ltfat_int k2 = 0; k2 < P; k2++)
        {
            gtmp[0] = g[k1]; gtmp[1] = g[k2]; atmp[0] = a[k1]; atmp[1] = a[k2];
            Mtmp[0] = M[k1]; Mtmp[1] = M[k2]; gltmp[0] = gl[k1]; gltmp[1] = gl[k2];

            CHECKSTATUS( LTFAT_NAME(dgtrealmp_kernel_init)( gtmp, gltmp,
                         atmp, Mtmp, L, p->params->kernrelthr,
                         p->params->hint == ltfat_dgtrealmp_allmods,
                         p->params->ptype,
                         &p->gramkerns[k1 + k2 * P]),
                         "dgtrealmp_kernel_init failed");
        }
    }



    /* #ifndef NDEBUG */
    /*     LTFAT_NAME(kerns)* kk = p->gramkerns[0]; */
    /*  */
    /*     for (ltfat_int m = 0; m < kk->size.height; m++ ) */
    /*     { */
    /*         for (ltfat_int n = 0; n < kk->size.width; n++ ) */
    /*         { */
    /*             printf("r=% 5.3e,i=% 5.3e ", ltfat_real(kk->kval[0][n * kk->size.height + m]), */
    /*                    ltfat_imag(kk->kval[0][n * kk->size.height + m])); */
    /*         } */
    /*         printf("\n"); */
    /*     } */
    /* #endif */

    CHECKSTATUS( LTFAT_NAME(dgtrealmpiter_init)(a, M, P, L, &p->iterstate),
                 "dgtrealmpiter_init failed" );

    if (p->params->alg == ltfat_dgtrealmp_alg_locomp)
    {
        // TODO: Compute the right size
        ltfat_int kernSizeAccum = 0;
        for (ltfat_int k = 0; k < P; k++)
            kernSizeAccum += p->gramkerns[k]->size.width * p->gramkerns[k]->size.height;

        CHECKMEM( p->iterstate->gramBuf =
                      LTFAT_NAME_COMPLEX(calloc)( 2 * kernSizeAccum * kernSizeAccum));
        CHECKMEM( p->iterstate->cvalBuf =
                      LTFAT_NAME_COMPLEX(calloc)( 2 * kernSizeAccum));
        CHECKMEM( p->iterstate->cvalinvBuf =
                      LTFAT_NAME_COMPLEX(calloc)( 2 * kernSizeAccum));
        CHECKMEM( p->iterstate->cvalBufPos =
                      LTFAT_NEWARRAY( kpoint, 2 * kernSizeAccum ));

        CHECKSTATUS( LTFAT_NAME_COMPLEX(hermsystemsolver_init)(
                         2 * kernSizeAccum, &p->iterstate->hplan),
                     "hermsystemsolver_init failed" );
    }

    if (p->params->alg == ltfat_dgtrealmp_alg_cyclicmp)
    {
        CHECKMEM( p->iterstate->pointBuf =
                      LTFAT_NEWARRAY( kpoint, p->params->maxatoms) );
    }

    *pout = p;
    return LTFATERR_SUCCESS;
error:
    if (p) LTFAT_NAME(dgtrealmp_done)(&p);
    if (dgtparams) ltfat_dgtreal_params_free(dgtparams);
    *pout = NULL;
    return status;
}


LTFAT_API int
LTFAT_NAME(dgtrealmp_reset)(LTFAT_NAME(dgtrealmp_plan)* p, const LTFAT_REAL* f)
{
    int status = LTFATERR_SUCCESS;

    LTFAT_NAME(dgtrealmpiter_state)* istate = NULL;

    CHECKNULL(p); CHECKNULL(f);
    istate = p->iterstate;

    istate->currit = 0;
    istate->curratoms = 0;
    istate->err = 0.0;

    for (ltfat_int l = 0; l < p->L; l++)
        istate->err += f[l] * f[l];

    istate->fnorm2 = istate->err;
    p->params->errtoladj = powl((long double)10.0,
                                p->params->errtoldb / 10.0) * p->iterstate->fnorm2;

    for (ltfat_int k = 0; k < p->P; k++)
    {
        LTFAT_REAL* sEl = istate->s[k];
        LTFAT_COMPLEX* cEl = istate->c[k];

        CHECKSTATUS(
            LTFAT_NAME(dgtreal_execute_ana_newarray)(p->dgtplans[k], f, cEl),
            "dgtreal_execute failed");

        /* LTFAT_NAME_COMPLEX(dgtreal2dgt)(cEl, p->M[k], p->N[k], cEl); */

        for (size_t l = 0; l < (size_t) ( p->M2[k] * p->N[k]); l++ )
            sEl[l] = ltfat_norm(cEl[l]);

        for (ltfat_int n = 0; n < p->N[k]; n++)
        {
            LTFAT_NAME(maxtree_reset)(istate->fmaxtree[k][n], sEl + n * p->M2[k]);
            LTFAT_NAME(maxtree_findmax)(istate->fmaxtree[k][n],
                                        &istate->maxcols[k][n],
                                        &istate->maxcolspos[k][n]);
        }

        LTFAT_NAME(maxtree_reset)(istate->tmaxtree[k], istate->maxcols[k]);

        memset( p->iterstate->suppind[k], 0 ,
                p->M2[k] * p->N[k] * sizeof * p->iterstate->suppind[k] );
    }

error:
    return status;
}

#define NLOOP \
    for ( ltfat_int nidx = n2start, knidx = kstart2.n; \
          knidx < k->size.width; \
          nidx = ++nidx>=p->N[w2]? nidx - p->N[w2]: nidx, knidx+=k->astep )

#define NLOOPUNBOUNDED \
    for ( ltfat_int nidx = n2start, nidx2 = pos.n2-kmid2.wmid, knidx = kstart2.n; \
          knidx < k->size.width; \
          nidx = ++nidx>=p->N[w2]? nidx-p->N[w2]: nidx, nidx2++, knidx+=k->astep)

#define MLOOP \
    for ( ltfat_int midx = m2start, kmidx = kstart2.m;\
          kmidx < k->size.height;\
          midx = ++midx>=p->M[w2]? midx - p->M[w2]: midx, kmidx += k->Mstep)

#define LTFAT_DGTREALMP_APPLYKERNEL(ctmp, SIGN)\
NLOOP{\
    LTFAT_COMPLEX* currcCol = s->c[w2] + nidx * p->M2[w2];\
    LTFAT_REAL*    currsCol = s->s[w2] + nidx * p->M2[w2];\
    LTFAT_COMPLEX* kcurrCol = kvals + knidx * k->size.height;\
MLOOP{\
    if(!(midx >=0 && midx < p->M2[w2])) continue;\
    currcCol[midx] = currcCol[midx] SIGN ctmp * kcurrCol[kmidx];\
    currsCol[midx] = ltfat_norm( currcCol[midx]);\
}}

#define LTFAT_DGTREALMP_MARKMODIFIED \
NLOOP{\
    LTFAT_NAME(maxtree_setdirty)(s->fmaxtree[w2][nidx], m2start, m2start + kdim2.height);}\
    LTFAT_NAME(maxtree_setdirty)(s->tmaxtree[w2],       n2start, n2start + kdim2.width);

LTFAT_API int
LTFAT_NAME(dgtrealmp_execute_niters)(
    LTFAT_NAME(dgtrealmp_plan)* p, ltfat_int itno, LTFAT_COMPLEX** cout)
{
    int status = LTFAT_DGTREALMP_STATUS_CANCONTINUE;

    LTFAT_NAME(dgtrealmpiter_state)* s = p->iterstate;

    for (ltfat_int iter = 0;
         iter < itno && status == LTFAT_DGTREALMP_STATUS_CANCONTINUE;
         iter++)
    {
        s->currit++;

        kpoint origpos;
        LTFAT_NAME(dgtrealmp_execute_findmaxatom)(p, &origpos);

        // Only increase the number of used atoms if the atom has not been selected
        // previously
        if ( !s->suppind[PTOI(origpos)] ) s->curratoms++;

        int uniquenyquest = p->M[origpos.w] % 2 == 0;
        /* int do_conj = !( origpos.m == 0 || (origpos.m == p->M2[w] - 1 && uniquenyquest)); */

        switch ( p->params->alg)
        {
        case ltfat_dgtrealmp_alg_mp:
        {
            LTFAT_REAL cvalenergy =
                LTFAT_NAME(dgtrealmp_execute_mp)( p, origpos, cout);
            s->suppind[PTOI(origpos)] += cvalenergy;
            s->err -= cvalenergy;
            break;
        }
        case ltfat_dgtrealmp_alg_locomp:
        {
            int count = s->suppind[PTOI(origpos)];
            DEBUG("\n*****************\n Count %d \n ****************", count);
            if (count > 10)
            {
                /* status = LTFAT_DGTREALMP_STATUS_LOCOMP_ORTHFAILED; break; */
            }

            // STEP 1: Find all active atoms around the current one
            ltfat_int cvalNo = 0;

            LTFAT_NAME(kerns)* k1 = p->gramkerns[origpos.w + s->P * origpos.w];
            int do_conjat = 0;
            if (  (k1->mid.hmid + 2 * origpos.m) < k1->size.height )
                do_conjat = 1;

            s->cvalBuf[cvalNo] = s->c[PTOI(origpos)];
            s->cvalBufPos[cvalNo] = origpos;
            cvalNo++;


            for (ltfat_int w2 = 0; w2 < s->P; w2++)
            {
                kpoint pos = kpoint_init(0, 0, w2);
                ltfat_int m2start, n2start;
                ksize   kdim2; kanchor kmid2; kpoint  kstart2;

                LTFAT_NAME(dgtrealmp_execute_indices)(
                    p, origpos, &pos, &m2start, &n2start,
                    &kdim2, &kmid2, &kstart2);

                LTFAT_NAME(kerns)* k = p->gramkerns[origpos.w + s->P * w2];

                NLOOPUNBOUNDED
                {
                    LTFAT_REAL* suppCol = s->suppind[w2] + nidx * p->M2[w2];

                    MLOOP
                    {
                        if ( midx == origpos.m && nidx == origpos.n && w2 == origpos.w) continue;
                        if ( midx >= 0 && midx < p->M2[w2] && suppCol[midx])
                        {
                            kpoint cvalPos = kpoint_init2(midx, nidx, nidx2, w2);
                            s->cvalBuf[cvalNo] = s->c[PTOI(cvalPos)];
                            s->cvalBufPos[cvalNo] = cvalPos;
                            cvalNo++;

                            if ( midx > 0 && do_conjat )
                            {
                                /*     cvalPos.m =  -midx; */
                                /*     s->cvalBuf[cvalNo] = conj(s->cvalBuf[cvalNo - 1]); */
                                /*     s->cvalBufPos[cvalNo] = cvalPos; */
                                /*     cvalNo++; */
                                DEBUGNOTE("PRD***********************************");
                            }
                        }
                    }
                }
            }

#ifndef NDEBUG
            for (ltfat_int cidx = 0; cidx < cvalNo; cidx++)
            {
                kpoint cvalPos = s->cvalBufPos[cidx];
                LTFAT_COMPLEX cval =  s->cvalBuf[cidx];
                DEBUG("m=%td,n=%td,w=%td, r=% 5.3e,i=% 5.3e",
                      cvalPos.m, cvalPos.n, cvalPos.w, ltfat_real(cval), ltfat_imag(cval));
            }
            DEBUGNOTE("--------------------");
#endif

            /* cvalNo = 1; */
            /* kpoint cvalPos; cvalPos.m = m; cvalPos.n = n; cvalPos.w = w; */
            /* s->cvalBuf[0] = s->c[PTOI(cvalPos)]; */
            /* s->cvalBufPos[0] = cvalPos; */

            //memset(s->gramBuf, 0, cvalNo * cvalNo * sizeof * s->gramBuf);
            // STEP 2: Construct the Gram matrix
            for (ltfat_int cidx1 = 0; cidx1 < cvalNo; cidx1++)
            {
                kpoint cvalPos            = s->cvalBufPos[cidx1];
                LTFAT_COMPLEX* gramBufCol = s->gramBuf + cidx1 * cvalNo;

                gramBufCol[cidx1] = 1;

                for (ltfat_int cidx2 = cidx1 + 1; cidx2 < cvalNo; cidx2++)
                {
                    kpoint cvalPos2      = s->cvalBufPos[cidx2];
                    kpoint pos           = cvalPos2;
                    LTFAT_NAME(kerns)* k =
                        p->gramkerns[cvalPos.w + s->P * cvalPos2.w];

                    LTFAT_COMPLEX* kvals =
                        LTFAT_NAME(dgtrealmp_execute_pickkernel)(
                            k, cvalPos.m, cvalPos.n, p->params->ptype);

                    ltfat_int m2start, n2start;
                    ksize   kdim2; kanchor kmid2; kpoint kstart2;

                    LTFAT_NAME(dgtrealmp_execute_indices)(
                        p, cvalPos, &pos, &m2start, &n2start, &kdim2,
                        &kmid2, &kstart2);

                    ltfat_int muse = cvalPos2.m  - pos.m  + kmid2.hmid;
                    ltfat_int nuse = cvalPos2.n2 - pos.n2 + kmid2.wmid;

                    if ( muse >= 0 && muse < kdim2.height &&
                         nuse >= 0 && nuse < kdim2.width )
                        gramBufCol[cidx2] =
                            (kvals[k->size.height * (kstart2.n + k->astep * nuse) +
                                   k->Mstep * muse + kstart2.m]);
                    else
                        gramBufCol[cidx2] = 0;
                }
            }

            /* #ifndef NDEBUG */
            /*             printf("\n"); */
            /*             for (ltfat_int m = 0; m < cvalNo; m++ ) */
            /*             { */
            /*                 for (ltfat_int n = 0; n < cvalNo; n++ ) */
            /*                 { */
            /*                     printf("r=% 2.3e,i=% 2.3e ", ltfat_real(s->gramBuf[n * cvalNo + m]), */
            /*                            ltfat_imag(s->gramBuf[n * cvalNo + m])); */
            /*                 } */
            /*                 printf("\n"); */
            /*             } */
            /* #endif */

            memcpy(s->cvalinvBuf, s->cvalBuf, cvalNo * sizeof * s->cvalinvBuf);
            // STEP 3: Invert that S**T
            if ( LTFAT_NAME_COMPLEX(hermsystemsolver_execute)(
                     s->hplan, s->gramBuf, cvalNo, s->cvalinvBuf) )
            {
                /* status = LTFATERR_NOTPOSDEFMATRIX; */
                status = LTFAT_DGTREALMP_STATUS_LOCOMP_NOTHERM;
                break;
            }

#ifndef NDEBUG
            for (ltfat_int cidx = 0; cidx < cvalNo; cidx++)
            {
                kpoint cvalPos = s->cvalBufPos[cidx];
                LTFAT_COMPLEX cval =  s->cvalinvBuf[cidx];
                DEBUG("m=%td,n=%td,n2=%td, r=% 2.3e,i=% 2.3e",
                      cvalPos.m, cvalPos.n, cvalPos.n2, ltfat_real(cval), ltfat_imag(cval));
            }
            DEBUGNOTE("------------+-------");
#endif

            // STEP 4: Update result and the residuum
            for (ltfat_int cidx = 0; cidx < cvalNo; cidx++)
            {
                LTFAT_COMPLEX cvalinv = s->cvalinvBuf[cidx];
                LTFAT_COMPLEX cval = s->cvalBuf[cidx];
                kpoint cvalPos     = s->cvalBufPos[cidx];
                if (cvalPos.m < 0) continue;

                int do_conj = !( cvalPos.m == 0 || (cvalPos.m == p->M2[cvalPos.w] - 1
                                                    && uniquenyquest));

                LTFAT_REAL fac =
                    LTFAT_NAME(dgtrealmp_execute_normfac)( p, cvalPos, cvalinv);
                cvalinv *= fac;

                cout[PTOI(cvalPos)] += cvalinv;

                /* LTFAT_REAL cvalabs = ltfat_norm(cvalinv); */
                /*  */
                /* s->err -= cvalabs / fac; */
                /* if (do_conj) s->err -= cvalabs / fac; */

                LTFAT_NAME(dgtrealmp_execute_updateresiduum)( p,
                        cvalPos, cvalinv,  1);

            }




            break;
        }

        case ltfat_dgtrealmp_alg_cyclicmp:
        {
            LTFAT_REAL cvalenergy =
                LTFAT_NAME(dgtrealmp_execute_mp)( p, origpos, cout);
            s->suppind[PTOI(origpos)] += cvalenergy;
            s->err -= cvalenergy;

            /* DEBUG("Start: %Lf", s->err); */

            s->pointBufNo = 0;
            for (ltfat_int w2 = 0; w2 < s->P; w2++)
            {
                ltfat_int m2start, n2start;
                ksize   kdim2; kanchor kmid2; kpoint  kstart2;
                kpoint pos; pos.w = w2;

                LTFAT_NAME(dgtrealmp_execute_indices)(
                    p, origpos, &pos, &m2start, &n2start,
                    &kdim2, &kmid2, &kstart2);

                LTFAT_NAME(kerns)* k = p->gramkerns[origpos.w + s->P * w2];

                NLOOP
                {
                    LTFAT_REAL* suppCol = s->suppind[w2] + nidx * p->M2[w2];
                    MLOOP
                    if ( !(midx == origpos.m && nidx == origpos.n && w2 == origpos.w))
                        if ( midx >= 0 && midx < p->M2[w2] && suppCol[midx])
                            s->pointBuf[s->pointBufNo++] = kpoint_init(midx, nidx, w2);
                }
            }

            /* break; */
            if (s->pointBufNo == 0 ) break;
            /* s->pointBuf[s->pointBufNo++] = origpos; */
            /* s->pointBuf[0] = origpos; */
            /* s->pointBufNo = 1; */


            for ( ltfat_int cycl = 0; cycl < 1; cycl++ )
            {
                for (size_t pIdx = 0; pIdx < s->pointBufNo ; pIdx++)
                {
                    kpoint* pos = &s->pointBuf[pIdx];
                    LTFAT_COMPLEX cval = cout[PTOI((*pos))];
                    /* s->err += s->suppind[PTOI((*pos))]; */

                    s->curratoms -= 1;
                    cout[PTOI((*pos))]       = 0;


                    LTFAT_REAL cvale = ltfat_norm(cval);
                    int do_conj = !( pos->m == 0 || (pos->m == p->M2[pos->w] - 1
                                                     && uniquenyquest));

                    LTFAT_REAL fac =
                        LTFAT_NAME(dgtrealmp_execute_normfac)( p, *pos, cval);

                    cvale /= fac;
                    if (do_conj) cvale *= 2.0;
                    s->suppind[PTOI((*pos))] -=  cvale;
                    s->err += cvale;

                    /* DEBUG("Middle %Lf", s->err); */

                    LTFAT_NAME(dgtrealmp_execute_updateresiduum)(
                        p, *pos, cval, 0);

                    LTFAT_NAME(dgtrealmp_execute_findmaxatom)(p, pos);

                    if ( !s->suppind[PTOI((*pos))] ) s->curratoms++;

                    cvalenergy = LTFAT_NAME(dgtrealmp_execute_mp)( p, *pos, cout);
                    s->suppind[PTOI((*pos))] += cvalenergy;
                    s->err -= cvalenergy;
                    /* DEBUG("End %Lf", s->err); */
                }
            }



            /* s->err -= cvalenergy; */

            /* for (ltfat_int w2 = 0; w2 < s->P; w2++) */
            /* { */
            /*     ltfat_int m2start, n2start, Mstep, astep, m2, n2; */
            /*     ksize   kdim2; kanchor kmid2; kpoint  kstart2; */
            /*  */
            /*     LTFAT_NAME(dgtrealmp_execute_indices)( */
            /*         p, m, n, w, w2, &m2, &n2, &m2start, &n2start, */
            /*         &Mstep, &astep, &kdim2, &kmid2, &kstart2); */
            /*  */
            /*     LTFAT_NAME(kerns)* k = p->gramkerns[w + s->P * w2]; */
            /*  */
            /*     NLOOP */
            /*     { */
            /*         int* suppCol = s->suppind[w2] + nidx * p->M2[w2]; */
            /*  */
            /*         MLOOP */
            /*         { */
            /*             if ( midx >= 0 && midx < p->M2[w2] && suppCol[midx]) */
            /*             { */
            /*                 for ( ltfat_int cycl = 0; cycl < 2; cycl++ ) */
            /*                 { */
            /*                     kpoint cvalPos; */
            /*                     cvalPos.m = midx; cvalPos.n = nidx; cvalPos.w = w2; */
            /*                     kpoint* pos = &cvalPos; */
            /*                     s->suppind[PTOI((*pos))] = 0; */
            /*                     LTFAT_COMPLEX* cval = &cout[PTOI((*pos))]; */
            /*                     LTFAT_REAL cvalabs = ltfat_norm(*cval); */
            /*                     do_conj = !( pos->m == 0 || (pos->m == p->M2[pos->w] - 1 */
            /*                     && uniquenyquest)); */
            /*                     #<{(| s->err += cvalabs; |)}># */
            /*                     #<{(| if (do_conj) s->err += cvalabs; |)}># */
            /*  */
            /*                     LTFAT_NAME(dgtrealmp_execute_updateresiduum)( */
            /*                         p, *cval, pos->m, pos->n, pos->w, do_conj, 0); */
            /*  */
            /*                     LTFAT_NAME(dgtrealmp_execute_findmaxatom)(p, &pos->m, &pos->n, &pos->w); */
            /*  */
            /*                     *cval = 0; */
            /*                     LTFAT_NAME(dgtrealmp_execute_mp)( p, pos->m, pos->n, pos->w, cout); */
            /*                     s->suppind[PTOI((*pos))] = 1; */
            /*                 } */
            /*  */
            /*             } */
            /*         } */
            /*     } */
            /* } */


            break;
        }
        }

        /* s->err = s->fnorm2; */
        /* for (ltfat_int w = 0; w < p->P; w++ ) */
        /* { */
        /*     for (ltfat_int n = 0; n < p->N[w]; n++) */
        /*     { */
        /*         for (ltfat_int m = 0; m < p->M2[w]; m++) */
        /*         { */
        /*             int do_conj = !( m == 0 || (m == p->M2[w] - 1 */
        /*                                         && uniquenyquest)); */
        /*  */
        /*                 LTFAT_COMPLEX cval = cout[w][m + n * p->M2[w]]; */
        /*                 LTFAT_REAL cvalabs = ltfat_norm(cval); */
        /*                 LTFAT_REAL fac = */
        /*                     LTFAT_NAME(dgtrealmp_execute_normfac)( p, kpoint_init(m, n, w), cval); */
        /*                 s->err -= cvalabs / fac; */
        /*                 if (do_conj) */
        /*                     s->err -= cvalabs / fac; */
        /*         } */
        /*     } */
        /* } */

        DEBUG("ERR=%2.3Lf", s->err);

        if (s->err < 0)
        {
            status = LTFAT_DGTREALMP_STATUS_STALLED;
            break;
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


LTFAT_API int
LTFAT_NAME(dgtrealmp_revert)(
    LTFAT_NAME(dgtrealmp_plan)* p, LTFAT_COMPLEX** cout)
{

    for (ltfat_int w = 0; w < p->P; w++)
        for (ltfat_int n = 0; n < p->N[w]; n++)
            for (ltfat_int m = 0; m < p->M2[w]; m++)
            {
                kpoint pos = kpoint_init(m, n, w);
                if (cout[PTOI(pos)])
                {
                    LTFAT_REAL cvalenergy =
                        LTFAT_NAME(dgtrealmp_execute_invmp)(
                            p, pos, cout);

                    p->iterstate->suppind[PTOI(pos)] -= cvalenergy;
                    p->iterstate->err += cvalenergy;
                }
            }
    return 0;
}

LTFAT_REAL
LTFAT_NAME(dgtrealmp_execute_invmp)( LTFAT_NAME(dgtrealmp_plan)* p,
                                     kpoint pos, LTFAT_COMPLEX** cout)
{
    /* LTFAT_NAME(dgtrealmpiter_state)* s = p->iterstate; */

    int uniquenyquest = p->M[pos.w] % 2 == 0;
    int do_conj = !( pos.m == 0 || (pos.m == p->M2[pos.w] - 1
                                    && uniquenyquest));

    LTFAT_COMPLEX cval = cout[PTOI(pos)];
    cout[PTOI(pos)] = 0.0;

    LTFAT_REAL fac =
        LTFAT_NAME(dgtrealmp_execute_normfac)( p, pos, cval);

    cval *= fac;

    LTFAT_REAL cvalenergy = ltfat_norm(cval) / fac;
    if (do_conj) cvalenergy *= 2.0;
    /* s->err -= cvalenergy / fac; */
    /* if (do_conj) s->err -= cvalenergy / fac; */

    LTFAT_NAME(dgtrealmp_execute_updateresiduum)(
        p, pos, cval, 0);

    return cvalenergy;
}

LTFAT_REAL
LTFAT_NAME(dgtrealmp_execute_mp)( LTFAT_NAME(dgtrealmp_plan)* p,
                                  kpoint pos, LTFAT_COMPLEX** cout)
{
    LTFAT_NAME(dgtrealmpiter_state)* s = p->iterstate;

    int uniquenyquest = p->M[pos.w] % 2 == 0;
    int do_conj = !( pos.m == 0 || (pos.m == p->M2[pos.w] - 1
                                    && uniquenyquest));

    // Max complex atom
    LTFAT_COMPLEX cval = s->c[PTOI(pos)];

    LTFAT_REAL fac =
        LTFAT_NAME(dgtrealmp_execute_normfac)( p, pos, cval);

    cval *= fac;

    cout[PTOI(pos)] += cval;

    LTFAT_REAL cvalenergy = ltfat_norm(cval) / fac;
    if (do_conj) cvalenergy *= 2.0;
    /* s->err -= cvalenergy / fac; */
    /* if (do_conj) s->err -= cvalenergy / fac; */

    LTFAT_NAME(dgtrealmp_execute_updateresiduum)(
        p, pos, cval, 1);

    return cvalenergy;
}


LTFAT_REAL
LTFAT_NAME(dgtrealmp_execute_normfac)(
    LTFAT_NAME(dgtrealmp_plan)* p,
    kpoint pos, LTFAT_COMPLEX cval)
{
    LTFAT_NAME(dgtrealmpiter_state)* s = p->iterstate;
    int uniquenyquest = p->M[pos.w] % 2 == 0;
    int do_conj = !( pos.m == 0 || (pos.m == p->M2[pos.w] - 1 && uniquenyquest));

    LTFAT_NAME(kerns)* k = p->gramkerns[pos.w + s->P * pos.w];
    LTFAT_COMPLEX* kvals  =
        LTFAT_NAME(dgtrealmp_execute_pickkernel)( k, pos.m, pos.n, p->params->ptype);

    /* Check whether the kernel overflows */
    /* If it does, adjust cval by the inner product */
    ltfat_int posinkern  = k->mid.hmid + 2 * pos.m;
    /* ltfat_int posinkern2 = k->mid.hmid - 2 * m; */
    if (do_conj &&
        (posinkern < k->size.height ) )
    {
        LTFAT_COMPLEX cvalphase = exp( I * ((LTFAT_REAL)2.0) * ltfat_arg(cval));
        LTFAT_COMPLEX atinprod =
            kvals[k->size.height * k->mid.wmid + posinkern];

        return ( 1.0 / (1.0 + ltfat_real(cvalphase * conj( atinprod) )));
        //*cval *= *fac;
    }
    else
    {
        return 1.0;
    }
    return 0;
}

int
LTFAT_NAME(dgtrealmp_execute_updateresiduum)(
    LTFAT_NAME(dgtrealmp_plan)* p,
    kpoint origpos, LTFAT_COMPLEX cval,
    int do_substract)
{
    ltfat_int  m2start, n2start;
    kpoint pos;
    ksize   kdim2; kanchor kmid2; kpoint  kstart2;

    int uniquenyquest = p->M[origpos.w] % 2 == 0;
    int do_conj = !( origpos.m == 0 ||
                     (origpos.m == p->M2[origpos.w] - 1 && uniquenyquest));

    LTFAT_NAME(dgtrealmpiter_state)* s = p->iterstate;
    LTFAT_COMPLEX cval2 = conj(cval);
    /* ltfat_int Morig  = p->M[pos.w]; */
    /* ltfat_int M2orig  = p->M2[w]; */
    kpoint origposconj = origpos;
    origposconj.m = p->M[origpos.w] - origpos.m;


    /* This loop is trivially pararelizable */
    for (ltfat_int w2 = 0; w2 < s->P; w2++)
    {
        LTFAT_NAME(kerns)* k     = p->gramkerns[origpos.w + s->P * w2];
        LTFAT_COMPLEX*     kvals = LTFAT_NAME(dgtrealmp_execute_pickkernel)(
                                       k, origpos.m, origpos.n, p->params->ptype);

        pos.w = w2;
        LTFAT_NAME(dgtrealmp_execute_indices)(
            p, origpos, &pos, &m2start, &n2start, &kdim2, &kmid2, &kstart2);

        if (do_substract)
        {
            LTFAT_DGTREALMP_APPLYKERNEL(cval, -)
        }
        else
        {
            LTFAT_DGTREALMP_APPLYKERNEL(cval, +)
        }
        LTFAT_DGTREALMP_MARKMODIFIED

        ltfat_int posinkern  = kmid2.hmid + 2 * pos.m;
        /* posinkern2 = kmid2.hmid - 2 * m2; */
        if (do_conj && (posinkern < kdim2.height ))
        {
            kvals = LTFAT_NAME(dgtrealmp_execute_pickkernel)(
                        k, origposconj.m, origpos.n, p->params->ptype);

            LTFAT_NAME(dgtrealmp_execute_indices)(
                p, origposconj, &pos, &m2start, &n2start, &kdim2, &kmid2, &kstart2);

            if (do_substract)
            {
                LTFAT_DGTREALMP_APPLYKERNEL(cval2, -)
            }
            else
            {
                LTFAT_DGTREALMP_APPLYKERNEL(cval2, +)
            }
        }
    }
    return 0;
}

inline LTFAT_COMPLEX*
LTFAT_NAME(dgtrealmp_execute_pickkernel)(
    LTFAT_NAME(kerns)* k, ltfat_int m, ltfat_int n,
    ltfat_phaseconvention pconv)
{
    if (pconv == LTFAT_FREQINV)
        return k->kval[ltfat_positiverem( n, k->kNo)];
    else if (pconv == LTFAT_TIMEINV)
        return k->kval[ltfat_positiverem( m, k->kNo)];
    else
        return NULL;
}

/* inline int */
/* LTFAT_NAME(dgtrealmp_execute_kpos)(LTFAT_NAME(dgtrealmp_plan)* p, */
/*                                    kpoint pos1, kpoint pos2, */
/*                                    ltfat_int* m2, ltfat_int* n2, */
/*                                    ltfat_int* Mstep, ltfat_int* astep, */
/*                                    ksize* kdim2, kanchor* kmid2, kpoint* kstart2) */
/* { */
/*  */
/*  */
/*     return 0; */
/* } */

inline int
LTFAT_NAME(dgtrealmp_execute_indices)(LTFAT_NAME(dgtrealmp_plan)* p,
                                      kpoint origpos, kpoint* pos,
                                      /* ltfat_int m, ltfat_int n, ltfat_int w,   */
                                      /* ltfat_int w2, ltfat_int* m2, ltfat_int* n2,  */
                                      ltfat_int* m2start, ltfat_int* n2start,
                                      ksize* kdim2, kanchor* kmid2, kpoint* kstart2)
{
    LTFAT_NAME(kerns)* k     = p->gramkerns[origpos.w + p->P * pos->w];

    pos->n    = (ltfat_int) ltfat_round( origpos.n / k->arat);
    pos->m    = (ltfat_int) ltfat_round( origpos.m / k->Mrat);

    ltfat_int n2off = origpos.n - (ltfat_int)(pos->n * k->arat);
    ltfat_int m2off = origpos.m - (ltfat_int)(pos->m * k->Mrat);

    *kdim2 = k->size;
    *kmid2 = k->mid;
    kmid2->hmid = (k->mid.hmid - m2off) / k->Mstep;
    kmid2->wmid = (k->mid.wmid - n2off) / k->astep;

    *m2start = pos->m - kmid2->hmid;
    *m2start = *m2start < 0 ? *m2start + p->M[pos->w] : *m2start;
    *n2start = pos->n - kmid2->wmid;
    *n2start = *n2start < 0 ? *n2start + p->N[pos->w] : *n2start;

    kstart2->m = k->mid.hmid - m2off - kmid2->hmid * k->Mstep;
    kstart2->n = k->mid.wmid - n2off - kmid2->wmid * k->astep;

    kdim2->height = (kdim2->height - kstart2->m) / k->Mstep;
    kdim2->width  = (kdim2->width  - kstart2->n) / k->astep;

    return 0;
}

inline int
LTFAT_NAME(dgtrealmp_execute_findmaxatom)(
    LTFAT_NAME(dgtrealmp_plan)* p, kpoint* pos)
/* ltfat_int* m, ltfat_int* n, ltfat_int* w) */
{
    LTFAT_NAME(dgtrealmpiter_state)* s = p->iterstate;
    LTFAT_REAL val = 0.0;
    for (ltfat_int k = 0; k < s->P; k++)
    {
        LTFAT_REAL valTmp; ltfat_int nTmp;
        ltfat_int dirtystart, dirtyend;
        LTFAT_NAME(maxtree_getdirty)(s->tmaxtree[k], &dirtystart, &dirtyend);

        ltfat_int N = p->N[k];
        /* DEBUG("start=%td, end=%td, N=%td",dirtystart, dirtyend, N); */
        dirtystart = ltfat_positiverem(dirtystart, N);
        dirtyend =   ltfat_positiverem(dirtyend,   N);

        /* DEBUG("start=%td, end=%td, N=%td",dirtystart, dirtyend, N); */

        for (ltfat_int nidx = dirtystart; nidx != dirtyend;
             nidx = ++nidx >= N ? nidx - N : nidx)
            LTFAT_NAME(maxtree_findmax)( s->fmaxtree[k][nidx],
                                         &s->maxcols[k][nidx],
                                         &s->maxcolspos[k][nidx]);

        LTFAT_NAME(maxtree_findmax)(s->tmaxtree[k], &valTmp, &nTmp);

        if ( valTmp > val )
        {
            val = valTmp; pos->m = s->maxcolspos[k][nTmp]; pos->n = nTmp; pos->w = k;
        }
    }


    return 0;
}

#undef LTFAT_DGTREALMP_APPLYKERNEL
#undef LTFAT_DGTREALMP_MARKMODIFIED
#undef NLOOP
#undef NLOOPUNBOUNDED
#undef MLOOP


LTFAT_API int
LTFAT_NAME(dgtrealmp_execute_compact)(LTFAT_NAME(dgtrealmp_plan)* p,
                                      const LTFAT_REAL* f,
                                      LTFAT_COMPLEX* cout,
                                      LTFAT_REAL* fout)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    for (ltfat_int k = 0, accum = 0; k < p->P; accum += p->N[k] * p->M2[k], k++)
        p->couttmp[k] = cout + accum;

    status = LTFAT_NAME(dgtrealmp_execute)( p, f, p->couttmp, fout);
error :
    return status;

}

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

    while ( LTFAT_DGTREALMP_STATUS_CANCONTINUE ==
            ( status2 = LTFAT_NAME(dgtrealmp_execute_niters)( p, iterstep, cout)))
    {
        if (p->params->verbose)
        {
            printf("Error %.2Lf, iteration %zu \n",
                   (long double) 10.0 * log10l( p->iterstate->err / p->iterstate->fnorm2 ),
                   p->iterstate->currit);
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
    ltfat_safefree(pp->couttmp);

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
    CHECKMEM( s->N = LTFAT_NEWARRAY(ltfat_int, P));
    CHECKMEM( s->suppind = LTFAT_NEWARRAY(LTFAT_REAL*, P));
    CHECKMEM( s->maxcols    =  LTFAT_NEWARRAY(LTFAT_REAL*, P));
    CHECKMEM( s->maxcolspos =  LTFAT_NEWARRAY(ltfat_int*, P));
    CHECKMEM( s->tmaxtree =  LTFAT_NEWARRAY( LTFAT_NAME(maxtree)*, P));
    CHECKMEM( s->fmaxtree =  LTFAT_NEWARRAY( LTFAT_NAME(maxtree)**, P));
    /* s->fmaxtree = NULL; */
    s->P = P;

    for (ltfat_int p = 0; p < P; p++)
    {
        ltfat_int N = L / a[p];
        s->N[p] = N;
        ltfat_int M2 = M[p] / 2 + 1;
        /* ltfat_int M2ext = M2 */
        CHECKMEM( s->s[p] = LTFAT_NAME_REAL(malloc)(N * M2) );
        CHECKMEM( s->c[p] = LTFAT_NAME_COMPLEX(malloc)(N * M2) );
        CHECKMEM( s->suppind[p] = LTFAT_NEWARRAY(LTFAT_REAL, N * M2 ));
        CHECKMEM( s->maxcols[p]    = LTFAT_NAME_REAL(malloc)(N) );
        CHECKMEM( s->maxcolspos[p] = LTFAT_NEWARRAY(ltfat_int, N) );
        CHECKSTATUS( LTFAT_NAME(maxtree_init)(N, N, 10, &s->tmaxtree[p]),
                     "maxtree_init failed" );

        CHECKMEM( s->fmaxtree[p] = LTFAT_NEWARRAY(LTFAT_NAME(maxtree)*, N));
        for (ltfat_int n = 0; n < N; n++ )
            CHECKSTATUS( LTFAT_NAME(maxtree_init)(M2, M[p], 8, &s->fmaxtree[p][n]),
                         "maxtree_init failed" );

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

    if (s->suppind)
    {
        for (ltfat_int p = 0; p < s->P; p++)
            ltfat_safefree(s->suppind[p]);

        ltfat_free(s->suppind);
    }

    if (s->maxcols)
    {
        for (ltfat_int p = 0; p < s->P; p++)
            ltfat_safefree(s->maxcols[p]);

        ltfat_free(s->maxcols);
    }

    if (s->tmaxtree)
    {
        for (ltfat_int p = 0; p < s->P; p++)
            LTFAT_NAME(maxtree_done)(&s->tmaxtree[p]);

        ltfat_free(s->tmaxtree);
    }

    if (s->fmaxtree)
    {
        for (ltfat_int p = 0; p < s->P; p++)
        {
            for (ltfat_int n = 0; n < s->N[p]; n++)
                LTFAT_NAME(maxtree_done)(&s->fmaxtree[p][n]);

            ltfat_free(s->fmaxtree[p]);
        }

        ltfat_free(s->fmaxtree);
    }

    ltfat_safefree(s->gramBuf);
    ltfat_safefree(s->cvalBuf);
    ltfat_safefree(s->cvalinvBuf);
    ltfat_safefree(s->cvalBufPos);
    ltfat_safefree(s->pointBuf);
    if (s->hplan) LTFAT_NAME_COMPLEX(hermsystemsolver_done)(&s->hplan);
    ltfat_safefree(s->N);
    ltfat_free(s);
    *state = NULL;
error:
    return status;
}

int
LTFAT_NAME(dgtrealmp_kernel_init)(
    const LTFAT_REAL* g[], ltfat_int gl[], ltfat_int a[], ltfat_int M[],
    ltfat_int L, LTFAT_REAL reltol, int do_allmods, ltfat_phaseconvention ptype,
    LTFAT_NAME(kerns)** pout)
{
    ltfat_int modNo, amin, Mmax, lefttail0, righttail0, lefttail1, righttail1,
              Lshort, Nshort;
    LTFAT_REAL* g0tmp = NULL, *g1tmp = NULL;
    LTFAT_COMPLEX* kernlarge = NULL;
    LTFAT_NAME(kerns)* ktmp = NULL;
    int status = LTFATERR_SUCCESS;

    CHECKMEM( ktmp = LTFAT_NEW(LTFAT_NAME(kerns)) );
    ktmp->arat = 1.0; ktmp->Mrat = 1.0;


    ktmp->Mrat = ((double) M[0]) / M[1];
    ktmp->arat = ((double) a[1]) / a[0];

    amin = a[0] < a[1] ? a[0] : a[1];
    Mmax = M[0] > M[1] ? M[0] : M[1];

    LTFAT_NAME(dgtrealmp_essentialsupport)(g[0], gl[0], 1e-6, &lefttail0,
                                           &righttail0);
    LTFAT_NAME(dgtrealmp_essentialsupport)(g[1], gl[1], 1e-6, &lefttail1,
                                           &righttail1);

    /* DEBUG("lefttail=%td, righttail=%td", lefttail0, righttail0); */

    Lshort =
        (lefttail0 > righttail0 ? 2 * lefttail0 : 2 * righttail0) +
        (lefttail1 > righttail1 ? 2 * lefttail1 : 2 * righttail1);

    Lshort = ltfat_dgtlength(Lshort > L ? L : Lshort, amin, Mmax);

    /* DEBUG("Lshort=%td", Lshort); */

    Nshort = Lshort / amin;

    CHECKMEM(g0tmp     = LTFAT_NAME_REAL(malloc)(Lshort));
    CHECKMEM(g1tmp     = LTFAT_NAME_REAL(malloc)(Lshort));
    CHECKMEM(kernlarge = LTFAT_NAME_COMPLEX(malloc)(Mmax * Nshort));

    LTFAT_NAME(middlepad)(g[0], gl[0], LTFAT_WHOLEPOINT, Lshort, g0tmp);
    LTFAT_NAME(middlepad)(g[1], gl[1], LTFAT_WHOLEPOINT, Lshort, g1tmp);

    LTFAT_NAME_REAL(dgt_long)(g0tmp, g1tmp, Lshort, 1, amin, Mmax,
                              ptype, kernlarge);

    LTFAT_NAME(dgtrealmp_kernel_findsmallsize)(kernlarge, Mmax, Nshort,
            reltol, &ktmp->size, &ktmp->mid);

    LTFAT_NAME_COMPLEX(circshift2)(kernlarge, Mmax, Nshort,
                                   ktmp->mid.hmid, ktmp->mid.wmid, kernlarge);

    modNo = ltfat_lcm(amin, Mmax) / amin;

    ktmp->Mstep = ktmp->Mrat > 1 ? (ltfat_int) ktmp->Mrat : 1;
    ktmp->astep = ktmp->arat > 1 ? (ltfat_int) ktmp->arat : 1;

    ltfat_int kernskip = 1;
    if (ptype == LTFAT_FREQINV)
    {
        if (a[0] > a[1])
        {
            modNo = ltfat_lcm(amin, Mmax) / a[0];
            kernskip = a[0] / a[1];
        }
    }
    else if (ptype == LTFAT_TIMEINV)
    {
        if (M[1] > M[0])
        {
            modNo = ltfat_lcm(amin, Mmax) / ktmp->Mrat;
            kernskip = M[1] / M[0];
        }
    }


    ktmp->kNo = modNo;

    /* DEBUG("h=%td,w=%td,hmid=%td,wmid=%td,kno=%td", */
    /*       ktmp->size.height, ktmp->size.width, ktmp->mid.hmid, ktmp->mid.wmid, modNo ); */

    CHECKMEM( ktmp->kval = LTFAT_NEWARRAY(LTFAT_COMPLEX*, ktmp->kNo));

    if (do_allmods)
    {
        for (ltfat_int k = 0; k < ktmp->kNo; k++)
            CHECKMEM( ktmp->kval[k] =
                          LTFAT_NAME_COMPLEX(malloc)( ktmp->size.height * ktmp->size.width));
    }
    else
    {
        CHECKMEM( ktmp->kval[0] =
                      LTFAT_NAME_COMPLEX(malloc)( ktmp->size.height * ktmp->size.width));
    }

    // Copy the zero-th kernel
    for (ltfat_int n = 0; n < ktmp->size.width; n++)
        memcpy(ktmp->kval[0] + n * ktmp->size.height,
               kernlarge + n * Mmax,
               ktmp->size.height * sizeof * kernlarge);

    if (do_allmods)
    {
        // Compute kernel modulations
        if (ptype == LTFAT_FREQINV)
            for (ltfat_int n = 1; n < ktmp->kNo; n++)
                LTFAT_NAME(dgtrealmp_kernel_modfi)(
                    ktmp->kval[0], ktmp->size, ktmp->mid,
                    n * kernskip , amin, Mmax, ktmp->kval[n]);
        else if (ptype == LTFAT_TIMEINV)
            for (ltfat_int m = 1; m < ktmp->kNo; m++)
                LTFAT_NAME(dgtrealmp_kernel_modti)(
                    ktmp->kval[0], ktmp->size, ktmp->mid,
                    m * kernskip, amin, Mmax, ktmp->kval[m]);
    }
    else
    {

    }

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
LTFAT_NAME(dgtrealmp_kernel_modfi)(
    LTFAT_COMPLEX* kfirst, ksize size, kanchor mid, ltfat_int n, ltfat_int a,
    ltfat_int M, LTFAT_COMPLEX* kmod)
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
LTFAT_NAME(dgtrealmp_kernel_modti)(
    LTFAT_COMPLEX* kfirst, ksize size, kanchor mid, ltfat_int m, ltfat_int a,
    ltfat_int M, LTFAT_COMPLEX* kmod)
{
    for (ltfat_int nn = 0; nn < size.width; nn++)
    {
        const LTFAT_COMPLEX* kfirstCol = kfirst + nn * size.height;
        LTFAT_COMPLEX* kmodCol         = kmod + nn * size.height;

        ltfat_int xval = nn - mid.wmid;
        LTFAT_REAL exparg = 2.0 * M_PI * m * xval * a / ((double) M);
        LTFAT_COMPLEX expmul = exp( I * exparg );

        for (ltfat_int mm = 0; mm < size.height; mm++ )
            kmodCol[mm] = expmul * kfirstCol[mm];
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
    CHECK(LTFATERR_NOTPOSARG, maxatoms > 0, "maxatoms must be greater than 0");
    p->params->maxatoms = maxatoms;
    p->params->maxit = 2 * maxatoms;
error:
    return status;
}

LTFAT_API int
LTFAT_NAME(dgtrealmp_set_errtoldb)(LTFAT_NAME(dgtrealmp_plan)* p,
                                   double errtoldb)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);
    CHECK(LTFATERR_BADARG, errtoldb <= 0, "errtoldb must be lower than 0");
    p->params->errtoldb = errtoldb;

error:
    return status;
}



LTFAT_API LTFAT_NAME(dgtreal_plan)**
LTFAT_NAME(dgtrealmp_getdgtrealplan)(LTFAT_NAME(dgtrealmp_plan)* p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);

    return p->dgtplans;
error:
    return NULL;

}

LTFAT_API int
LTFAT_NAME(dgtrealmp_getresidualcoef_compact)(
    LTFAT_NAME(dgtrealmp_plan)* p, LTFAT_COMPLEX* c)
{

    int status = LTFATERR_SUCCESS;
    CHECKNULL(p); CHECKNULL(c);

    ltfat_int accum = 0;
    for (ltfat_int k = 0; k < p->P; k++)
    {
        ltfat_int L = p->N[k] * p->M2[k];
        memcpy(c + accum, p->iterstate->c[k], L * sizeof * c);
        accum += L;
    }

error:
    return status;
}

/* LTFAT_API LTFAT_COMPLEX** */
/* LTFAT_NAME(dgtrealmp_getresidualcoef)(LTFAT_NAME(dgtrealmp_plan)* p) */
/* { */
/*     int status = LTFATERR_SUCCESS; */
/*     CHECKNULL(p); */
/*  */
/*     return p->iterstate->c; */
/* error: */
/*     return NULL; */
/*  */
/* } */

LTFAT_API ltfat_dgtrealmp_params*
LTFAT_NAME(dgtrealmp_getparams)(LTFAT_NAME(dgtrealmp_plan)* p)
{
    int status = LTFATERR_SUCCESS;
    CHECKNULL(p);

    return p->params;
error:
    return NULL;
}
