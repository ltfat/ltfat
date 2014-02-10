/* NOT PROCESSED DIRECTLY, see ltfat_complexindependent.c */
#ifdef LTFAT_TYPE

#ifndef GGA_UNROLL
#   define GGA_UNROLL 8
#endif

struct LTFAT_NAME(gga_plan_struct)
{
    const LTFAT_REAL* cos_term;
    const LTFAT_COMPLEX* cc_term;
    const LTFAT_COMPLEX* cc2_term;
    const ltfatInt M;
    const ltfatInt L;
};

struct LTFAT_NAME(chzt_plan_struct)
{
    LTFAT_COMPLEX* fbuffer;
    LTFAT_COMPLEX* W2;
    LTFAT_COMPLEX* Wo;
    LTFAT_COMPLEX* chirpF;
    const LTFAT_FFTW(plan) plan;
    const LTFAT_FFTW(plan) plan2;
    const ltfatInt L;
    const ltfatInt K;
    const ltfatInt Lfft;
};


LTFAT_EXTERN LTFAT_NAME(gga_plan)
LTFAT_NAME(gga_init)(const LTFAT_REAL *indVecPtr, const ltfatInt M, const ltfatInt L)
{
    LTFAT_REAL* cos_term = ltfat_malloc(M*sizeof*cos_term);
    LTFAT_COMPLEX* cc_term = ltfat_malloc(M*sizeof*cc_term);
    LTFAT_COMPLEX* cc2_term = ltfat_malloc(M*sizeof*cc2_term);

    LTFAT_REAL pik_term_pre = (LTFAT_REAL) (2.0*PI/((double) L));
    LTFAT_REAL _Complex cc2_pre = (-1.0*_Complex_I*((double)(L-1)));
    LTFAT_REAL _Complex cc_pre =  (-1.0*_Complex_I*((double)(L)));

    for(ltfatInt m=0; m<M; m++)
    {
        LTFAT_REAL pik_term = pik_term_pre*indVecPtr[m];
        cos_term[m] = (LTFAT_REAL) cos(pik_term)*2.0;
        cc_term[m] = (LTFAT_COMPLEX) cexp(cc_pre*pik_term);
        cc2_term[m] = (LTFAT_COMPLEX) cexp(cc2_pre*pik_term);
    }

    // This is workaround for defining constant elements of the struct.
    struct LTFAT_NAME(gga_plan_struct) plan_tmp =
    {.cos_term=cos_term,.cc_term=cc_term,.cc2_term=cc2_term,.M=M,.L=L};

    LTFAT_NAME(gga_plan) plan = ltfat_malloc(sizeof*plan);
    memcpy(plan,&plan_tmp,sizeof*plan);

    return plan;
}

LTFAT_EXTERN
void LTFAT_NAME(gga_done)(LTFAT_NAME(gga_plan) plan)
{
    LTFAT_SAFEFREEALL((void*)plan->cos_term,
                      (void*)plan->cc_term,
                      (void*)plan->cc2_term);
    ltfat_free(plan);
}


LTFAT_EXTERN
void LTFAT_NAME(gga_execute)(LTFAT_NAME(gga_plan) p,
                             const LTFAT_TYPE *fPtr,
                             const ltfatInt W,
                             LTFAT_COMPLEX *cPtr)
{
#ifndef GGA_UNROLL

    for(ltfatInt w=0; w<W; w++)
    {
        LTFAT_COMPLEX *cPtrTmp = (LTFAT_COMPLEX*) cPtr+w*p.M;

        for(ltfatInt m=0; m<p.M; m++)
        {
            LTFAT_TYPE s0 =  0.0;
            LTFAT_TYPE s1 =  0.0;
            LTFAT_TYPE s2 =  0.0;
            LTFAT_TYPE *fPtrTmp = (LTFAT_TYPE*) fPtr+w*p.L;

            for(ltfatInt ii=0; ii<p.L-1; ii++)
            {
                s0 = *fPtrTmp++ + p.cos_term[m]*s1 - s2;
                s2=s1;
                s1=s0;
            }
            s0 = *fPtrTmp + p.cos_term[m]*s1 - s2;

            *cPtrTmp++ = (s0*p.cc2_term[m] - s1*p.cc_term[m]);
        }
    }
#else
    for(ltfatInt w=0; w<W; w++)
    {
        LTFAT_COMPLEX *cPtrTmp = (LTFAT_COMPLEX*) cPtr+w*p->M;
        ltfatInt unrollRem = p->M%GGA_UNROLL;

        const LTFAT_REAL* cos_term = p->cos_term;
        const LTFAT_COMPLEX* cc_term = p->cc_term;
        const LTFAT_COMPLEX* cc2_term = p->cc2_term;
//#pragma omp parallel for
        for(ltfatInt m=0; m<p->M-unrollRem; m+=GGA_UNROLL)
        {
            LTFAT_TYPE s0[GGA_UNROLL];
            LTFAT_TYPE s1[GGA_UNROLL];
            LTFAT_TYPE s2[GGA_UNROLL];

            for(ltfatInt un=0; un<GGA_UNROLL; un++)
            {
                s0[un] = 0.0;
                s1[un] = 0.0;
                s2[un] = 0.0;
            }

            LTFAT_TYPE *fPtrTmp = (LTFAT_TYPE*) fPtr+w*p->L;

            for(ltfatInt ii=0; ii<p->L-1; ii++)
            {
                for(ltfatInt un=0; un<GGA_UNROLL; un++)
                {
                    s0[un] = *fPtrTmp + cos_term[un]*s1[un] - s2[un];
                    s2[un]=s1[un];
                    s1[un]=s0[un];
                }
                fPtrTmp++;
            }
            for(ltfatInt un=0; un<GGA_UNROLL; un++)
            {
                s0[un] = *fPtrTmp + cos_term[un]*s1[un] - s2[un];
                cPtrTmp[m+un] = (s0[un]*cc2_term[un] - s1[un]*cc_term[un]);
            }
            cos_term+=GGA_UNROLL;
            cc_term+=GGA_UNROLL;
            cc2_term+=GGA_UNROLL;
        }

        ltfatInt m= p->M-unrollRem;


        LTFAT_TYPE s0[GGA_UNROLL];
        LTFAT_TYPE s1[GGA_UNROLL];
        LTFAT_TYPE s2[GGA_UNROLL];

        for(ltfatInt un=0; un<unrollRem; un++)
        {
            s0[un] = 0.0;
            s1[un] = 0.0;
            s2[un] = 0.0;
        }

        LTFAT_TYPE *fPtrTmp = (LTFAT_TYPE*) fPtr+w*p->L;

        for(ltfatInt ii=0; ii<p->L-1; ii++)
        {
            for(ltfatInt un=0; un<unrollRem; un++)
            {
                s0[un] = *fPtrTmp + cos_term[un]*s1[un] - s2[un];
                s2[un]=s1[un];
                s1[un]=s0[un];
            }
            fPtrTmp++;
        }

        for(ltfatInt un=0; un<unrollRem; un++)
        {
            s0[un] = *fPtrTmp + cos_term[un]*s1[un] - s2[un];
            cPtrTmp[m+un] = (s0[un]*cc2_term[un] - s1[un]*cc_term[un]);
        }

    }
#endif
}



LTFAT_EXTERN
void LTFAT_NAME(gga)(const LTFAT_TYPE *fPtr, const LTFAT_REAL *indVecPtr,
                     const ltfatInt L, const ltfatInt W, const ltfatInt M, LTFAT_COMPLEX *cPtr)
{
    LTFAT_NAME(gga_plan) p = LTFAT_NAME(gga_init)(indVecPtr, M, L);
    LTFAT_NAME(gga_execute)(p,fPtr,W,cPtr);
    LTFAT_NAME(gga_done)(p);
}



LTFAT_EXTERN void
LTFAT_NAME(chzt)(const LTFAT_TYPE *fPtr, const ltfatInt L, const ltfatInt W,
                 const ltfatInt K, const LTFAT_REAL deltao, const LTFAT_REAL o,
                 LTFAT_COMPLEX *cPtr)
{
    LTFAT_NAME(chzt_plan) p = LTFAT_NAME(chzt_init)(K,L,deltao, o,
                              FFTW_ESTIMATE,
                              CZT_NEXTFASTFFT);

    LTFAT_NAME(chzt_execute)(p, fPtr, W, cPtr);

    LTFAT_NAME(chzt_done)(p);
}


LTFAT_EXTERN void
LTFAT_NAME(chzt_execute)(LTFAT_NAME(chzt_plan) p, const LTFAT_TYPE *fPtr,
                         const ltfatInt W, LTFAT_COMPLEX *cPtr)
{

    ltfatInt L = p->L;
    ltfatInt K = p->K;
    ltfatInt Lfft = p->Lfft;
    LTFAT_COMPLEX* fbuffer = p->fbuffer;
    LTFAT_FFTW(plan) plan_f = p->plan;
    LTFAT_FFTW(plan) plan_fi = p->plan2;
    LTFAT_COMPLEX* W2 = p->W2;
    LTFAT_COMPLEX* Wo = p->Wo;
    LTFAT_COMPLEX* chirpF = p->chirpF;

    for(ltfatInt w = 0; w<W; w++)
    {
        memset(fbuffer,0,Lfft*sizeof*fbuffer);
        LTFAT_NAME(array2complex)((LTFAT_TYPE *)(fPtr+w*L),fbuffer,L);

        //1) Premultiply by a chirp

        for(ltfatInt ii=0; ii<L; ii++)
        {
            fbuffer[ii]*=Wo[ii];
        }

        // 2) FFT of input
        LTFAT_FFTW(execute)(plan_f);

        // Frequency domain filtering
        for(ltfatInt ii=0; ii<Lfft; ii++)
        {
            fbuffer[ii]*=chirpF[ii];
        }


        // Inverse FFT
        LTFAT_COMPLEX *fPtrTmp = fbuffer;
        LTFAT_FFTW(execute)(plan_fi);

        // Final chirp multiplication and normalization
        LTFAT_COMPLEX *cPtrTmp = cPtr + w*K;
        for(ltfatInt ii=0; ii<K; ii++)
        {
            cPtrTmp[ii]=fPtrTmp[ii]*W2[ii];
        }

    }

}

LTFAT_EXTERN LTFAT_NAME(chzt_plan)
LTFAT_NAME(chzt_init)(const ltfatInt K, ltfatInt L, const LTFAT_REAL deltao,
                      const LTFAT_REAL o, const unsigned fftw_flags,
                      czt_ffthint hint)
{
    ltfatInt Lfft = L+K-1;

    if(hint == CZT_NEXTPOW2)
        Lfft = nextPow2(Lfft);
    else
        Lfft = nextfastfft(Lfft);

    LTFAT_COMPLEX* fbuffer = ltfat_malloc(Lfft*sizeof*fbuffer);
    LTFAT_FFTW(plan) plan_f =  LTFAT_FFTW(plan_dft_1d)(Lfft, fbuffer, fbuffer,
                               FFTW_FORWARD, fftw_flags);
    LTFAT_FFTW(plan) plan_fi =  LTFAT_FFTW(plan_dft_1d)(Lfft, fbuffer, fbuffer,
                                FFTW_BACKWARD, fftw_flags);

    // Pre and post chirp
    ltfatInt N = L>K?L:K;
    LTFAT_COMPLEX* W2 = ltfat_malloc(Lfft*sizeof*W2);
    LTFAT_COMPLEX* chirpF = ltfat_malloc(Lfft*sizeof*chirpF);
    LTFAT_COMPLEX* Wo = ltfat_malloc(L*sizeof*Wo);


    for(ltfatInt ii=0; ii<N; ii++)
    {
        W2[ii] = (LTFAT_COMPLEX) cexp(-1.0*I*deltao*ii*ii/2.0);
    }

    for(ltfatInt ii=0; ii<L; ii++)
    {
        Wo[ii] = (LTFAT_COMPLEX) cexp(-1.0*I*o*ii)*W2[ii];
    }
    // Set the rest to zero
    memset(W2+N,0,(Lfft-N)*sizeof*W2);

    LTFAT_NAME_COMPLEX(conjugate_array)(W2,chirpF,K);
    LTFAT_NAME_COMPLEX(conjugate_array)(W2+1,chirpF+Lfft-L+1,L-1);
    LTFAT_NAME_COMPLEX(reverse_array)(chirpF+Lfft-L+1,chirpF+Lfft-L+1,L-1);
    memset(chirpF+K,0,(Lfft-(L+K-1))*sizeof*chirpF);

    LTFAT_FFTW(execute_dft)(plan_f,chirpF,chirpF);


    for(ltfatInt ii=0; ii<K; ii++)
    {
        W2[ii] = (LTFAT_COMPLEX) cexp(-1.0*I*deltao*ii*ii/2.0)/((LTFAT_REAL) Lfft);
    }



    /*
    We could have shrinked the W2 to length K here.
    */
    struct LTFAT_NAME(chzt_plan_struct) p_struct =
    {
        .fbuffer = fbuffer, .plan=plan_f,.plan2=plan_fi,.L=L, .K=K, .W2 = W2,
        .Wo=Wo,.chirpF=chirpF,.Lfft=Lfft
    };

    LTFAT_NAME(chzt_plan) p = ltfat_malloc(sizeof*p);
    memcpy(p,&p_struct,sizeof*p);

    return  p;
}

LTFAT_EXTERN
void LTFAT_NAME(chzt_done)(LTFAT_NAME(chzt_plan) p)
{
    LTFAT_SAFEFREEALL(p->fbuffer,p->W2,p->Wo,p->chirpF);
    LTFAT_FFTW(destroy_plan)(p->plan);
    LTFAT_FFTW(destroy_plan)(p->plan2);
    ltfat_free(p);
}


LTFAT_EXTERN void
LTFAT_NAME(chzt_fac)(const LTFAT_TYPE *fPtr, const ltfatInt L,
                     const ltfatInt W, const ltfatInt K, const LTFAT_REAL deltao,
                     const LTFAT_REAL o, LTFAT_COMPLEX *cPtr)
{
    LTFAT_NAME(chzt_plan) p = LTFAT_NAME(chzt_fac_init)(K,L,deltao, o,FFTW_ESTIMATE, CZT_NEXTFASTFFT);

    LTFAT_NAME(chzt_fac_execute)(p, fPtr, W, cPtr);

    LTFAT_NAME(chzt_done)(p);
}

LTFAT_EXTERN void
LTFAT_NAME(chzt_fac_execute)(LTFAT_NAME(chzt_plan) p, const LTFAT_TYPE *fPtr,
                             const ltfatInt W, LTFAT_COMPLEX *cPtr)
{
    ltfatInt L = p->L;
    ltfatInt K = p->K;
    ltfatInt Lfft = p->Lfft;
    LTFAT_COMPLEX* fbuffer = p->fbuffer;
    LTFAT_FFTW(plan) plan_f = p->plan;
    LTFAT_FFTW(plan) plan_fi = p->plan2;
    LTFAT_COMPLEX* W2 = p->W2;
    LTFAT_COMPLEX* Wo = p->Wo;
    LTFAT_COMPLEX* chirpF = p->chirpF;

    LTFAT_COMPLEX* fBufTmp;
    ltfatInt q = (ltfatInt) ceil(((double)L)/((double)K));

    ltfatInt lastK = (L/q);

    for(ltfatInt w = 0; w<W; w++)
    {
        // *********************************
        // 1) Read and reorganize input data
        // *********************************
        memset(fbuffer,0,q*Lfft*sizeof*fbuffer);
        LTFAT_TYPE* fPtrTmp = ((LTFAT_TYPE *)fPtr) + w*L;

        for(ltfatInt k=0; k<lastK; k++)
        {
            LTFAT_TYPE* fTmp = fPtrTmp + k*q;
            fBufTmp = fbuffer + k;
            for(ltfatInt jj=0; jj<q; jj++)
            {
                *fBufTmp = (LTFAT_COMPLEX) fTmp[jj];
                fBufTmp+=Lfft;
            }
        }

        LTFAT_TYPE* fTmp = fPtrTmp + lastK*q;
        fBufTmp = fbuffer + lastK;
        for(ltfatInt jj=0; jj<L-lastK*q; jj++)
        {
            *fBufTmp = (LTFAT_COMPLEX) fTmp[jj];
            fBufTmp+=Lfft;
        }

        // *********************************
        // 2) Premultiply
        // *********************************

        fBufTmp = fbuffer;
        for(ltfatInt jj=0; jj<q; jj++)
        {
            for(ltfatInt k=0; k<K; k++)
            {
                fBufTmp[k]*=W2[k];
            }
            fBufTmp+=Lfft;
        }

        // *********************************
        // 3) q ffts of length Lfft
        // *********************************
        LTFAT_FFTW(execute)(plan_f);

        // *********************************
        // 4) Filter
        // *********************************
        fBufTmp = fbuffer;
        // Frequency domain filtering
        for(ltfatInt jj=0; jj<q; jj++)
        {
            for(ltfatInt ii=0; ii<Lfft; ii++)
            {
                fBufTmp[ii]*=chirpF[ii];
            }
            fBufTmp+=Lfft;
        }

        // *********************************
        // 5) q iffts of length Lfft
        // *********************************
        LTFAT_FFTW(execute)(plan_fi);


        // *********************************
        // 6) Postmultiply
        // *********************************
        fBufTmp = fbuffer;
        LTFAT_COMPLEX* Wotmp = Wo;
        for(ltfatInt jj=0; jj<q; jj++)
        {
            for(ltfatInt k=0; k<K; k++)
            {
                fBufTmp[k]*=Wotmp[k];
            }
            fBufTmp+=Lfft;
            Wotmp+=K;
        }

        // *********************************
        // 7) Sum cols
        // *********************************
        LTFAT_COMPLEX *cPtrTmp = cPtr + w*K;
        for(ltfatInt k=0; k<K; k++)
        {
            fBufTmp = fbuffer + k;
            cPtrTmp[k] = (LTFAT_COMPLEX) 0.0;
            for(ltfatInt jj=0; jj<q; jj++)
            {
                cPtrTmp[k]+= *fBufTmp;
                fBufTmp+=Lfft;
            }
        }

    }
}

LTFAT_EXTERN LTFAT_NAME(chzt_plan)
LTFAT_NAME(chzt_fac_init)(const ltfatInt K, const ltfatInt L,
                          const LTFAT_REAL deltao, const LTFAT_REAL o,
                          const unsigned fftw_flags, czt_ffthint hint)
{

    ltfatInt Lfft=2*K-1;
    if(hint == CZT_NEXTPOW2)
        Lfft = nextPow2(Lfft);
    else
        Lfft = nextfastfft(Lfft);

    ltfatInt q = (ltfatInt) ceil(((double)L)/((double)K));

    LTFAT_COMPLEX* fbuffer = ltfat_malloc(q*Lfft*sizeof*fbuffer);

    const LTFAT_FFTW(iodim) dims = {.n = Lfft, .is = 1, .os = 1};
    const LTFAT_FFTW(iodim) howmany_dims = {.n = q,.is = Lfft, .os = Lfft};
    LTFAT_FFTW(plan) plan_f =  LTFAT_FFTW(plan_guru_dft)(1, &dims, 1, &howmany_dims,
                               (LTFAT_COMPLEX*)fbuffer, (LTFAT_COMPLEX*) fbuffer,
                               FFTW_FORWARD, fftw_flags);

    LTFAT_FFTW(plan) plan_fi =  LTFAT_FFTW(plan_guru_dft)(1, &dims, 1, &howmany_dims,
                                (LTFAT_COMPLEX*)fbuffer, (LTFAT_COMPLEX*) fbuffer,
                                FFTW_BACKWARD, fftw_flags);

    LTFAT_COMPLEX* W2 = ltfat_malloc(K*sizeof*W2);
    LTFAT_COMPLEX* chirpF = ltfat_malloc(Lfft*sizeof*chirpF);
    LTFAT_COMPLEX* Wo = ltfat_malloc(q*K*sizeof*Wo);

    LTFAT_FFTW(plan) plan_chirpF =  LTFAT_FFTW(plan_dft_1d)(Lfft, (LTFAT_COMPLEX*)chirpF, (LTFAT_COMPLEX*) chirpF,
                                    FFTW_FORWARD, fftw_flags);

    for(ltfatInt k=0; k<K; k++)
    {
        W2[k] = (LTFAT_COMPLEX) cexp(-1.0*q*I*deltao*k*k/2.0);
    }

    LTFAT_NAME_COMPLEX(conjugate_array)(W2,chirpF,K);
    LTFAT_NAME_COMPLEX(conjugate_array)(W2+1,chirpF+Lfft-K+1,K-1);
    LTFAT_NAME_COMPLEX(reverse_array)(chirpF+Lfft-K+1,chirpF+Lfft-K+1,K-1);
    memset(chirpF+K,0,(Lfft-(2*K-1))*sizeof*chirpF);
    LTFAT_FFTW(execute)(plan_chirpF);
    LTFAT_FFTW(destroy_plan)(plan_chirpF);

    LTFAT_REAL oneoverLfft = 1.0/((LTFAT_REAL) Lfft);

    for(ltfatInt jj=0; jj<q; jj++)
    {
        LTFAT_COMPLEX* Wotmp = Wo + jj*K;
        for(ltfatInt k=0; k<K; k++)
        {
            Wotmp[k] = (LTFAT_COMPLEX) cexp(-1.0*I*jj*(k*deltao+o))*W2[k]*oneoverLfft;
        }
    }

    for(ltfatInt k=0; k<K; k++)
    {
        W2[k] *= (LTFAT_COMPLEX) cexp(-1.0*I*k*q*o);
    }


    struct LTFAT_NAME(chzt_plan_struct) p_struct =
    {.fbuffer = fbuffer, .plan=plan_f, .plan2=plan_fi,.L=L,.K=K, .W2 = W2,
     .Wo=Wo,.chirpF=chirpF,.Lfft=Lfft};

   LTFAT_NAME(chzt_plan) p = ltfat_malloc(sizeof*p);
   memcpy(p,&p_struct,sizeof*p);
   return  p;
}



#endif // LTFAT_TYPE

