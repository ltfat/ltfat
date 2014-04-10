typedef struct
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt L;
    ltfatInt W;
    ltfatInt c;
    ltfatInt h_a;
    dgt_phasetype ptype;
    LTFAT_FFTW(plan) p_before;
    LTFAT_FFTW(plan) p_after;
    LTFAT_FFTW(plan) p_veryend;
    LTFAT_REAL *sbuf;
    const LTFAT_COMPLEX *f;
    LTFAT_COMPLEX *gf;
    LTFAT_COMPLEX *cout;
    LTFAT_REAL *ff, *cf;

} LTFAT_NAME(dgt_long_plan);



LTFAT_EXTERN void
LTFAT_NAME_COMPLEX(dgt_long)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                             const ltfatInt L, const ltfatInt W,  const ltfatInt a,
                             const ltfatInt M, const dgt_phasetype ptype,
                             LTFAT_COMPLEX *cout);

LTFAT_EXTERN void
LTFAT_NAME(dgt_long)(const LTFAT_REAL *f, const LTFAT_REAL *g,
                     const ltfatInt L, const ltfatInt W,  const ltfatInt a,
                     const ltfatInt M, const dgt_phasetype ptype,
                     LTFAT_COMPLEX *cout);


LTFAT_EXTERN LTFAT_NAME(dgt_long_plan)
LTFAT_NAME(dgt_long_init)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                          const ltfatInt L, const ltfatInt W, const ltfatInt a,
                          const ltfatInt M, LTFAT_COMPLEX *cout,
                          const dgt_phasetype ptype, unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgt_long_execute)(const LTFAT_NAME(dgt_long_plan) plan);

LTFAT_EXTERN void
LTFAT_NAME(dgt_long_done)(LTFAT_NAME(dgt_long_plan) plan);


LTFAT_EXTERN void
LTFAT_NAME(dgt_walnut_plan)(LTFAT_NAME(dgt_long_plan) plan);
