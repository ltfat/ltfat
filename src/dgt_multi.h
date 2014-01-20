

typedef struct
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt L;
    ltfatInt Lg;
    ltfatInt W;
    ltfatInt lt1;
    ltfatInt lt2;

    LTFAT_COMPLEX *f;
    LTFAT_COMPLEX *c_scratch;
    LTFAT_COMPLEX *cout;

    LTFAT_COMPLEX *mwin;
    LTFAT_COMPLEX *c_rect;

    LTFAT_COMPLEX *mod;

    LTFAT_NAME(dgt_long_plan) *rect_plan_array;

} LTFAT_NAME(dgt_multi_plan);

LTFAT_EXTERN LTFAT_NAME(dgt_multi_plan)
LTFAT_NAME(dgt_multi_init)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                           const ltfatInt L, const ltfatInt Lg, const ltfatInt W, const ltfatInt a, const ltfatInt M,
                           const ltfatInt lt1, const ltfatInt lt2,
                           LTFAT_COMPLEX *c,unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgt_multi_execute)(const LTFAT_NAME(dgt_multi_plan) plan);

LTFAT_EXTERN void
LTFAT_NAME(dgt_multi_done)(LTFAT_NAME(dgt_multi_plan) plan);

LTFAT_EXTERN void
LTFAT_NAME(nonsepwin2multi)(const LTFAT_COMPLEX *g,
                            const ltfatInt L, const ltfatInt Lg, const ltfatInt a, const ltfatInt M,
                            const ltfatInt lt1, const ltfatInt lt2,
                            LTFAT_COMPLEX *mwin);



LTFAT_EXTERN void
LTFAT_NAME(dgt_multi)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                      const ltfatInt L, const ltfatInt Lg, const ltfatInt W, const ltfatInt a, const ltfatInt M,
                      const ltfatInt lt1, const ltfatInt lt2,
                      LTFAT_COMPLEX *c);
