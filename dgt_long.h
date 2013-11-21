LTFAT_EXTERN void
LTFAT_NAME(dgt_long)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                      const int L, const int W,  const int a,
                      const int M, LTFAT_COMPLEX *cout);

typedef struct
{
  int a;
  int M;
  int L;
  int W;
  int c;
  int h_a;
  LTFAT_FFTW(plan) p_before;
  LTFAT_FFTW(plan) p_after;
  LTFAT_FFTW(plan) p_veryend;
  LTFAT_REAL *sbuf;
  const LTFAT_COMPLEX *f;
  LTFAT_COMPLEX *gf;
  LTFAT_COMPLEX *cout;
  LTFAT_REAL *ff, *cf;

} LTFAT_NAME(dgt_long_plan);


LTFAT_EXTERN LTFAT_NAME(dgt_long_plan)
LTFAT_NAME(dgt_long_init)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
                       const int L, const int W, const int a,
                       const int M, LTFAT_COMPLEX *cout,
                       unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgt_long_execute)(const LTFAT_NAME(dgt_long_plan) plan);

LTFAT_EXTERN void
LTFAT_NAME(dgt_long_done)(LTFAT_NAME(dgt_long_plan) plan);


LTFAT_EXTERN void
LTFAT_NAME(dgt_walnut_plan)(LTFAT_NAME(dgt_long_plan) plan);
