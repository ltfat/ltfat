LTFAT_EXTERN void 
LTFAT_H_NAME(dgt_long)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *g,
                      const int L, const int W,  const int a,
                      const int M, LTFAT_H_COMPLEX *cout);

typedef struct
{
  int a;
  int M;
  int L;
  int W;
  int c;
  int h_a;
  LTFAT_H_FFTW(plan) p_before; 
  LTFAT_H_FFTW(plan) p_after;
  LTFAT_H_FFTW(plan) p_veryend;
  LTFAT_H_REAL *sbuf;
  const LTFAT_H_COMPLEX *f;
  LTFAT_H_COMPLEX *gf;
  LTFAT_H_COMPLEX *cout;
  LTFAT_H_REAL *ff, *cf;

} LTFAT_H_NAME(dgt_long_plan);


LTFAT_EXTERN LTFAT_H_NAME(dgt_long_plan)
LTFAT_H_NAME(dgt_long_init)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *g,
                       const int L, const int W, const int a,
                       const int M, LTFAT_H_COMPLEX *cout,
                       unsigned flags);

LTFAT_EXTERN void 
LTFAT_H_NAME(dgt_long_execute)(const LTFAT_H_NAME(dgt_long_plan) plan);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_long_done)(LTFAT_H_NAME(dgt_long_plan) plan);


LTFAT_EXTERN void
LTFAT_H_NAME(dgt_walnut_plan)(LTFAT_H_NAME(dgt_long_plan) plan);
