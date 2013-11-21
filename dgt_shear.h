typedef struct
{
   int a;
   int M;
   int L;
   int W;
   int s0;
   int s1;
   int br;

   LTFAT_COMPLEX *p0;
   LTFAT_COMPLEX *p1;

   LTFAT_COMPLEX *fwork;
   LTFAT_COMPLEX *gwork;
   LTFAT_COMPLEX *c_rect;

   LTFAT_COMPLEX *finalmod;

   LTFAT_FFTW(plan) f_plan;
   LTFAT_FFTW(plan) g_plan;


   LTFAT_NAME(dgt_long_plan) rect_plan;

   const LTFAT_COMPLEX *f;
   LTFAT_COMPLEX *cout;

} LTFAT_NAME(dgt_shear_plan);


LTFAT_EXTERN LTFAT_NAME(dgt_shear_plan)
LTFAT_NAME(dgt_shear_init)(
   const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
   const int L, const int W, const int a,
   const int M, const int s0, const int s1, const int br,
   LTFAT_COMPLEX *cout,
   unsigned flags);

LTFAT_EXTERN void
LTFAT_NAME(dgt_shear_execute)(const LTFAT_NAME(dgt_shear_plan) plan);

LTFAT_EXTERN void
LTFAT_NAME(dgt_shear_done)(LTFAT_NAME(dgt_shear_plan) plan);

LTFAT_EXTERN void
LTFAT_NAME(dgt_shear)(const LTFAT_COMPLEX *f, const LTFAT_COMPLEX *g,
			const int L, const int W, const int a, const int M,
			const int s0, const int s1, const int br,
			LTFAT_COMPLEX *c);

LTFAT_EXTERN void
LTFAT_NAME(pchirp)(const long L, const long n, LTFAT_COMPLEX *g);
