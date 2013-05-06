typedef struct
{
   int a;
   int M;
   int L;
   int W;
   int s0;
   int s1;
   int br;

   LTFAT_H_COMPLEXH *p0;
   LTFAT_H_COMPLEXH *p1;

   LTFAT_H_COMPLEXH *fwork;
   LTFAT_H_COMPLEXH *gwork;
   LTFAT_H_COMPLEXH *c_rect;

   LTFAT_H_COMPLEXH *finalmod;

   LTFAT_H_FFTW(plan) f_plan;
   LTFAT_H_FFTW(plan) g_plan;


   LTFAT_H_NAME(dgt_long_plan) rect_plan;

   const LTFAT_H_COMPLEXH *f;
   LTFAT_H_COMPLEXH *cout;

} LTFAT_H_NAME(dgt_shear_plan);


LTFAT_EXTERN LTFAT_H_NAME(dgt_shear_plan)
LTFAT_H_NAME(dgt_shear_init)(
   const LTFAT_H_COMPLEXH *f, const LTFAT_H_COMPLEXH *g,
   const int L, const int W, const int a,
   const int M, const int s0, const int s1, const int br,
   LTFAT_H_COMPLEXH *cout,
   unsigned flags);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_shear_execute)(const LTFAT_H_NAME(dgt_shear_plan) plan);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_shear_done)(LTFAT_H_NAME(dgt_shear_plan) plan);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_shear)(const LTFAT_H_COMPLEXH *f, const LTFAT_H_COMPLEXH *g,
			const int L, const int W, const int a, const int M,
			const int s0, const int s1, const int br,
			LTFAT_H_COMPLEXH *c);

LTFAT_EXTERN void
LTFAT_H_NAME(pchirp)(const long L, const long n, LTFAT_H_COMPLEXH *g);
