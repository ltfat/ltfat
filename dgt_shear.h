typedef struct
{
   int a;
   int M;
   int L;
   int W;
   int s0;
   int s1;   
   int br;

   LTFAT_H_COMPLEX *p0;
   LTFAT_H_COMPLEX *p1;

   LTFAT_H_COMPLEX *fwork;
   LTFAT_H_COMPLEX *gwork;
   LTFAT_H_COMPLEX *c_rect;

   LTFAT_H_COMPLEX *finalmod;

   LTFAT_H_FFTW(plan) f_plan; 
   LTFAT_H_FFTW(plan) g_plan;


   LTFAT_H_NAME(dgt_long_plan) rect_plan;

   const LTFAT_H_COMPLEX *f;
   LTFAT_H_COMPLEX *cout;
   
} LTFAT_H_NAME(dgt_shear_plan);


LTFAT_EXTERN LTFAT_H_NAME(dgt_shear_plan)
LTFAT_H_NAME(dgt_shear_init)(
   const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *g,
   const int L, const int W, const int a,
   const int M, const int s0, const int s1, const int br,
   LTFAT_H_COMPLEX *cout,
   unsigned flags);
				  
LTFAT_EXTERN void 
LTFAT_H_NAME(dgt_shear_execute)(const LTFAT_H_NAME(dgt_shear_plan) plan);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_shear_done)(LTFAT_H_NAME(dgt_shear_plan) plan);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_shear)(const LTFAT_H_COMPLEX *f, const LTFAT_H_COMPLEX *g,
			const int L, const int W, const int a, const int M,
			const int s0, const int s1, const int br,  
			LTFAT_H_COMPLEX *c);

LTFAT_EXTERN void
LTFAT_H_NAME(pchirp)(const int L, const int n, LTFAT_H_COMPLEX *g);
