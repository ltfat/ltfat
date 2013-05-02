

typedef struct
{
   int a;
   int M;
   int L;
   int Lg;
   int W;
   int lt1;
   int lt2;

   LTFAT_H_COMPLEXH *f;
   LTFAT_H_COMPLEXH *c_scratch;
   LTFAT_H_COMPLEXH *cout;

   LTFAT_H_COMPLEXH *mwin;
   LTFAT_H_COMPLEXH *c_rect;

   LTFAT_H_COMPLEXH *mod;

   LTFAT_H_NAME(dgt_long_plan) *rect_plan_array;

} LTFAT_H_NAME(dgt_multi_plan);

LTFAT_EXTERN LTFAT_H_NAME(dgt_multi_plan)
LTFAT_H_NAME(dgt_multi_init)(const LTFAT_H_COMPLEXH *f, const LTFAT_H_COMPLEXH *g,
			     const int L, const int Lg, const int W, const int a, const int M,
			     const int lt1, const int lt2,
			     LTFAT_H_COMPLEXH *c,unsigned flags);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_multi_execute)(const LTFAT_H_NAME(dgt_multi_plan) plan);

LTFAT_EXTERN void
LTFAT_H_NAME(dgt_multi_done)(LTFAT_H_NAME(dgt_multi_plan) plan);

LTFAT_EXTERN void
LTFAT_H_NAME(nonsepwin2multi)(const LTFAT_H_COMPLEXH *g,
			      const int L, const int Lg, const int a, const int M,
			      const int lt1, const int lt2,
			      LTFAT_H_COMPLEXH *mwin);



LTFAT_EXTERN void
LTFAT_H_NAME(dgt_multi)(const LTFAT_H_COMPLEXH *f, const LTFAT_H_COMPLEXH *g,
			const int L, const int Lg, const int W, const int a, const int M,
			const int lt1, const int lt2,
			LTFAT_H_COMPLEXH *c);
