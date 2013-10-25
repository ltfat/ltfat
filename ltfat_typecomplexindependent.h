#include "wavelets.h"

LTFAT_EXTERN void
LTFAT_H_NAME(col2diag)(const LTFAT_H_TYPE *cin, const int L,
		       LTFAT_H_TYPE *cout);

LTFAT_EXTERN
void LTFAT_H_NAME(gabdual_long)(const LTFAT_H_TYPE *g,
				const int L, const int R, const int a,
				const int M, LTFAT_H_TYPE *gd);

LTFAT_EXTERN
void LTFAT_H_NAME(gabtight_long)(const LTFAT_H_TYPE *g,
				 const int L, const int R, const int a,
				 const int M, LTFAT_H_TYPE *gd);


/*
Goertzel algorithm
*/
typedef struct
{
    const LTFAT_H_REAL* cos_term;
    const LTFAT_H_COMPLEXH* cc_term;
    const LTFAT_H_COMPLEXH* cc2_term;
    const size_t M;
    const size_t L;
} LTFAT_H_NAME(gga_plan);

LTFAT_EXTERN
LTFAT_H_NAME(gga_plan) LTFAT_H_NAME(create_gga_plan)(const double *indVecPtr, const size_t M, const size_t L);

LTFAT_EXTERN
void LTFAT_H_NAME(destroy_gga_plan)(LTFAT_H_NAME(gga_plan) plan);


LTFAT_EXTERN
void LTFAT_H_NAME(gga)(const LTFAT_H_TYPE *fPtr, const double *indVecPtr,
                  const size_t L, const size_t W, const size_t M,
		          LTFAT_H_COMPLEXH *cPtr);

LTFAT_EXTERN
void LTFAT_H_NAME(gga_with_plan)(LTFAT_H_NAME(gga_plan) p,
                                 const LTFAT_H_TYPE *fPtr,
                                 LTFAT_H_COMPLEXH *cPtr,
                                 const size_t W );

/*
Chirped Z transform
*/

typedef struct
{
    LTFAT_H_COMPLEXH* fbuffer;
    LTFAT_H_COMPLEXH* W2;
    LTFAT_H_COMPLEXH* Wo;
    LTFAT_H_COMPLEXH* chirpF;
    const LTFAT_H_FFTW(plan) plan;
    const LTFAT_H_FFTW(plan) plan2;
    size_t L;
    size_t K;
    size_t Lfft;
} LTFAT_H_NAME(chzt_plan);

LTFAT_EXTERN
void LTFAT_H_NAME(chzt)(const LTFAT_H_TYPE *fPtr, const size_t L, const size_t W, const size_t K,
                        const double deltao, const double o, LTFAT_H_COMPLEXH *cPtr);

LTFAT_EXTERN
void LTFAT_H_NAME(chzt_with_plan)(LTFAT_H_NAME(chzt_plan) p, const LTFAT_H_TYPE *fPtr, const size_t W,
                        const double deltao, const double o, LTFAT_H_COMPLEXH *cPtr);

LTFAT_EXTERN
LTFAT_H_NAME(chzt_plan) LTFAT_H_NAME(create_chzt_plan)(const size_t K, const size_t L, const double deltao, const double o, const unsigned fftw_flags, const unsigned czt_flags);

LTFAT_EXTERN
void LTFAT_H_NAME(destroy_chzt_plan)(LTFAT_H_NAME(chzt_plan) p);





LTFAT_EXTERN
void LTFAT_H_NAME(chzt_fact)(const LTFAT_H_TYPE *fPtr, const size_t L, const size_t W, const size_t K,
                        const double deltao, const double o, LTFAT_H_COMPLEXH *cPtr);

LTFAT_EXTERN
void LTFAT_H_NAME(chzt_with_plan_fact)(LTFAT_H_NAME(chzt_plan) p, const LTFAT_H_TYPE *fPtr, const size_t W,
                        const double deltao, const double o, LTFAT_H_COMPLEXH *cPtr);

LTFAT_EXTERN
LTFAT_H_NAME(chzt_plan) LTFAT_H_NAME(create_chzt_plan_fact)(const size_t K, const size_t L, const double deltao, const double o, const unsigned fftw_flags, const unsigned czt_flags);




/*
Inplace friendly functions.
in==out is possibe.
*/

LTFAT_EXTERN
void LTFAT_H_NAME(circshift)(LTFAT_H_TYPE *in, LTFAT_H_TYPE *out, const ptrdiff_t L, const ptrdiff_t shift);

LTFAT_EXTERN
void LTFAT_H_NAME(reverse_array)(LTFAT_H_TYPE *in, LTFAT_H_TYPE *out, const size_t L);

LTFAT_EXTERN
void LTFAT_H_NAME(conjugate_array)(LTFAT_H_TYPE *in, LTFAT_H_TYPE *out, const size_t L);

LTFAT_EXTERN
void LTFAT_H_NAME(array2complex)(LTFAT_H_TYPE *in, LTFAT_H_COMPLEXH *out, const size_t L);






