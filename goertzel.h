/*
Goertzel algorithm
*/
typedef struct
{
    const LTFAT_REAL* cos_term;
    const LTFAT_COMPLEX* cc_term;
    const LTFAT_COMPLEX* cc2_term;
    const size_t M;
    const size_t L;
} LTFAT_NAME(gga_plan);

LTFAT_EXTERN
LTFAT_NAME(gga_plan) LTFAT_NAME(create_gga_plan)(const double *indVecPtr, const size_t M, const size_t L);

LTFAT_EXTERN
void LTFAT_NAME(destroy_gga_plan)(LTFAT_NAME(gga_plan) plan);


LTFAT_EXTERN
void LTFAT_NAME(gga)(const LTFAT_TYPE *fPtr, const double *indVecPtr,
                  const size_t L, const size_t W, const size_t M,
		          LTFAT_COMPLEX *cPtr);

#define CZT_NEXTPOW2 0x0000
#define CZT_NEXTFASTFFT 0x0001

LTFAT_EXTERN
void LTFAT_NAME(gga_with_plan)(LTFAT_NAME(gga_plan) p,
                                 const LTFAT_TYPE *fPtr,
                                 LTFAT_COMPLEX *cPtr,
                                 const size_t W );

/*
Chirped Z transform
*/

typedef struct
{
    LTFAT_COMPLEX* fbuffer;
    LTFAT_COMPLEX* W2;
    LTFAT_COMPLEX* Wo;
    LTFAT_COMPLEX* chirpF;
    const LTFAT_FFTW(plan) plan;
    const LTFAT_FFTW(plan) plan2;
    size_t L;
    size_t K;
    size_t Lfft;
} LTFAT_NAME(chzt_plan);

LTFAT_EXTERN
void LTFAT_NAME(chzt)(const LTFAT_TYPE *fPtr, const size_t L, const size_t W, const size_t K,
                        const double deltao, const double o, LTFAT_COMPLEX *cPtr);

LTFAT_EXTERN
void LTFAT_NAME(chzt_with_plan)(LTFAT_NAME(chzt_plan) p, const LTFAT_TYPE *fPtr, const size_t W,
                         LTFAT_COMPLEX *cPtr);

LTFAT_EXTERN
LTFAT_NAME(chzt_plan) LTFAT_NAME(create_chzt_plan)(const size_t K, const size_t L, const double deltao, const double o, const unsigned fftw_flags, const unsigned czt_flags);

LTFAT_EXTERN
void LTFAT_NAME(destroy_chzt_plan)(LTFAT_NAME(chzt_plan) p);





LTFAT_EXTERN
void LTFAT_NAME(chzt_fact)(const LTFAT_TYPE *fPtr, const size_t L, const size_t W, const size_t K,
                        const double deltao, const double o, LTFAT_COMPLEX *cPtr);

LTFAT_EXTERN
void LTFAT_NAME(chzt_with_plan_fact)(LTFAT_NAME(chzt_plan) p, const LTFAT_TYPE *fPtr, const size_t W,
                         LTFAT_COMPLEX *cPtr);

LTFAT_EXTERN
LTFAT_NAME(chzt_plan) LTFAT_NAME(create_chzt_plan_fact)(const size_t K, const size_t L, const double deltao, const double o, const unsigned fftw_flags, const unsigned czt_flags);
