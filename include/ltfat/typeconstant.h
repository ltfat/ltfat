#ifndef _LTFAT_TYPECONSTANT
#define _LTFAT_TYPECONSTANT
#include "memalloc.h"
#include "dgt_common.h"

typedef struct
{
    ltfatInt quot;
    ltfatInt rem;
} ltfat_div_t;

/* -------- Define routines that do not change between single/double-- */
LTFAT_EXTERN ltfat_div_t
ltfat_idiv(const ltfatInt a, const ltfatInt b);

LTFAT_EXTERN ltfatInt
ltfat_gcd(const ltfatInt a, const ltfatInt b, ltfatInt *r, ltfatInt *s );

LTFAT_EXTERN void
ltfat_fftindex(const ltfatInt N, ltfatInt *indexout);

LTFAT_EXTERN
ltfatInt makelarger(const ltfatInt L, const ltfatInt K);

LTFAT_EXTERN
ltfatInt ltfat_imax(const ltfatInt a, const ltfatInt b);

LTFAT_EXTERN
ltfatInt ltfat_imin(const ltfatInt a, const ltfatInt b);

/** \addtogroup utils
 * @{
 */
/** Find least common multiple of a and b
 */
LTFAT_EXTERN
ltfatInt ltfat_lcm(const ltfatInt a, const ltfatInt b);

/** Find next suitable L for signal length Ls and Gabor lattice parameters a and M
 */
LTFAT_EXTERN ltfatInt
ltfat_dgtlength(const ltfatInt Ls, const ltfatInt a, const ltfatInt M);
/** @}*/

LTFAT_EXTERN
void gabimagepars(const ltfatInt Ls, const ltfatInt x, const ltfatInt y,
                  ltfatInt *a, ltfatInt *M, ltfatInt *L, ltfatInt *N, ltfatInt *Ngood);

LTFAT_EXTERN
ltfatInt wfacreal_size(const ltfatInt L, const ltfatInt a, const ltfatInt M);

LTFAT_EXTERN ltfatInt
ltfat_nextfastfft(const ltfatInt x);

LTFAT_EXTERN ltfatInt
ltfat_pow2(const ltfatInt x);

LTFAT_EXTERN ltfatInt
ltfat_nextpow2(const ltfatInt x);

LTFAT_EXTERN ltfatInt
ltfat_modpow2(const ltfatInt x, const ltfatInt pow2var);

LTFAT_EXTERN ltfatInt
ltfat_round(const double x);

LTFAT_EXTERN ltfatInt
ltfat_positiverem(const ltfatInt a, const ltfatInt b);

LTFAT_EXTERN ltfatInt
ltfat_rangelimit(const ltfatInt a, const ltfatInt amin, const ltfatInt amax);


// Custom headers are down here
#include "reassign_typeconstant.h"

#endif /* _LTFAT_TYPECONSTANT */
