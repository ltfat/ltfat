#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"


struct LTFAT_NAME(dgtreal_long_plan)
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt L;
    ltfatInt W;
    ltfatInt c;
    ltfatInt h_a;
    ltfat_phaseconvention ptype;
    LTFAT_FFTW(plan) p_before;
    LTFAT_FFTW(plan) p_after;
    LTFAT_FFTW(plan) p_veryend;
    LTFAT_REAL* sbuf;
    LTFAT_COMPLEX* cbuf;
    const LTFAT_REAL* f;
    LTFAT_COMPLEX* gf;
    LTFAT_REAL* cwork;
    LTFAT_COMPLEX* cout;
    LTFAT_REAL* ff, *cf;
};
