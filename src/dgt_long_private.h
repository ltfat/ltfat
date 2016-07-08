#include "ltfat.h"
#include "ltfat/types.h"
#include "ltfat/macros.h"

struct LTFAT_NAME_REAL(dgt_long_plan)
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt L;
    ltfatInt W;
    ltfatInt c;
    ltfatInt h_a;
    dgt_phasetype ptype;
    LTFAT_FFTW(plan) p_before;
    LTFAT_FFTW(plan) p_after;
    LTFAT_FFTW(plan) p_veryend;
    LTFAT_REAL* sbuf;
    const LTFAT_REAL* f;
    LTFAT_COMPLEX* gf;
    LTFAT_COMPLEX* cout;
    LTFAT_REAL* ff, *cf;
};

struct LTFAT_NAME_COMPLEX(dgt_long_plan)
{
    ltfatInt a;
    ltfatInt M;
    ltfatInt L;
    ltfatInt W;
    ltfatInt c;
    ltfatInt h_a;
    dgt_phasetype ptype;
    LTFAT_FFTW(plan) p_before;
    LTFAT_FFTW(plan) p_after;
    LTFAT_FFTW(plan) p_veryend;
    LTFAT_REAL* sbuf;
    const LTFAT_COMPLEX* f;
    LTFAT_COMPLEX* gf;
    LTFAT_COMPLEX* cout;
    LTFAT_REAL* ff, *cf;
};

