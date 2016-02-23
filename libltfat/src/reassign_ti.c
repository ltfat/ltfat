#include "ltfat.h"
#include "ltfat_types.h"


LTFAT_EXTERN void
LTFAT_NAME(filterbankphasegrad)(const LTFAT_COMPLEX* c [],
                                const LTFAT_COMPLEX* ch[],
                                const LTFAT_COMPLEX* cd[],
                                const ltfatInt          M,
                                const ltfatInt        N[],
                                const ltfatInt          L,
                                const LTFAT_REAL   minlvl,
                                LTFAT_REAL*        tgrad[],
                                LTFAT_REAL*        fgrad[],
                                LTFAT_REAL*           cs[])
{
#define FOREACHCOEF \
    for(ltfatInt m=0;m<M;++m){\
        for(ltfatInt ii=0;ii<N[m];++ii){

#define ARRAYEL(c) ((c)[m][ii])
#define ENDFOREACHCOEF }}

LTFAT_REAL minlvlAlt = LTFAT_COMPLEXH(cabs)(c[0][0]);

// Compute spectrogram from coefficients
// Keep max value
FOREACHCOEF
LTFAT_REAL en = LTFAT_COMPLEXH(cabs)(ARRAYEL(c))*LTFAT_COMPLEXH(cabs)(ARRAYEL(c));
ARRAYEL(cs) = en;
if(en>minlvlAlt)
    minlvlAlt = en;
ENDFOREACHCOEF

// Adjust minlvl 
minlvlAlt *= minlvl;

// Force spectrogram values less tha minLvlAlt to minlvlAlt
FOREACHCOEF
LTFAT_REAL csEl = ARRAYEL(cs);
if(csEl<minlvlAlt)
    ARRAYEL(cs) = minlvlAlt;
ENDFOREACHCOEF

// Instantaneous frequency
FOREACHCOEF
LTFAT_REAL tgradEl = LTFAT_COMPLEXH(creal)(
                         ARRAYEL(cd)*LTFAT_COMPLEXH(conj)(ARRAYEL(c))/ARRAYEL(cs)
                                          )/L*2;
ARRAYEL(tgrad) = fabs(tgradEl)<=2?tgradEl:0.0f;
ENDFOREACHCOEF


FOREACHCOEF
ARRAYEL(fgrad) = LTFAT_COMPLEXH(cimag)(
                        ARRAYEL(ch)*LTFAT_COMPLEXH(conj)(ARRAYEL(c))/ARRAYEL(cs)
                                      );
ENDFOREACHCOEF

#undef FOREACHCOEF
#undef ENDFOREACHCOEF
#undef ARRAYEL
}
