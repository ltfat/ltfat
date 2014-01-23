/* NOT PROCESSED DIRECTLY, see ltfat_complexindependent.c */
#ifdef LTFAT_TYPE
#include "ltfat.h"

LTFAT_EXTERN void
LTFAT_NAME(col2diag)(const LTFAT_TYPE *cin, const ltfatInt L,
                     LTFAT_TYPE *cout)
{
    ltfatInt ii;

    LTFAT_TYPE *pcout;
    const LTFAT_TYPE *pcin;

    pcout=cout;
    const ltfatInt Lp1=L+1;
    for (ltfatInt jj=0; jj<L; jj++)
    {
        pcin=cin+(L-jj)*L;
        for (ii=0; ii<jj; ii++)
        {
            (*pcout) = (*pcin);
            pcout++;
            pcin+=Lp1;
        }
        pcin-=L*L;
        for (ii=jj; ii<L; ii++)
        {
            (*pcout) = (*pcin);
            pcout++;
            pcin+=Lp1;
        }
    }

}

#endif
