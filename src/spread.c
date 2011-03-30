#include "config.h"
#include <math.h>
#include "ltfat.h"

LTFAT_EXTERN void
LTFAT_NAME(col2diag)(const LTFAT_COMPLEX *cin, const int L,
			  LTFAT_COMPLEX *cout)
{
   int ii;
   LTFAT_REAL *pcout;
   const LTFAT_REAL *pcin;

   pcout=(LTFAT_REAL*)cout;
   const int Lp1=2*(L+1);
   for (int jj=0;jj<L;jj++)
   {
      pcin=(const LTFAT_REAL*)cin+2*(L-jj)*L;
      for (ii=0;ii<jj;ii++)
      {
	 pcout[0] = pcin[0];
	 pcout[1] = pcin[1];
	 pcout+=2;
	 pcin+=Lp1;
      }
      pcin-=2*L*L;
      for (ii=jj;ii<L;ii++)
      {
	 pcout[0] = pcin[0];
	 pcout[1] = pcin[1];
	 pcout+=2;
	 pcin+=Lp1;
      }
   }   
}

LTFAT_EXTERN void
LTFAT_NAME(col2diag_r)(const LTFAT_REAL *cin, const int L,
			    LTFAT_REAL *cout)
{
   int ii;
 
   LTFAT_REAL *pcout;
   const LTFAT_REAL *pcin;
   
   pcout=cout;
   const int Lp1=L+1;
   for (int jj=0;jj<L;jj++)
   {
      pcin=cin+(L-jj)*L;
      for (ii=0;ii<jj;ii++)
      {
	 (*pcout) = (*pcin);
	 pcout++;
	 pcin+=Lp1;
      }
      pcin-=L*L;
      for (ii=jj;ii<L;ii++)
      {
	 (*pcout) = (*pcin);
	 pcout++;
	 pcin+=Lp1;
      }
   }
   
}
