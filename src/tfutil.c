#include "config.h"

int LTFAT_NAME(complexprod)(LTFAT_COMPLEX *c, const LTFAT_COMPLEX a, const LTFAT_COMPLEX b)
{
#ifdef HAVE_COMPLEX_H
  (*c)=a*b;
#else
  
  (*c)[0]=a[0]*b[0]-a[1]*b[1];
  (*c)[1]=a[1]*b[0]+a[0]*b[1];

#endif
  return (0);
}


/* --- usefull for debugging --------
void print_z(const int N, LTFAT_COMPLEX *p)
{
  int i;
  
  for(i=0;i<N;i++)
  {
#ifdef HAVE_COMPLEX_H
    printf("%.16lf  %.16lf\n",creal(p[i]),cimag(p[i]));
#else
    printf("%.16lf  %.16lf\n",p[i][0],p[i][1]);
#endif
  }
  fflush(NULL);
  
}

*/
