/* NOT PROCESSED DIRECTLY, see ltfat_complexindependent.c */
#ifdef LTFAT_TYPE
#include <complex.h>
#include "config.h"
#include "ltfat.h"

LTFAT_EXTERN
void LTFAT_NAME(circshift)(LTFAT_TYPE *in, LTFAT_TYPE *out, const ptrdiff_t L, const ptrdiff_t shift)
{
   ptrdiff_t shiftMod = shift%L;

   if(in==out)
   {

       if(1)
       {
          LTFAT_TYPE *inTmp = (LTFAT_TYPE *)ltfat_malloc(L*sizeof(LTFAT_TYPE));
          memcpy(inTmp,in,L*sizeof(LTFAT_TYPE));
          LTFAT_NAME(circshift)(inTmp,out,L,shift);
          ltfat_free(inTmp);
       }
       else
       {
          int m,count,ii,jj;
          for(m=0,count=0;count!=L;m++)
          {
              LTFAT_TYPE t = in[m];
              for(ii=m,jj=m+shiftMod;
                  jj!=m;
                  ii=jj,jj=jj+shiftMod<L?jj+shiftMod:jj+shiftMod-L,count++)
              {
                  in[ii]=in[jj];
              }
              in[ii]=t;
              count++;
          }
       }



       return;
   }



   if(shiftMod<0)
   {
       memcpy(out,in-shiftMod,(L+shiftMod)*sizeof(LTFAT_TYPE));
       memcpy(out+(L+shiftMod),in,-shiftMod*sizeof(LTFAT_TYPE));
   }
   else if(shiftMod>0)
   {
       memcpy(out+shiftMod,in,(L-shiftMod)*sizeof(LTFAT_TYPE));
       memcpy(out,in+L-shiftMod,shiftMod*sizeof(LTFAT_TYPE));
   }
   else
   {
       memcpy(out,in,L*sizeof(LTFAT_TYPE));
   }
}

LTFAT_EXTERN
void LTFAT_NAME(reverse_array)(LTFAT_TYPE *in, LTFAT_TYPE *out,const size_t L)
{

   if(in==out)
   {
      LTFAT_TYPE tmpVar = (LTFAT_TYPE) 0.0;
      for(size_t ii=0;ii<L/2;ii++)
      {
        tmpVar = in[L-1-ii];
        in[L-1-ii] = in[ii];
        in[ii] = tmpVar;
      }
   }
   else
   {
      for(size_t ii=0;ii<L;ii++)
      {
        out[ii] = in[L-1-ii];
      }
   }
}

LTFAT_EXTERN
void LTFAT_NAME(conjugate_array)(LTFAT_TYPE *in, LTFAT_TYPE *out,const size_t L)
{
#ifdef LTFAT_COMPLEXTYPE
for(size_t ii=0;ii<L;ii++)
{
    out[ii] = LTFAT_COMPLEXH_NAME(conj)(in[ii]);
}
#else
if(in==out)
{
    return;
}
else
{
    memcpy(out,in,L*sizeof(LTFAT_TYPE));
}
#endif

}


LTFAT_EXTERN
void LTFAT_NAME(array2complex)(LTFAT_TYPE *in, LTFAT_COMPLEXH *out, const size_t L)
{
#ifdef LTFAT_COMPLEXTYPE
if(in==(LTFAT_TYPE*)out)
{
    return;
}
else
{
    memcpy(out,in,L*sizeof(LTFAT_COMPLEXH));
}
#else
if(in==(LTFAT_TYPE*)out)
{
    // This should produce an error
}
else
{
   LTFAT_REAL (*outTmp)[2] = (LTFAT_REAL(*)[2])  out;
   for(size_t ii=0;ii<L;ii++)
   {
      outTmp[ii][0] = in[ii];
      outTmp[ii][1] = (LTFAT_TYPE) 0.0;
   }
}
#endif
}


#endif // LTFAT_TYPE
