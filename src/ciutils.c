/* NOT PROCESSED DIRECTLY, see ltfat_complexindependent.c */
#ifdef LTFAT_TYPE
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
void LTFAT_NAME(reverse_array)(LTFAT_TYPE *in, LTFAT_TYPE *out, size_t L)
{

}


#endif // LTFAT_TYPE
