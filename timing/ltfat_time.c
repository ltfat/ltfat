#include <sys/time.h>
#include <stdlib.h>

double ltfat_time()
{
   struct timeval tv;
   
   gettimeofday(&tv, NULL);

   /*printf("  %i  %i  \n",tv.tv_sec,tv.tv_usec);*/
   
   return ((double)tv.tv_sec+(double)tv.tv_usec/1000000.0);
}
