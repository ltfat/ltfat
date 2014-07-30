/* The 2nd row of this file will be used as additional parameters to gcc
../../thirdparty/Playrec/ltfatresample.c -lm -pedantic -std=c99
*/
#include "../../thirdparty/Playrec/ltfatresample.h"
#include "dbg.h"
#include "minunit.h"



void
fillrand(SAMPLE* a, size_t L)
{
   size_t m;
   srand (time(NULL));
   for(m=0;m<L;m++)
   {
      a[m]=((SAMPLE) (rand()))/RAND_MAX;
   }
}


char* test_filter()
{
   size_t L,ii;
   SAMPLE *in, *out, *out2,err;
   EMQFfilters ef;
   size_t bufNo = 10;
   size_t bufLen = 65;
   L = bufNo*bufLen;

   in = malloc(L * sizeof * in); fillrand(in,L);
   out = malloc(L * sizeof * out);
   out2 = malloc(L * sizeof * out);


   ef = emqffilters_init(0.1);

   /* Filter by blocks */
   for(ii=0;ii<bufNo;ii++)
   {
      emqffilters_dofilter(ef,in+ii*bufLen,bufLen,out+ii*bufLen);
   }
  
   emqffilters_done(&ef);
   ef = emqffilters_init(0.1);
   /* Filter */
   emqffilters_dofilter(ef,in,L,out2);


   emqffilters_done(&ef);

   err = 0;
   for(ii=0;ii<L;ii++)
   {
      err += abs(out[ii]-out2[ii]);
   }
   free(in);
   free(out);
   free(out2);

   mu_assert(err<1e-10,"FILT BY BLOCKS")

   return NULL;
}

char* test_filter_signal()
{
   size_t L,ii;
   SAMPLE *in, *out, *out2,err;
   EMQFfilters ef;
   size_t bufNo = 10;
   size_t bufLen = 65;
   L = bufNo*bufLen;

   in = malloc(L * sizeof * in); fillrand(in,L);
   out = malloc(L * sizeof * out);
   out2 = malloc(L * sizeof * out);


   ef = emqffilters_init(0.1);

   /* Filter by blocks */
   for(ii=0;ii<bufNo;ii++)
   {
      emqffilters_dofilter(ef,in+ii*bufLen,bufLen,out+ii*bufLen);
   }
  
   emqffilters_done(&ef);
   ef = emqffilters_init(0.1);
   /* Filter */
   emqffilters_dofilter(ef,in,L,out2);


   emqffilters_done(&ef);

   err = 0;
   for(ii=0;ii<L;ii++)
   {
      err += abs(out[ii]-out2[ii]);
   }
   free(in);
   free(out);
   free(out2);

   mu_assert(err<1e-10,"FILT BY BLOCKS")

   return NULL;

}

char *all_tests()
{
   mu_suite_start();

   mu_run_test(test_filter);
   mu_run_test(test_filter_signal);

   return NULL;
}

RUN_TESTS(all_tests)

