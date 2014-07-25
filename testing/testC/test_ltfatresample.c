/* The 2nd row of this file will be used as additional parameters to gcc
  ../../thirdparty/Playrec/ltfatresample.c -lm
*/
#include "dbg.h"
#include "minunit.h"
#include "../../thirdparty/Playrec/ltfatresample.h"

char* test_coefload()
{
   size_t L;
   SAMPLE *in, *out;
   EMQFfilters ef;
   L = 2048;
   in = malloc(L * sizeof * in);
   out = malloc(L * sizeof * out);
   ef = emqffilters_init(0.1);
   mu_assert()

   emqffilters_dofilter(ef,in,L,out);

   

   emqffilters_done(&ef);
   free(in);
   free(out);
   return NULL;
}


char *all_tests()
{
   mu_suite_start();

   mu_run_test(test_coefload);

   return NULL;
}

RUN_TESTS(all_tests);

