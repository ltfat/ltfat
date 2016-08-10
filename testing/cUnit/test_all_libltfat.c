#include "ltfat.h"
#include "ltfat/errno.h"
#include "ltfat/macros.h"
#include "minunit.h"


void all_tests()
{
    mu_suite_start();

    mu_run_test_singledoublecomplex(test_circshift);
    mu_run_test_singledoublecomplex(test_fftshift);
    mu_run_test_singledoublecomplex(test_ifftshift);
    mu_run_test_singledoublecomplex(test_fir2long);
    mu_run_test_singledoublecomplex(test_long2fir);
    mu_run_test_singledoublecomplex(test_normalize);
    mu_run_test_singledouble(test_fftcircshift);
    mu_run_test_singledouble(test_fftfftshift);
    mu_run_test_singledouble(test_fftifftshift);
    mu_run_test_singledouble(test_fftrealcircshift);
    mu_run_test_singledouble(test_fftrealfftshift);
    mu_run_test_singledouble(test_fftrealifftshift);

    mu_suite_stop();
}

RUN_TESTS(all_tests)
