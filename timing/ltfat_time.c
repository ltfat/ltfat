#include <stdlib.h>
#include <time.h>

/*
gettimeofday in MinGW has too coarse resolution. See
https://sourceforge.net/p/mingw/bugs/1821/
*/
#if defined(_WIN32) || defined(__WIN32__)
#include <windows.h>

double ltfat_time()
{
    LARGE_INTEGER frequency;
    LARGE_INTEGER t1;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&t1);

    return (t1.QuadPart) * 1000.0 / frequency.QuadPart;
}
#else
#if __STDC_VERSION__ >= 201112L
#include <sys/time.h>

double ltfat_time()
{
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);

    return (1000.0 * ((double)tv.tv_sec) + ((double)tv.tv_nsec) / 1000000.0);
}

#else
#include <sys/time.h>

double ltfat_time()
{
    struct timeval tv;

    gettimeofday(&tv, NULL);

    return (1000.0 * ((double)tv.tv_sec) + ((double)tv.tv_usec) / 1000.0);
}
#endif
#endif

/*
Fills array with pseudorandom values in range [0-1]
*/
void fillRand_d(double* f, int L)
{
    srand (time(NULL));

    for (int m = 0; m < L; m++)
    {
        f[m] = ((double) (rand())) / RAND_MAX;
    }
}

void fillRand_s(float* f, int L)
{
    srand (time(NULL));

    for (int m = 0; m < L; m++)
    {
        f[m] = ((float) (rand())) / RAND_MAX;
    }
}

void fillRand_cd(double _Complex* f, int L)
{
    double (*fTmp)[2]  = (double (*)[2]) f;
    srand (time(NULL));

    for (int m = 0; m < L; m++)
    {
        fTmp[m][0] = ((double) (rand())) / RAND_MAX;
        fTmp[m][1] = ((double) (rand())) / RAND_MAX;
    }
}

void fillRand_cs(float _Complex* f, int L)
{
    float (*fTmp)[2] = (float (*)[2]) f;
    srand (time(NULL));

    for (int m = 0; m < L; m++)
    {
        fTmp[m][0] = ((float) (rand())) / RAND_MAX;
        fTmp[m][1] = ((float) (rand())) / RAND_MAX;
    }
}


