#if !(defined(LTFAT_DOUBLE) || defined(LTFAT_SINGLE))
#define LTFAT_DOUBLE
#endif

#include "../include/ltfat.h"
#include "../include/ltfat/types.h"
#include "../include/ltfat/macros.h"

#ifdef __cplusplus
#include <chrono>
#include <iostream>

using Clock = std::chrono::high_resolution_clock;
using namespace std;
#else
#endif
