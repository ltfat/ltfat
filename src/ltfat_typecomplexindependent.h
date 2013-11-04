#include "wavelets.h"
#include "goertzel.h"
#include "ciutils.h"

LTFAT_EXTERN void
LTFAT_H_NAME(col2diag)(const LTFAT_H_TYPE *cin, const int L,
		       LTFAT_H_TYPE *cout);

LTFAT_EXTERN
void LTFAT_H_NAME(gabdual_long)(const LTFAT_H_TYPE *g,
				const int L, const int R, const int a,
				const int M, LTFAT_H_TYPE *gd);

LTFAT_EXTERN
void LTFAT_H_NAME(gabtight_long)(const LTFAT_H_TYPE *g,
				 const int L, const int R, const int a,
				 const int M, LTFAT_H_TYPE *gd);


