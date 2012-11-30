#ifndef WAVELETS_H
#define WAVELETS_H

#include <string.h>
#include "config.h"
#include "fftw3.h"
#include "ltfat.h"


#define ONEOVERSQRT2 0.707106781186548 \

//#define WORDBITSM1 (sizeof(int)*8)-1

#define LTFAT_H_REAL double
#define LTFAT_H_COMPLEX ltfat_complex
#define LTFAT_H_NAME(name) name
#define LTFAT_H_FFTW(name) fftw_ ## name  

// common basic routines

LTFAT_EXTERN
void LTFAT_H_NAME(extend_left)(const LTFAT_H_REAL *in,int inLen, LTFAT_H_REAL *buffer, int buffLen, int filtLen, int type);

LTFAT_EXTERN
void LTFAT_H_NAME(extend_right)(const LTFAT_H_REAL *in,int inLen, LTFAT_H_REAL *buffer, int filtLen, int type);

LTFAT_EXTERN
void LTFAT_H_NAME(conv_td_sub)(const LTFAT_H_REAL *in, int inLen, LTFAT_H_REAL *out[], const int outLen, const LTFAT_H_REAL *filts[], int fLen, int noOfFilts, int sub, int skip, int ext, int filtUps);


LTFAT_EXTERN
void LTFAT_H_NAME(up_conv_td)(const LTFAT_H_REAL *in[], int inLen, LTFAT_H_REAL *out, const int outLen, const LTFAT_H_REAL *filts[], int fLen, int noOfFilts, int up, int skip, int ext, int filtUps);

// additional common basic routines

LTFAT_EXTERN
void LTFAT_H_NAME(up_conv_sub)(const LTFAT_H_REAL *in, int inLen, LTFAT_H_REAL *out, int outLen, LTFAT_H_REAL *filt, int fLen, int up, int sub, int skip, int ext);

LTFAT_EXTERN
void LTFAT_H_NAME(up_conv_sub_1toN)(const LTFAT_H_REAL *in, int inLen, LTFAT_H_REAL *out[], const int outLen, const LTFAT_H_REAL *filts[], int fLen, int noOfFilts, int sub, int skip, int ext);

LTFAT_EXTERN
void LTFAT_H_NAME(up_conv_sub_Nto1)(const LTFAT_H_REAL *in[], int inLen, LTFAT_H_REAL *out, const int outLen, const LTFAT_H_REAL *filts[], int fLen, int noOfFilts, int up, int skip, int ext);


// execution routines
// FORWARD TRANSFORMS
LTFAT_EXTERN
void LTFAT_H_NAME(dyadic_dwt_exp)(const LTFAT_H_REAL *in, int inLen, LTFAT_H_REAL *out[], const int outLen[], const LTFAT_H_REAL *filts[], int fLen, const int J, int ext);

LTFAT_EXTERN
void LTFAT_H_NAME(dyadic_dwt_per)(const LTFAT_H_REAL *in, int inLen, LTFAT_H_REAL *out[], const int outLen[], const LTFAT_H_REAL *filts[], int fLen, const int J);

LTFAT_EXTERN
void LTFAT_H_NAME(undec_dwt_exp)(const LTFAT_H_REAL *in, int inLen, LTFAT_H_REAL *out[], const int outLen[], const LTFAT_H_REAL *filts[], int fLen, const int J, int ext);

LTFAT_EXTERN
void LTFAT_H_NAME(undec_dwt_per)(const LTFAT_H_REAL *in, int inLen, LTFAT_H_REAL *out[], const int outLen, const LTFAT_H_REAL *filts[], int fLen, const int J);

// execution routines
// INVERSE TRANSFORMS
LTFAT_EXTERN
void LTFAT_H_NAME(dyadic_idwt_exp)(const LTFAT_H_REAL *in[], int inLen[], LTFAT_H_REAL *out, const int outLen, const LTFAT_H_REAL *filts[], int fLen, const int J);

LTFAT_EXTERN
void LTFAT_H_NAME(dyadic_idwt_per)(const LTFAT_H_REAL *in[], int inLen[], LTFAT_H_REAL *out, const int outLen, const LTFAT_H_REAL *filts[], int fLen, const int J);

LTFAT_EXTERN
void LTFAT_H_NAME(undec_idwt_exp)(const LTFAT_H_REAL *in[], int inLen[], LTFAT_H_REAL *out, const int outLen, const LTFAT_H_REAL *filts[], int fLen, const int J);

LTFAT_EXTERN
void LTFAT_H_NAME(undec_idwt_per)(const LTFAT_H_REAL *in[], int inLen, LTFAT_H_REAL *out, const int outLen, const LTFAT_H_REAL *filts[], int fLen, const int J);

#undef LTFAT_H_REAL
#undef LTFAT_H_COMPLEX
#undef LTFAT_H_NAME
#undef LTFAT_H_FFTW


static inline int pow2(int x){
	return ((1)<<(x));
}

static inline int modPow2(int x,int pow2var){
	return ((x)&(pow2var-1));
}

static inline int nextPow2(int x){
	x--;
	(x) = ((x)>>1)  | (x);   
    (x) = ((x)>>2)  | (x);   
	(x) = ((x)>>4)  | (x);   
	(x) = ((x)>>8)  | (x);   
	(x) = ((x)>>16) | (x);   
	(x)++;   
	return x;
}

static inline int imin(int x,int y){
	//return x-(((x-y)>>(WORDBITSM1))&(x-y));
	return x>y?y:x;
}

static inline int imax(int x,int y){
	return x>y?x:y;
}


#endif







