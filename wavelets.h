#ifndef WAVELETS_H
#define WAVELETS_H
#include <string.h>

#define ONEOVERSQRT2 0.707106781186548 \

enum ltfatWavExtType
{
PER,
PERDEC,
PPD,
SYM,
EVEN,
SYMW,
ASYM,
ODD,
ASYMW,
SP0,
ZPD,
ZERO,
BAD_TYPE
};

inline static enum ltfatWavExtType ltfatExtStringToEnum(char* extType)
{
    if(strcmp(extType,"per")==0)
    {
       return PER;
    }
    else if(strcmp(extType,"perdec")==0)
    {
       return PERDEC;
    }
    else if(strcmp(extType,"ppd")==0)
    {
       return PPD;
    }
    else if(strcmp(extType,"sym")==0)
    {
       return SYM;
    }
    else if(strcmp(extType,"even")==0)
    {
       return EVEN;
    }
    else if(strcmp(extType,"symw")==0)
    {
       return SYMW;
    }
    else if(strcmp(extType,"odd")==0)
    {
       return ODD;
    }
    else if(strcmp(extType,"asymw")==0)
    {
       return ASYMW;
    }
    else if(strcmp(extType,"sp0")==0)
    {
       return SP0;
    }
    else if(strcmp(extType,"zpd")==0)
    {
       return ZPD;
    }
    else if(strcmp(extType,"zero")==0)
    {
       return ZERO;
    }
    else
    {
       return BAD_TYPE;
    }
}


// common basic routines

LTFAT_EXTERN
void LTFAT_H_NAME(extend_left)(const LTFAT_H_REAL *in,int inLen, LTFAT_H_REAL *buffer, int buffLen, int filtLen, enum ltfatWavExtType ext, int a);

LTFAT_EXTERN
void LTFAT_H_NAME(extend_right)(const LTFAT_H_REAL *in,int inLen, LTFAT_H_REAL *buffer, int filtLen, enum ltfatWavExtType ext, int a);

LTFAT_EXTERN
void LTFAT_H_NAME(convsub_td)(const LTFAT_H_REAL *in, int inLen, LTFAT_H_REAL *out, const int outLen, const LTFAT_H_REAL *filts, int fLen, int sub, int skip, enum ltfatWavExtType ext);

LTFAT_EXTERN
void LTFAT_H_NAME(conv_td_sub)(const LTFAT_H_REAL *in, int inLen, LTFAT_H_REAL *out[], const int outLen, const LTFAT_H_REAL *filts[], int fLen, int noOfFilts, int sub, int skip, enum ltfatWavExtType ext, int filtUps);


LTFAT_EXTERN
void LTFAT_H_NAME(up_conv_td)(LTFAT_H_REAL *in[], int inLen, LTFAT_H_REAL *out, const int outLen, const LTFAT_H_REAL *filts[], int fLen, int noOfFilts, int up, int skip, int ext, int filtUps);


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

static inline int isPow2(int x){
	return x==nextPow2(x);
}

static inline int ilog2(int x){
	int tmp = 0;
    while (x >>= 1) ++tmp;
	return tmp;
}

static inline int imin(int x,int y){
	//return x-(((x-y)>>(WORDBITSM1))&(x-y));
	return x>y?y:x;
}

static inline int imax(int x,int y){
	return x>y?x:y;
}

// integer power by squaring
static inline int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}


#endif







