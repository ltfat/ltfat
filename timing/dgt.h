#ifndef DGT_H
#define DGT_H 1

/* #include <common.h> */

#include "fftw3.h"


#define FFTW_OPTITYPE FFTW_ESTIMATE

/* BEGIN_C_DECLS */

/* FFTW_ESTIMATE is the only one to work, as the other types destroys input */

extern int gcd(const int a, const int b, int *r, int *s );

extern void wfac(fftw_complex *g, const int L, const int R,
	  const int a, const int M, fftw_complex *gf);

extern void wfac_r(double *g, const int L, const int R,
	  const int a, const int M, fftw_complex *gf);

extern void iwfac(fftw_complex *gf, const int L, const int R,
	   const int a, const int M, fftw_complex *g);

extern void iwfac_r(fftw_complex *gf, const int L, const int R,
	   const int a, const int M, double *g);

extern void dgt_fac(fftw_complex *f, fftw_complex *gf, const int L, const int W,
		    const int R, const int a, const int M, fftw_complex *cout, int dotime);

extern void idgt_fac(fftw_complex *c, fftw_complex *gf, const int L,
		    const int W,const int R,const int a, const int M,
		    fftw_complex *f);

extern void candual_fac(fftw_complex *g, const int L, const int R,
	  const int a, const int M, fftw_complex *gdualf);

extern void dgt_fb(fftw_complex *f, fftw_complex *g,
		   const int L, const int Lg,
		   const int W, const int R, const int a, const int M, 
		   fftw_complex *cout);

extern void idgt_fb(fftw_complex *cin, fftw_complex *g, const int L, const int gl,
		    const int W, const int R, const int a, const int M, 
		    fftw_complex *f);

/* This one is for debugging */
extern void print_z(const int N, fftw_complex *p);


/* END_C_DECLS */

#endif /* !DGT_H */
