#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_wfac // change to filename
#define OCTFILEHELP "Computes window factorization.\n\
                    Usage: c=comp_wfac(g,a,M);\n\
                    Yeah."

#include "ltfat_oct_template_helper.h"
// octave_idx_type 32 or 64 bit signed integer
/*
  dgt_fb forwarders
*/

static inline void fwd_comp_wfac(const Complex *g,
                              const octave_idx_type L, const octave_idx_type R,
		                      const octave_idx_type a, const octave_idx_type M,
							  Complex *cout)
{
   wfac_d(   reinterpret_cast<const double (*)[2]>(g),
             L,R,a,M,
			 reinterpret_cast<double (*)[2]>(cout));
}

static inline void fwd_comp_wfac(const FloatComplex *g,
                              const octave_idx_type L, const octave_idx_type R,
		                      const octave_idx_type a, const octave_idx_type M,
							  FloatComplex *cout)
{
   wfac_s(   reinterpret_cast<const float (*)[2]>(g),
             L,R,a,M,
			 reinterpret_cast<float (*)[2]>(cout));
}

static inline void fwd_comp_wfac(const double *g,
                              const octave_idx_type L, const octave_idx_type R,
		                      const octave_idx_type a, const octave_idx_type M,
							  Complex *cout)
{
   wfac_r_d(reinterpret_cast<const double*>(g),
            L,R,a,M,
			reinterpret_cast<double (*)[2]>(cout));
}

static inline void fwd_comp_wfac(const float *g,
                              const octave_idx_type L, const octave_idx_type R,
		                      const octave_idx_type a, const octave_idx_type M,
							  FloatComplex *cout)
{
   wfac_r_s(reinterpret_cast<const float*>(g),
            L,R,a,M,
			reinterpret_cast<float (*)[2]>(cout));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(0));
	 const octave_idx_type a = args(1).int_value();
     const octave_idx_type M = args(2).int_value();
	 const octave_idx_type L = g.rows();
     const octave_idx_type R = g.columns();
     
     const octave_idx_type b = L/M;
	 
     int h_a, h_m;    
     const octave_idx_type c=gcd(a, M,&h_a, &h_m);
     const octave_idx_type p=a/c;
     const octave_idx_type q=M/c;
     const octave_idx_type d=b/p;

     MArray<LTFAT_COMPLEX> cout(dim_vector(p*q*R,c*d)); 
     cout.fill(0);
	 
	 fwd_comp_wfac(g.data(),L,R,a,M,cout.fortran_vec());
	 
     return octave_value(cout);
}

/*
#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_wfac, args, ,
  "Computes window factorization.\n\
  Usage: c=comp_wfac(g,a,M);\n\
  Yeah.")
{

  const bool g_is_complex  = args(0).is_complex_type();

  const int a = args(1).int_value();
  const int M = args(2).int_value();

  int h_a, h_m;
   
  if (g_is_complex)
  {
     const ComplexMatrix g = args(0).complex_matrix_value();
     const int L = g.rows();
     const int R = g.columns();
     
     const int b = L/M;
         
     const int c=gcd(a, M,&h_a, &h_m);
     const int p=a/c;
     const int q=M/c;
     const int d=b/p;
     
     ComplexMatrix gf(p*q*R,c*d);
     
     wfac((const ltfat_complex*)g.data(),L,R,a,M,
	  (ltfat_complex*)gf.fortran_vec());
     
     return octave_value (gf);

  }
  else
  {

     const Matrix g = args(0).matrix_value();
     const int L = g.rows();
     const int R = g.columns();
     const int b = L/M;
         
     const int c=gcd(a, M,&h_a, &h_m);
     const int p=a/c;
     const int q=M/c;
     const int d=b/p;
     
     ComplexMatrix gf(p*q*R,c*d);
     
     wfac_r(g.data(),L,R,a,M,(ltfat_complex*)gf.fortran_vec());
     
     return octave_value (gf);
  }
}
*/