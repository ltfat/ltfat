#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_dgt_long // change to filename
#define OCTFILEHELP "This function calls the C-library\n c=comp_dgt_long(f,g,a,M);\n Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type is 32 or 64 bit signed integer
/*
  dgt_long forwarders
*/

static inline void fwd_dgt_long(const Complex *f, const Complex *g,
                              const octave_idx_type L, const octave_idx_type W,
							  const octave_idx_type a, const octave_idx_type M,
							  Complex *cout)
{
   dgt_long_d(reinterpret_cast<const double (*)[2]>(f),
              reinterpret_cast<const double (*)[2]>(g),
              L,W,a,M,
			  reinterpret_cast<double (*)[2]>(cout));
}

static inline void fwd_dgt_long(const FloatComplex *f, const FloatComplex *g,
                              const octave_idx_type L, const octave_idx_type W, 
							  const octave_idx_type a, const octave_idx_type M,
							  FloatComplex *cout)
{
   dgt_long_s(reinterpret_cast<const float (*)[2]>(f),
              reinterpret_cast<const float (*)[2]>(g),
              L,W,a,M,
			  reinterpret_cast<float (*)[2]>(cout));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 const octave_idx_type a = args(2).int_value();
     const octave_idx_type M = args(3).int_value();
	 
	 MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
	 MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(1));
	 const octave_idx_type L  = f.rows();
     const octave_idx_type W  = f.columns();
     const octave_idx_type N = L/a;
	 
	 dim_vector dims_out(M,N,W);  
     dims_out.chop_trailing_singletons();

     MArray<LTFAT_TYPE> cout(dims_out); 
     cout.fill(0);
	 
	 fwd_dgt_long(f.data(),g.data(),L,W,a,M,cout.fortran_vec());
	 
     return octave_value(cout);
}

/* OLD CODE
#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"


DEFUN_DLD (comp_dgt_long, args, ,
  "This function calls the C-library\n\
  c=comp_dgt_long(f,g,a,M);\n\
  Yeah.")
{

  //const bool f_is_complex  = args(0).is_complex_type();

  const int a = args(2).int_value();
  const int M = args(3).int_value();

  const ComplexMatrix f = args(0).complex_matrix_value();
  const ComplexMatrix g = args(1).complex_matrix_value();
  
  const int L = f.rows();
  const int W = f.columns();
  const int N = L/a;

  dim_vector dims_out(M,N,W);  
  dims_out.chop_trailing_singletons();

  ComplexNDArray cout(dims_out);  
  
  dgt_long((ltfat_complex*)f.data(),(ltfat_complex*)g.data(),
	   L, W, a, M,(ltfat_complex*)cout.data());
  
  return octave_value (cout);



  if (f_is_complex)
  {


  }
  else
  {

    const Matrix f = args(0).matrix_value();
    ltfat_complex *gf;

    const int L = f.rows();
    const int W = f.columns();
    const int N = L/a;

    gf = (ltfat_complex*)ltfat_malloc(L*sizeof(ltfat_complex));
    
    if (args(1).is_complex_type())
    {
       const ComplexMatrix g = args(1).complex_matrix_value();
       wfac((ltfat_complex*)g.data(),L,1,a,M,gf);
    }
    else
    {
       const Matrix g = args(1).matrix_value();
       wfac_r((double*)g.data(),L,1,a,M,gf);
    }
    
    ComplexMatrix cout(M,N*W);  
    
    dgt_fac_r((double*)f.data(),gf,L, W, a, M,
	      (ltfat_complex*)cout.data());
    
    ltfat_free(gf);

    return octave_value (cout);

  }

  

}
*/