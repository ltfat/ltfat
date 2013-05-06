#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_gabdual_long // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                    gd=comp_gabdual_long(g,a,M);\n\
                    Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type is 32 or 64 bit signed integer
/*
  gabdual_long forwarders
*/

static inline void fwd_gabdual_long(const Complex *g,
                                  const octave_idx_type L, 
							      const octave_idx_type R, 
								  const octave_idx_type a, 
							      const octave_idx_type M, 
								  Complex *gd)
{
   gabdual_long_d(reinterpret_cast<const double (*)[2]>(g),
				L,R,a,M,
				reinterpret_cast<double (*)[2]>(gd));
}

static inline void fwd_gabdual_long(const FloatComplex *g, 
                                  const octave_idx_type L, 
							      const octave_idx_type R, 
								  const octave_idx_type a, 
							      const octave_idx_type M, 
								  FloatComplex *gd)
{
   gabdual_long_s(reinterpret_cast<const float (*)[2]>(g),
				L,R,a,M,
				reinterpret_cast<float (*)[2]>(gd));
}

static inline void fwd_gabdual_long(const double *g,
                                  const octave_idx_type L, 
							      const octave_idx_type R, 
								  const octave_idx_type a, 
							      const octave_idx_type M, 
								  double *gd)
{
   gabdualreal_long_d(g,L,R,a,M,gd);
}

static inline void fwd_gabdual_long(const float *g,
                                  const octave_idx_type L, 
							      const octave_idx_type R, 
								  const octave_idx_type a, 
							      const octave_idx_type M, 
								  float *gd)
{
   gabdualreal_long_s(g,L,R,a,M,gd);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(0));
	 const octave_idx_type L = g.rows();
     const octave_idx_type R = g.columns();
     const octave_idx_type a = args(1).int_value();
     const octave_idx_type M = args(2).int_value();
  
     MArray<LTFAT_TYPE> gd(dim_vector(L,R)); 
     gd.fill(0);
	 
	 fwd_gabdual_long(g.data(),L, R, a, M,gd.fortran_vec());
	 
     return octave_value(gd);
}

/*
#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_gabdual_long, args, ,
  "This function calls the C-library\n\
  gd=comp_gabdual_long(g,a,M);\n\
  Yeah.")
{

  if (args(0).is_complex_type())
  {

     const ComplexMatrix g = args(0).complex_matrix_value();
     const int L = g.rows();
     const int R = g.columns();
     const int a = args(1).int_value();
     const int M = args(2).int_value();

     ComplexMatrix gd(L,R);
     
     gabdual_long((const ltfat_complex*)g.fortran_vec(),
		  L,R,a,M,
		  (ltfat_complex*)gd.data());

     return octave_value (gd);

  }
  else
  {

     const Matrix g = args(0).matrix_value();
     const int L = g.rows();
     const int R = g.columns();
     const int a = args(1).int_value();
     const int M = args(2).int_value();

     Matrix gd(L,R);
  
     gabdualreal_long((double*)g.fortran_vec(),
		      L,R,a,M,
		      (double*)gd.fortran_vec());

     return octave_value (gd);

  }
}
*/