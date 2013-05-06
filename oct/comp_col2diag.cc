#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_col2diag // change to filename
#define OCTFILEHELP "Computes spreading permutation.\n Usage: cout=comp_col2diag(cin);\n Yeah."

#include "ltfat_oct_template_helper.h"

/*
  col2diag forwarders
*/

static inline void fwd_col2diag(const Complex *cin, const octave_idx_type L,Complex *cout)
{
   col2diag_cd(reinterpret_cast<const double _Complex*>(cin),L,reinterpret_cast<double _Complex*>(cout));
}

static inline void fwd_col2diag(const FloatComplex *cin, const octave_idx_type L, FloatComplex *cout)
{
   col2diag_cs(reinterpret_cast<const float _Complex*>(cin),L,reinterpret_cast<float _Complex*>(cout));
}

static inline void fwd_col2diag(const double *cin, const octave_idx_type L,double *cout)
{
   col2diag_d(cin,L,cout);
}

static inline void fwd_col2diag(const float *cin, const octave_idx_type L, float *cout)
{
   col2diag_s(cin,L,cout);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
     MArray<LTFAT_TYPE> cin = ltfatOctArray<LTFAT_TYPE>(args(0));
     MArray<LTFAT_TYPE> cout(cin.dims());
     cout.fill(0);

     const octave_idx_type L = cin.rows();

     fwd_col2diag(cin.data(),L,cout.fortran_vec());

     return octave_value(cout);
}

/* // OLD CODE
#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_col2diag, args, ,
  )
{

  const bool cin_is_complex  = args(0).is_complex_type();

  if (cin_is_complex)
  {
     const ComplexMatrix cin = args(0).complex_matrix_value();
     const int L = cin.rows();

     ComplexMatrix cout(L,L);

     col2diag((ltfat_complex*)cin.data(),L,(ltfat_complex*)cout.data());

     return octave_value (cout);

  }
  else
  {

     const Matrix cin = args(0).matrix_value();
     const int L = cin.rows();

     Matrix cout(L,L);

     col2diag_r((double*)cin.data(),L,(double*)cout.data());

     return octave_value (cout);



  }
}
*/


