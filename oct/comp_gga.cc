#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_gga // change to filename
#define OCTFILEHELP "This function calls the C-library\n c=gga(...);\n Yeah."
//#define _DEBUG


#include "ltfat_oct_template_helper.h"
// octave_idx_type 32 or 64 bit signed integer
 

static inline void fwd_gga(const Complex *fPtr, const double*  indVecPtr,
                           const octave_idx_type L, const octave_idx_type W,
						   const octave_idx_type M, Complex *cPtr )
{
      gga_cd(reinterpret_cast<const double _Complex *>(fPtr),
	         indVecPtr,L,W,M,
			 reinterpret_cast<double _Complex *>(cPtr));
}

static inline void fwd_gga(const FloatComplex *fPtr, const double*  indVecPtr,
                           const octave_idx_type L, const octave_idx_type W,
						   const octave_idx_type M, FloatComplex *cPtr )
{
      gga_cs(reinterpret_cast<const float _Complex *>(fPtr),
	         indVecPtr,L,W,M,
			 reinterpret_cast<float _Complex *>(cPtr));
}

static inline void fwd_gga(const double *fPtr, const double*  indVecPtr,
                           const octave_idx_type L, const octave_idx_type W,
						   const octave_idx_type M, Complex *cPtr )
{
      gga_d(fPtr, indVecPtr,L,W,M,
			reinterpret_cast<double _Complex *>(cPtr));
}

static inline void fwd_gga(const float *fPtr, const double*  indVecPtr,
                           const octave_idx_type L, const octave_idx_type W,
						   const octave_idx_type M, FloatComplex *cPtr )
{
      gga_s(fPtr, indVecPtr,L,W,M,
			reinterpret_cast<float _Complex *>(cPtr));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
	 //DEBUGINFO;
	 // Input data
	 MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
	 MArray<double> indVec = ltfatOctArray<double>(args(1));
	 
	 // Input length
	 const octave_idx_type L  = f.rows();
	 // Number of channels
     const octave_idx_type W  = f.columns();
	 // Number of coefficients
     const octave_idx_type M = indVec.nelem();

     //dims_out.chop_trailing_singletons();
	 MArray<LTFAT_COMPLEX> c(dim_vector(M,W));
	 c.fill(0);
	 
	 fwd_gga(f.data(),indVec.data(),L,W,M,c.fortran_vec());

     return octave_value(c);
}
