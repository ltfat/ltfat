#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_iatrousfilterbank_td // change to filename
#define OCTFILEHELP "This function calls the C-library\n c=upconv_td(...);\n Yeah."
#define _DEBUG


#include "ltfat_oct_template_helper.h"
// octave_idx_type 32 or 64 bit signed integer


static inline void fwd_atrousupconv_td(const Complex *in, const octave_idx_type inLen,
                                  Complex *out, const octave_idx_type outLen,
								  const Complex *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   atrousupconv_td_cd(reinterpret_cast<const double _Complex *>(in),inLen,
                 reinterpret_cast<double _Complex *>(out),outLen,
			     reinterpret_cast<const double _Complex *>(filt),fLen,
                 sub,skip,ext);
}

static inline void fwd_atrousupconv_td(const FloatComplex *in, const octave_idx_type inLen,
                                  FloatComplex *out, const octave_idx_type outLen,
								  const FloatComplex *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   atrousupconv_td_cs(reinterpret_cast<const float _Complex *>(in),inLen,
                 reinterpret_cast<float _Complex *>(out),outLen,
			     reinterpret_cast<const float _Complex *>(filt),fLen,
                 sub,skip,ext);
}

static inline void fwd_atrousupconv_td(const double *in, const octave_idx_type inLen,
                                  double *out, const octave_idx_type outLen,
								  const double *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   atrousupconv_td_d(in,inLen,out,outLen,filt,fLen,sub,skip,ext);
}

static inline void fwd_atrousupconv_td(const float *in, const octave_idx_type inLen,
                                  float *out, const octave_idx_type outLen,
								  const float *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   atrousupconv_td_s(in,inLen,out,outLen,filt,fLen,sub,skip,ext);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
	 //DEBUGINFO;
	 // Input data
	 MArray<LTFAT_TYPE> c = ltfatOctArray<LTFAT_TYPE>(args(0));
	 // Cell aray containing impulse responses
	 MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(1));
	 // Subsampling factors
	 Matrix a = args(2).matrix_value();
	 Matrix offset = args(3).matrix_value();

	 const octave_idx_type L = c.dim1();
	 const octave_idx_type M = c.dim2();
	 octave_idx_type W = 1;
	 if(c.ndims()>2)
	 {
	  W = c.dim3();
	 }

     const octave_idx_type filtLen = g.rows();

	 //octave_stdout << L << ", "<< M << ", "<< W << ", "<< filtLen << "\n";

	 OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, gPtrs, M);

	 for(octave_idx_type m=0;m<M;m++)
     {
	    gPtrs[m] = g.data() +  m*filtLen;
     }

	 MArray<LTFAT_TYPE> f(dim_vector(L,W));
	 f.fill(0);

	 for(octave_idx_type m =0; m<M; m++)
        {
          for(octave_idx_type w =0; w<W; w++)
          {
           // Obtain pointer to w-th column in input
           LTFAT_TYPE *fPtrCol = f.fortran_vec() + w*L;
		   const LTFAT_TYPE *cPtrPlane = c.data() + w*L*M;
           // Obtaing pointer to w-th column in m-th element of output cell-array
           const LTFAT_TYPE *cPtrCol = cPtrPlane + m*L;
           fwd_atrousupconv_td(cPtrCol,L,fPtrCol,L,gPtrs[m],filtLen,a(0),-offset(m),ltfatExtStringToEnum("per"));
          }
        }

     return octave_value(f);
}
