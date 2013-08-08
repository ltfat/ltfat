#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_ifilterbank_td // change to filename
#define OCTFILEHELP "This function calls the C-library\n c=upconv_td(...);\n Yeah."
#define _DEBUG


#include "ltfat_oct_template_helper.h"
// octave_idx_type 32 or 64 bit signed integer


static inline void fwd_upconv_td(const Complex *in, const octave_idx_type inLen,
                                  Complex *out, const octave_idx_type outLen,
								  const Complex *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   upconv_td_cd(reinterpret_cast<const double _Complex *>(in),inLen,
                 reinterpret_cast<double _Complex *>(out),outLen,
			     reinterpret_cast<const double _Complex *>(filt),fLen,
                 sub,skip,ext);
}

static inline void fwd_upconv_td(const FloatComplex *in, const octave_idx_type inLen,
                                  FloatComplex *out, const octave_idx_type outLen,
								  const FloatComplex *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   upconv_td_cs(reinterpret_cast<const float _Complex *>(in),inLen,
                 reinterpret_cast<float _Complex *>(out),outLen,
			     reinterpret_cast<const float _Complex *>(filt),fLen,
                 sub,skip,ext);
}

static inline void fwd_upconv_td(const double *in, const octave_idx_type inLen,
                                  double *out, const octave_idx_type outLen,
								  const double *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   upconv_td_d(in,inLen,out,outLen,filt,fLen,sub,skip,ext);
}

static inline void fwd_upconv_td(const float *in, const octave_idx_type inLen,
                                  float *out, const octave_idx_type outLen,
								  const float *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   upconv_td_s(in,inLen,out,outLen,filt,fLen,sub,skip,ext);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
	 //DEBUGINFO;
	 // Input data
	 Cell c = args(0).cell_value();
	 // Cell aray containing impulse responses
	 Cell g = args(1).cell_value();
	 // Subsampling factors
	 Matrix a = args(2).matrix_value();
	 // Skips
	 const octave_idx_type Ls = args(3).int_value();
	 Matrix offset = args(4).matrix_value();

	 charMatrix ext = args(5).char_matrix_value ();
	 // Number of filters
     const octave_idx_type M = g.nelem();

	 // Allocating temporary arrays
	 // Filter lengts
	 OCTAVE_LOCAL_BUFFER (octave_idx_type, filtLen, M);
	 // Output subband lengths
	 OCTAVE_LOCAL_BUFFER (octave_idx_type, Lc, M);
	 // Impulse responses pointers
	 OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, gPtrs, M);
	 // Output subbands pointers
	 OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, cPtrs, M);
	 // Input cell elements array,
	 OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, c_elems, M);
	 //
	 OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, g_elems, M);

	 for(octave_idx_type m=0;m<M;m++)
     {
	    g_elems[m] = ltfatOctArray<LTFAT_TYPE>(g.elem(m));
		c_elems[m] = ltfatOctArray<LTFAT_TYPE>(c.elem(m));
	    gPtrs[m] = g_elems[m].data();
		cPtrs[m] = c_elems[m].data();
        filtLen[m] = g_elems[m].nelem();
		Lc[m] = c_elems[m].rows();
     }

	 const octave_idx_type W  = c_elems[0].columns();

	 MArray<LTFAT_TYPE> f(dim_vector(Ls,W));
	 f.fill(0);
	 LTFAT_TYPE* fPtr = f.fortran_vec();

	 for(octave_idx_type m =0; m<M; m++)
        {
          for(octave_idx_type w =0; w<W; w++)
          {
           // Obtain pointer to w-th column in input
           LTFAT_TYPE *fPtrCol = fPtr + w*Ls;
           // Obtaing pointer to w-th column in m-th element of output cell-array
           const LTFAT_TYPE *cPtrCol = cPtrs[m] + w*Lc[m];
           fwd_upconv_td(cPtrCol,Lc[m],fPtrCol,Ls,gPtrs[m],filtLen[m],a(m),-offset(m),ltfatExtStringToEnum(ext.row_as_string(0).c_str()));
          }
        }

     return octave_value(f);
}
