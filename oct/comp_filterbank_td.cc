#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_filterbank_td // change to filename
#define OCTFILEHELP "This function calls the C-library\n c=convsub_td(...);\n Yeah."
#define _DEBUG


#include "ltfat_oct_template_helper.h"
// octave_idx_type 32 or 64 bit signed integer


static inline void fwd_convsub_td(const Complex *in, const octave_idx_type inLen,
                                  Complex *out, const octave_idx_type outLen,
								  const Complex *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   convsub_td_cd(reinterpret_cast<const double _Complex *>(in),inLen,
                 reinterpret_cast<double _Complex *>(out),outLen,
			     reinterpret_cast<const double _Complex *>(filt),fLen,
                 sub,skip,ext);
}

static inline void fwd_convsub_td(const FloatComplex *in, const octave_idx_type inLen,
                                  FloatComplex *out, const octave_idx_type outLen,
								  const FloatComplex *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   convsub_td_cs(reinterpret_cast<const float _Complex *>(in),inLen,
                 reinterpret_cast<float _Complex *>(out),outLen,
			     reinterpret_cast<const float _Complex *>(filt),fLen,
                 sub,skip,ext);
}

static inline void fwd_convsub_td(const double *in, const octave_idx_type inLen,
                                  double *out, const octave_idx_type outLen,
								  const double *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   convsub_td_d(in,inLen,out,outLen,filt,fLen,sub,skip,ext);
}

static inline void fwd_convsub_td(const float *in, const octave_idx_type inLen,
                                  float *out, const octave_idx_type outLen,
								  const float *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   convsub_td_s(in,inLen,out,outLen,filt,fLen,sub,skip,ext);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
	 // Input data
	 MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
	 // Cell aray containing impulse responses
	 Cell g = args(1).cell_value();
	 // Subsampling factors
	 Matrix a = args(2).matrix_value();
	 // Skips
	 Matrix offset = args(3).matrix_value();
	 charMatrix ext = args(4).char_matrix_value ();
	 // Input length
	 const octave_idx_type L  = f.rows();
	 // Number of channels
     const octave_idx_type W  = f.columns();
	 // Number of filters
     const octave_idx_type M = g.nelem();

	 // Allocating temporary arrays
	 // Filter lengts
	 OCTAVE_LOCAL_BUFFER (octave_idx_type, filtLen, M);
	 // Output subband lengths
	 OCTAVE_LOCAL_BUFFER (octave_idx_type, outLen, M);
	 // Impulse responses pointers
	 OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, gPtrs, M);
	 // Output subbands pointers
	 OCTAVE_LOCAL_BUFFER (LTFAT_TYPE*, cPtrs, M);
	 // Output cell elements array,
	 OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, c_elems, M);
	 //
	 OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, g_elems, M);

	 for(octave_idx_type m=0;m<M;m++)
     {
	    g_elems[m] = ltfatOctArray<LTFAT_TYPE>(g.elem(m));
	    gPtrs[m] = g_elems[m].data();
        filtLen[m] = g_elems[m].numel();
     }

	 if(ext=="per")
     {
        for(octave_idx_type m = 0; m < M; m++)
        {
           outLen[m] = (octave_idx_type) ceil( L/a(m) );
        }
     }
	 else if(ext=="valid")
	 {
	    for(octave_idx_type m = 0; m < M; m++)
        {
           outLen[m] = (octave_idx_type) ceil( (L-(filtLen[m]-1))/a(m) );
        }
	 }
     else
     {
        for(octave_idx_type m = 0; m < M; m++)
        {
           outLen[m] = (octave_idx_type) ceil( (L + filtLen[m] - 1 + offset(m) )/a(m) );
        }
     }


	 for(octave_idx_type m=0;m<M;++m)
     {
	    c_elems[m] = MArray<LTFAT_TYPE>(dim_vector(outLen[m], W));
		c_elems[m].fill(0);
		cPtrs[m] = c_elems[m].fortran_vec();
     }

	 for(octave_idx_type m =0; m<M; m++)
        {
          for(octave_idx_type w =0; w<W; w++)
          {
           // Obtain pointer to w-th column in input
           const LTFAT_TYPE *fPtrCol = f.data() + w*L;
           // Obtaing pointer to w-th column in m-th element of output cell-array
           LTFAT_TYPE *cPtrCol = cPtrs[m] + w*outLen[m];
           //conv_td_sub(fPtrCol,L,&cPtrCol,outLen[m],(const double**)&gPtrs[m],filtLen[m],1,a[m],skip[m],ltfatExtStringToEnum(ext),0);
           fwd_convsub_td(fPtrCol,L,cPtrCol,outLen[m],gPtrs[m],filtLen[m],a(m),-offset(m),ltfatExtStringToEnum(ext.row_as_string(0).c_str()));
          }
        }

	 Cell c(dim_vector(M,1));
	 for(octave_idx_type m=0;m<M;++m)
     {
	    c.elem(m) = c_elems[m];
	 }
     return octave_value(c);
}
