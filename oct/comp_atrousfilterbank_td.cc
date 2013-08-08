#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_atrousfilterbank_td // change to filename
#define OCTFILEHELP "This function calls the C-library\n c=convsub_td(...);\n Yeah."
#define _DEBUG


#include "ltfat_oct_template_helper.h"
// octave_idx_type 32 or 64 bit signed integer


static inline void fwd_atrousconvsub_td(const Complex *in, const octave_idx_type inLen,
                                  Complex *out, const octave_idx_type outLen,
								  const Complex *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   atrousconvsub_td_cd(reinterpret_cast<const double _Complex *>(in),inLen,
                 reinterpret_cast<double _Complex *>(out),outLen,
			     reinterpret_cast<const double _Complex *>(filt),fLen,
                 sub,skip,ext);
}

static inline void fwd_atrousconvsub_td(const FloatComplex *in, const octave_idx_type inLen,
                                  FloatComplex *out, const octave_idx_type outLen,
								  const FloatComplex *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   atrousconvsub_td_cs(reinterpret_cast<const float _Complex *>(in),inLen,
                 reinterpret_cast<float _Complex *>(out),outLen,
			     reinterpret_cast<const float _Complex *>(filt),fLen,
                 sub,skip,ext);
}

static inline void fwd_atrousconvsub_td(const double *in, const octave_idx_type inLen,
                                  double *out, const octave_idx_type outLen,
								  const double *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   atrousconvsub_td_d(in,inLen,out,outLen,filt,fLen,sub,skip,ext);
}

static inline void fwd_atrousconvsub_td(const float *in, const octave_idx_type inLen,
                                  float *out, const octave_idx_type outLen,
								  const float *filt, const octave_idx_type fLen,
                                  const octave_idx_type sub, const octave_idx_type skip,
		                          enum ltfatWavExtType ext)
{
   atrousconvsub_td_s(in,inLen,out,outLen,filt,fLen,sub,skip,ext);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
	 //DEBUGINFO;
	 // Input data
	 MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
	 MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(1));
	 Matrix a = args(2).matrix_value();
	 Matrix offset = args(3).matrix_value();

	 // Input length
	 const octave_idx_type L  = f.rows();
	 // Number of channels
     const octave_idx_type W  = f.columns();
	 // Number of filters
     const octave_idx_type M = g.columns();
	 const octave_idx_type filtLen = g.rows();

	 // Allocating temporary arrays
	 OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, gPtrs, M);
	 for(octave_idx_type m=0;m<M;m++)
     {
	    gPtrs[m] = g.data() +  m*filtLen;
     }

     dim_vector dims_out(L,M,W);
     dims_out.chop_trailing_singletons();
	 MArray<LTFAT_TYPE> c(dims_out);
	 c.fill(0);

	 for(octave_idx_type m =0; m<M; m++)
        {
          for(octave_idx_type w =0; w<W; w++)
          {
           // Obtain pointer to w-th column in input
           const LTFAT_TYPE *fPtrCol = f.data() + w*L;
		   LTFAT_TYPE *cPtrPlane = c.fortran_vec() + w*L*M;
           // Obtaing pointer to w-th column in m-th element of output cell-array
           LTFAT_TYPE *cPtrCol = cPtrPlane + m*L;
           //conv_td_sub(fPtrCol,L,&cPtrCol,outLen[m],(const double**)&gPtrs[m],filtLen[m],1,a[m],skip[m],ltfatExtStringToEnum(ext),0);
           fwd_atrousconvsub_td(fPtrCol,L,cPtrCol,L,gPtrs[m],filtLen,a(0),-offset(m),ltfatExtStringToEnum("per"));
          }
        }

     return octave_value(c);
}
