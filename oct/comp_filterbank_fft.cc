#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_filterbank_fft // change to filename
#define OCTFILEHELP "This function calls the C-library\n c=convsub_fft(...);\n Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type 32 or 64 bit signed integer


static inline void fwd_convsub_fft(const Complex *in,const Complex *filt, 
                                   const octave_idx_type inLen, const octave_idx_type sub,
                                   Complex *out)
{
   convsub_fft_d(reinterpret_cast<const double _Complex *>(in),
                 reinterpret_cast<const double _Complex *>(filt),inLen,sub,
                 reinterpret_cast<double _Complex *>(out));
}

static inline void fwd_convsub_fft(const FloatComplex *in,const FloatComplex *filt,
                                   const octave_idx_type inLen, const octave_idx_type sub,
                                   FloatComplex *out)
{
   convsub_fft_s(reinterpret_cast<const float _Complex *>(in),
                 reinterpret_cast<const float _Complex *>(filt),inLen,sub,
                 reinterpret_cast<float _Complex *>(out));
}


template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
   // Input data
   MArray<LTFAT_TYPE> F = ltfatOctArray<LTFAT_TYPE>(args(0));
   // Cell aray containing impulse responses
   Cell G = args(1).cell_value();
   // Subsampling factors
   Matrix a = args(2).matrix_value();

   // Input length
   const octave_idx_type L  = F.rows();
   // Number of channels
   const octave_idx_type W  = F.columns();
   // Number of filters
   const octave_idx_type M = G.nelem();

   // Allocating temporary arrays
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
      g_elems[m] = ltfatOctArray<LTFAT_TYPE>(G.elem(m));
      gPtrs[m] = g_elems[m].data();
      outLen[m] = (octave_idx_type) ceil( L/a(m) );
      c_elems[m] = MArray<LTFAT_TYPE>(dim_vector(outLen[m], W));
      c_elems[m].fill(0);
      cPtrs[m] = c_elems[m].fortran_vec();
   }


   for(octave_idx_type m =0; m<M; m++)
   {
      for(octave_idx_type w =0; w<W; w++)
      {
         // Obtain pointer to w-th column in input
         const LTFAT_TYPE *fPtrCol = F.data() + w*L;
         // Obtaing pointer to w-th column in m-th element of output cell-array
         LTFAT_TYPE *cPtrCol = cPtrs[m] + w*outLen[m];
         fwd_convsub_fft(fPtrCol,gPtrs[m],L,(octave_idx_type)a(m),cPtrCol);
      }
   }

   Cell c(dim_vector(M,1));
   for(octave_idx_type m=0;m<M;++m)
   {
      c.elem(m) = c_elems[m];
   }
   return octave_value(c);
}
