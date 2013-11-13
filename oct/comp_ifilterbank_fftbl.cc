#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_ifilterbank_fftbl // change to filename
#define OCTFILEHELP "This function calls the C-library\n c=upconv_fftbl(...);\n Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type 32 or 64 bit signed integer


static inline void fwd_upconv_fftbl(const Complex *in, const octave_idx_type inLen, 
                                  const Complex *filt, const octave_idx_type filtLen,
                                  const ptrdiff_t foff, const double afrac, const double realonly,
				  Complex *out)
{
   upconv_fftbl_d(reinterpret_cast<const double _Complex *>(in),inLen,
                  reinterpret_cast<const double _Complex *>(filt),filtLen,foff,afrac,realonly,
                  reinterpret_cast<double _Complex *>(out));
}

static inline void fwd_upconv_fftbl(const FloatComplex *in,const octave_idx_type inLen,
                                  const FloatComplex *filt, const octave_idx_type filtLen,
                                  const ptrdiff_t foff, const double afrac, const double realonly,
                                  FloatComplex *out)
{
   upconv_fftbl_s(reinterpret_cast<const float _Complex *>(in),inLen,
                  reinterpret_cast<const float _Complex *>(filt),filtLen,foff,afrac,realonly,
                  reinterpret_cast<float _Complex *>(out));
}


// Calling convention:
// F = comp_ifilterbank_fftbl(c,G,foff,a,realonly)
template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
   // Input data
   Cell c = args(0).cell_value();
   // Cell aray containing impulse responses
   Cell G = args(1).cell_value();
   // Subsampling factors
   const double* foff = args(2).matrix_value().data();
   const double* a = args(3).matrix_value().data();
   const double* realonly = args(4).matrix_value().data();

   octave_idx_type acols = args(3).matrix_value().columns();
   // Number of channels
   const octave_idx_type W  = c.elem(0).columns();
   // Number of filters
   const octave_idx_type M = G.nelem();

   OCTAVE_LOCAL_BUFFER (double, afrac, M);
   memcpy(afrac,a,M*sizeof(double));
   if(acols>1)
   {
      for(octave_idx_type m=0;m<M;m++)
      {
         afrac[m]/=a[M+m];
      }
   }

   // Allocating 
   // Output subband lengths
   OCTAVE_LOCAL_BUFFER (octave_idx_type, inLen, M);
   OCTAVE_LOCAL_BUFFER (octave_idx_type, filtLen, M);
   // Impulse responses pointers
   OCTAVE_LOCAL_BUFFER (const LTFAT_COMPLEX*, GPtrs, M);
   // Output subbands pointers
   OCTAVE_LOCAL_BUFFER (const LTFAT_COMPLEX*, cPtrs, M);
   // Output cell elements array,
   OCTAVE_LOCAL_BUFFER (MArray<LTFAT_COMPLEX>, c_elems, M);
   //
   OCTAVE_LOCAL_BUFFER (MArray<LTFAT_COMPLEX>, G_elems, M);

   for(octave_idx_type m=0;m<M;m++)
   {
      G_elems[m] = ltfatOctArray<LTFAT_COMPLEX>(G.elem(m));
      GPtrs[m] = G_elems[m].data();
      inLen[m] = c.elem(m).rows();
      c_elems[m] = ltfatOctArray<LTFAT_COMPLEX>(c.elem(m));
      cPtrs[m] = c_elems[m].data();
      filtLen[m] = G_elems[m].nelem();
   }

   // Output length
   const octave_idx_type L  = floor(afrac[0]*inLen[0]+0.5);
   // Output signal
   MArray<LTFAT_COMPLEX> F(dim_vector(L,W));
   F.fill(0);
  
   for(octave_idx_type m =0; m<M; m++)
   {
      for(octave_idx_type w =0; w<W; w++)
      {
         // Obtain pointer to w-th column in input
         LTFAT_COMPLEX *FPtrCol = F.fortran_vec() + w*L;
         // Obtaing pointer to w-th column in m-th element of output cell-array
         const LTFAT_COMPLEX *cPtrCol = cPtrs[m] + w*inLen[m];
         fwd_upconv_fftbl(cPtrCol,inLen[m],GPtrs[m],filtLen[m],(ptrdiff_t)foff[m],afrac[m],
                          realonly[m],FPtrCol);
      }
   }

   return octave_value(F);
}
