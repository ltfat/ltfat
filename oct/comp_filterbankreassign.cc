#define TYPEDEPARGS 0, 1, 2
#define SINGLEARGS
#define REALARGS
#define OCTFILENAME comp_filterbankreassign // change to filename
#define OCTFILEHELP "Computes reassigned filterbank spectrogram.\n\
                     Usage: sr=comp_filterbankreassign(s,itime,ifreq,a,cfreq);\n\
                     Yeah."


#include "ltfat_oct_template_helper.h"
/*
  gabreassign forwarders
*/

static inline void
fwd_filterbankreassign(const double *s[], const double *tgrad[], const double *fgrad[],
                       const ltfatInt N[], const double a[], const double cfreq[],
                       const ltfatInt M, double *sr[], fbreassOptOut* optout)
{
   filterbankreassign_d(s, tgrad, fgrad, N, a, cfreq, M, sr, REASS_DEFAULT, optout);
}

static inline void
fwd_filterbankreassign(const float *s[], const float *tgrad[], const float *fgrad[],
                       const ltfatInt N[], const double a[], const double cfreq[],
                       const ltfatInt M, float *sr[], fbreassOptOut* optout)
{
   filterbankreassign_s(s, tgrad, fgrad, N, a, cfreq, M, sr, REASS_DEFAULT, optout);
}




template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list
octFunction(const octave_value_list& args, int nargout)
{
   DEBUGINFO;
   Cell sCell      = args(0).cell_value();
   Cell tgradCell  = args(1).cell_value();
   Cell fgradCell  = args(2).cell_value();
   Matrix aMat     = args(3).matrix_value();
   Matrix cfreqMat = args(4).matrix_value();

   const ltfatInt M  = (ltfatInt) sCell.nelem();

   // Output cell
   Cell srCell(dim_vector(M, 1));

   OCTAVE_LOCAL_BUFFER (double, aPtr, M);

   OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, sPtr, M);
   OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, tgradPtr, M);
   OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, fgradPtr, M);
   OCTAVE_LOCAL_BUFFER (LTFAT_TYPE*, srPtr, M);
   OCTAVE_LOCAL_BUFFER (ltfatInt, NPtr, M);
   // The same thing but with
   /*OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, s, M);
   OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, tgrad, M);
   OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, fgrad, M);
   OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, sr, M);
   */

   for (octave_idx_type m = 0; m < M; ++m)
   {
      MArray<LTFAT_TYPE> stmp = ltfatOctArray<LTFAT_TYPE>(sCell.elem(m));
      sPtr[m] = stmp.data();
      NPtr[m] = stmp.rows();
      tgradPtr[m] = (ltfatOctArray<LTFAT_TYPE>(tgradCell.elem(m))).data();
      fgradPtr[m] = (ltfatOctArray<LTFAT_TYPE>(fgradCell.elem(m))).data();
      stmp = MArray<LTFAT_TYPE>(dim_vector(NPtr[m], 1));
      srPtr[m] = stmp.fortran_vec();
      srCell.elem(m) = stmp;
      aPtr[m] = aMat(m);
   }
   octave_value_list retlist;
   retlist(0) = srCell;

   fbreassOptOut* optout = NULL;
   Cell repos;
   octave_idx_type sumN = 0;
   if (nargout > 1)
   {
      for (octave_idx_type ii = 0; ii < M; ii++) sumN += NPtr[ii];

      repos = Cell(dim_vector(sumN, 1));
      optout = fbreassOptOut_init(sumN, 16 );
   }


   // Adjust a
   if (aMat.columns() > 1)
   {
      for (octave_idx_type m = 0; m < M; ++m)
      {
         aPtr[m] /= aMat(m + M);
      }
   }

   fwd_filterbankreassign(sPtr, tgradPtr, fgradPtr, NPtr, aPtr, cfreqMat.data(), M,
                          srPtr, optout);

   if (nargout > 1)
   {
      for (ltfatInt ii = 0; ii < M; ii++)
      {
         octave_idx_type l = optout->reposl[ii];
         MArray<double> cEl = MArray<double>(dim_vector(l, l > 0 ? 1 : 0));
         repos.elem(ii) = cEl;
         double* cElPtr = cEl.fortran_vec();
         for (octave_idx_type jj = 0; jj < optout->reposl[ii]; jj++)
         {
            cElPtr[jj] = (double) (optout->repos[ii][jj] + 1.0);
         }
      }

      retlist(1) = repos;
      fbreassOptOut_destroy(optout);
   }
   return retlist;
}
