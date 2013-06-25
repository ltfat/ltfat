#define TYPEDEPARGS 0
#define REALARGS
#define COMPLEXARGS
#define SINGLEARGS
#define OCTFILENAME comp_dct // change to filename
#define OCTFILEHELP "This function calls FFTW library.\n Yeah."

#include "ltfat_oct_template_helper.h"
#include "config.h"
// octave_idx_type is 32 or 64 bit signed integer

static inline void fwd_dct(double *f,
                           const octave_idx_type L,
						   const octave_idx_type W,
						   const fftw_r2r_kind *kind)
{
  fftw_iodim dims[1], howmanydims[1];
  fftw_plan p;

  dims[0].n = L;
  dims[0].is = 1;
  dims[0].os = 1;

  howmanydims[0].n = W;
  howmanydims[0].is = L;
  howmanydims[0].os = L;

  p = fftw_plan_guru_r2r(1, dims,
				   1, howmanydims,
				   f, f,
				   kind,
				   FFTW_OPTITYPE);

  // Real FFT.
  fftw_execute(p);

  fftw_destroy_plan(p);
}

static inline void fwd_dct(float *f,
                           const octave_idx_type L,
						   const octave_idx_type W,
						   const fftw_r2r_kind *kind)
{
  fftwf_iodim dims[1], howmanydims[1];
  fftwf_plan p;

  dims[0].n = L;
  dims[0].is = 1;
  dims[0].os = 1;

  howmanydims[0].n = W;
  howmanydims[0].is = L;
  howmanydims[0].os = L;

  p = fftwf_plan_guru_r2r(1, dims,
				   1, howmanydims,
				   f, f,
				   kind,
				   FFTW_OPTITYPE);

  // Real FFT.
  fftwf_execute(p);

  fftwf_destroy_plan(p);
}

static inline void fwd_dct(Complex *f,
                           const octave_idx_type L,
						   const octave_idx_type W,
						   const fftw_r2r_kind *kind)
{
  fftw_iodim dims[1], howmanydims[1];
  fftw_plan p;
  double *f_ptr = reinterpret_cast<double*>(f);

  dims[0].n = L;
  dims[0].is = 2;
  dims[0].os = 2;

  howmanydims[0].n = W;
  howmanydims[0].is = 2*L;
  howmanydims[0].os = 2*L;

  p = fftw_plan_guru_r2r(1, dims,
				   1, howmanydims,
				   f_ptr, f_ptr,
				   kind,
				   FFTW_OPTITYPE|FFTW_UNALIGNED);

  fftw_execute(p);
  fftw_execute_r2r(p,f_ptr+1,f_ptr+1);

  fftw_destroy_plan(p);



}

static inline void fwd_dct(FloatComplex *f,
                           const octave_idx_type L,
						   const octave_idx_type W,
						   const fftw_r2r_kind *kind)
{
  fftwf_iodim dims[1], howmanydims[1];
  fftwf_plan p;
  float *f_ptr = reinterpret_cast<float*>(f);

  dims[0].n = L;
  dims[0].is = 2;
  dims[0].os = 2;

  howmanydims[0].n = W;
  howmanydims[0].is = 2*L;
  howmanydims[0].os = 2*L;

  p = fftwf_plan_guru_r2r(1, dims,
				   1, howmanydims,
				   f_ptr, f_ptr,
				   kind,
				   FFTW_OPTITYPE|FFTW_UNALIGNED);

  // Real FFT.
  fftwf_execute(p);
  fftwf_execute_r2r(p,f_ptr+1,f_ptr+1);


  fftwf_destroy_plan(p);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 fftw_r2r_kind kind[1];
	 const octave_idx_type type = args(1).int_value();
	 MArray<LTFAT_TYPE> f = MArray<LTFAT_TYPE>(ltfatOctArray<LTFAT_TYPE>(args(0)));
	 const octave_idx_type L  = f.rows();
     const octave_idx_type W  = f.columns();
     octave_idx_type N = 2*L;
	 LTFAT_REAL sqrt2 = (LTFAT_REAL) sqrt(2.0);
	 LTFAT_REAL postScale = (LTFAT_REAL) 1.0/sqrt2;
     LTFAT_REAL scale = (LTFAT_REAL) sqrt2*(1.0/(double)N)*sqrt((double)L);

	 LTFAT_TYPE* f_ptr = f.fortran_vec();

	 // Pre-scaling
	 if(type==1||type==3)
     {
       for(octave_idx_type ii=0;ii<W;ii++)
       {
	      f_ptr[ii*L] *= sqrt2;
       }
	 }

     switch(type)
     {
	    case 1:
		   N -= 2;
		   for(octave_idx_type ii=0;ii<W;ii++)
           {
	          f_ptr[(ii+1)*L-1] *= sqrt2;
           }

		   scale = (LTFAT_REAL) sqrt2*(1.0/((double)N))*sqrt((double)L-1);
           kind[0] = FFTW_REDFT00;
        break;
	    case 2:
           kind[0] = FFTW_REDFT10;
        break;
	    case 3:
           kind[0] = FFTW_REDFT01;
        break;
	    case 4:
           kind[0] = FFTW_REDFT11;
        break;
        default:
		   error("Unknown type.");
     }

	 fwd_dct(f_ptr,L,W,kind);

	 // Post-scaling
	 for(int ii=0;ii<L*W;ii++)
     {
	    f_ptr[ii] *= scale;
     }

     if(type==1||type==2)
     {
        // Scale DC component
        for(int ii=0;ii<W;ii++)
        {
	       f_ptr[ii*L] *= postScale;
        }
     }

     if(type==1)
     {
        // Scale AC component
        for(int ii=0;ii<W;ii++)
        {
	       f_ptr[(ii+1)*L-1] *= postScale;
        }
     }

     return octave_value(f);
}


