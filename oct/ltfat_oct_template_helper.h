#if defined(OCTFILENAME) && defined(OCTFILEHELP)
#ifndef _LTFAT_OCT_TEMPLATE_HELPER_H
#define _LTFAT_OCT_TEMPLATE_HELPER_H
#include <octave/oct.h>

#include "ltfat.h"

#ifdef _DEBUG
#define DEBUGINFO  octave_stdout << __PRETTY_FUNCTION__ << "\n"
#else
#define DEBUGINFO
#endif

bool checkIsSingle(const octave_value& ov);
octave_value recastToSingle(const octave_value& ov);

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout);

template <class LTFAT_TYPE>
MArray<LTFAT_TYPE> ltfatOctArray(const octave_value& ov);

template <class LTFAT_TYPE>
MArray<LTFAT_TYPE> ltfatOctArray(const octave_value& ov)
{
   error("Casting to unknown type. Everything should be handled in the specialized functions.",__PRETTY_FUNCTION__);
   return MArray<LTFAT_TYPE>();
}

#define LTFAT_ARRAY_SPECIALIZED_REAL(real_type)               \
template <>                                                   \
MArray<real_type> ltfatOctArray(const octave_value& ov)           \
{                                                             \
    if(ov.is_complex_type())                                  \
    {                                                         \
        error("Casting from complex to real is forbidden (for now).");  \
    }                                                         \
    if(ov.is_single_type())                                   \
    {                                                         \
         return (ov.float_array_value());                     \
    }                                                         \
    else if(ov.is_double_type())                              \
    {                                                         \
         return (ov.array_value());                           \
    }                                                         \
    else                                                      \
    {                                                         \
       error("Unsupported data type..");                      \
    }                                                         \
    return MArray<real_type>();                                \
}

#define LTFAT_ARRAY_SPECIALIZED_COMPLEX(complex_type)               \
template <>                                                         \
MArray<complex_type> ltfatOctArray(const octave_value& ov)          \
{                                                                   \
    if(ov.is_single_type())                                         \
    {                                                               \
       if(ov.is_complex_type())                                     \
         return (ov.float_complex_array_value());                   \
       else                                                         \
         return (ov.float_array_value());                           \
    }                                                               \
    else if(ov.is_double_type())                                    \
    {                                                               \
       if(ov.is_complex_type())                                     \
          return (ov.complex_array_value());                        \
       else                                                         \
          return (ov.array_value());                                \
    }                                                               \
    else                                                            \
    {                                                               \
       error("Unsupported data type..");                            \
    }                                                               \
    return MArray<complex_type>();                                   \
}

LTFAT_ARRAY_SPECIALIZED_REAL(double)
LTFAT_ARRAY_SPECIALIZED_REAL(float)
LTFAT_ARRAY_SPECIALIZED_COMPLEX(Complex)
LTFAT_ARRAY_SPECIALIZED_COMPLEX(FloatComplex)

#undef LTFAT_ARRAY_SPECIALIZED_REAL
#undef LTFAT_ARRAY_SPECIALIZED_COMPLEX

bool checkIsSingle(const octave_value& ov)
{
   if(ov.is_cell())
   {
      Cell ov_cell = ov.cell_value();
      for(int jj=0;jj<ov_cell.numel();jj++)
      {
         if(!checkIsSingle(ov_cell.elem(jj)))
            return false;
      }
      return true;
   }
   return ov.is_single_type();
}

octave_value recastToSingle(const octave_value& ov)
{
   if(checkIsSingle(ov))
   {
      return ov;
   }

   if(ov.is_cell())
   {
      Cell ov_cell = ov.cell_value();
      Cell ovtmp_cell(ov.dims());
      for(int jj=0;jj<ovtmp_cell.nelem();jj++)
      {
         ovtmp_cell(jj) = recastToSingle(ov_cell.elem(jj));
      }
      return ovtmp_cell;
   }
   /*
   TODO: ov is struct
   */
   // just copy pointer if the element is not numeric
   if(!ov.is_numeric_type())
   {
      return ov;
   }

   if(ov.is_complex_type())
   {
      return ltfatOctArray<FloatComplex>(ov);
   }
   else
   {
      return ltfatOctArray<float>(ov);
   }
}


DEFUN_DLD (OCTFILENAME, args, nargout, OCTFILEHELP)
{
octave_value_list argsCopy(args);

bool isAnySingle = false;
bool isAnyComplex = false;
#ifndef TYPEDEPARGS
return octFunction<double,double,Complex>(argsCopy,nargout);
#else
int prhsToCheck[] = { TYPEDEPARGS };
int phrsToCheckLen = sizeof(prhsToCheck)/sizeof(*prhsToCheck);

octave_value_list typedepArgs;
for(int ii=0;ii<phrsToCheckLen;ii++)
{
   typedepArgs.append(argsCopy(prhsToCheck[ii]));
}

// WORKAROUND Incorrect detection of the single data type of complex diag. matrices
for(int ii=0;ii<typedepArgs.length();ii++)
    if(typedepArgs(ii).is_diag_matrix())
       typedepArgs(ii)= typedepArgs(ii).full_value();


for(int ii=0;ii<typedepArgs.length();ii++)
    if((isAnySingle=checkIsSingle(typedepArgs(ii)))) break;

for(int ii=0;ii<typedepArgs.length();ii++)
    if((isAnyComplex=typedepArgs(ii).is_complex_type())) break;

// just if the reference was changed
for(int ii=0;ii<phrsToCheckLen;ii++)
   argsCopy(prhsToCheck[ii])=typedepArgs(ii);

#if defined(REALARGS)&& !(defined(COMPLEXARGS) || defined(COMPLEXINDEPENDENT)) 
if(isAnyComplex)
{
   error("Only real inputs are accepted.");
   return octave_value_list();
}
#endif
 

 
#if defined(COMPLEXINDEPENDENT)
if(isAnySingle)
{   
   if(isAnyComplex)
   {
      for(int ii=0;ii<typedepArgs.length();ii++)
         typedepArgs(ii) = octave_value(recastToSingle(argsCopy(ii)));
		 
      return octFunction<FloatComplex,float,FloatComplex>(argsCopy,nargout);
   }
   else
   {
      for(int ii=0;ii<typedepArgs.length();ii++)
         typedepArgs(ii) = octave_value(recastToSingle(argsCopy(ii)));
		 
      return octFunction<float,float,FloatComplex>(argsCopy,nargout); 
   } 
}
else
{
   if(isAnyComplex)
   {
      for(int ii=0;ii<typedepArgs.length();ii++)
         typedepArgs(ii) = octave_value(ltfatOctArray<Complex>(argsCopy(ii)));
		 
      return octFunction<Complex,double,Complex>(argsCopy,nargout);
   }
   else
   {
      return octFunction<double,double,Complex>(argsCopy,nargout);      
   }
}
#elif defined(COMPLEXARGS) && !defined(REALARGS)
if(isAnySingle)
{   
   for(int ii=0;ii<typedepArgs.length();ii++)
      typedepArgs(ii) = octave_value(recastToSingle(argsCopy(ii)));
	 
   return octFunction<FloatComplex,float,FloatComplex>(argsCopy,nargout);

}
else
{
   for(int ii=0;ii<typedepArgs.length();ii++)
      typedepArgs(ii) = octave_value(ltfatOctArray<Complex>(argsCopy(ii)));
		 
   return octFunction<Complex,double,Complex>(argsCopy,nargout);
}
#elif defined(REALARGS)
if(isAnySingle)
{   
   #if defined(COMPLEXARGS)
    for(int ii=0;ii<typedepArgs.length();ii++)
      typedepArgs(ii) = octave_value(recastToSingle(argsCopy(ii)));
   #endif

   return octFunction<float,float,FloatComplex>(argsCopy,nargout);
}
else
{
   #if defined(COMPLEXARGS)
   for(int ii=0;ii<typedepArgs.length();ii++)
      typedepArgs(ii) = octave_value(ltfatOctArray<Complex>(argsCopy(ii)));
   #endif
   
   return octFunction<double,double,Complex>(argsCopy,nargout);
}
#else
error("Something wrong in the template system. My bad....\n");
#endif 



#endif // TYPEDEPARGS


error("Something fishy is going on in...\n");

return octave_value_list();


}


#endif // _LTFAT_OCT_TEMPLATE_HELPER_H
#endif // defined(OCTFILENAME) && defined(OCTFILEHELP)
