/*
This loads and executes MEX file from outside od the Matlab.
Based on
http://stackoverflow.com/questions/11220250/how-do-i-profile-a-mex-function-in-matlab/12405131#12405131
and
http://msdn.microsoft.com/en-us/library/ms810279.aspx
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "mat.h"
#include "matLoader.h"
#include "fftw3.h"

/*
Requires Matlab libraries to be installed -lmat and -lmx.

BEWARE! Any function from MEX API working with Matlab runtime is some way will crash this executer.
This is the case for functions with the mex prefix i.e. mexAtExit, mexPrintf etc.

Use saveargsfor.m to store data from Matlab in a format which can be read by this executer.

Calling convention:

mexExecuter m.mexw64 mat.mat
*/


/*
There is probably not a multiplatform way of hand loading of shared libraries, therefore 2 versions.
*/
#if defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

typedef void (*mexFunction_t)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

double compMSE_d(const mxArray* a,const mxArray* b);




int main(int argc,char* argv[])
{

   if(argc<2)
   {
      fprintf(stderr,"Correct parameters: NAMEOFMEX.MEX ARGS.mat\n");
      return -1;
   }

   void* tmp = fftw_malloc(1);
   fftw_free(tmp);
   tmp = fftwf_malloc(1);
   fftwf_free(tmp);

   mexFunction_t mexfunction = NULL;
#if defined(_WIN32) || defined(__WIN32__)
// NOTE: When specifying a path, be sure to use backslashes (\), not forward slashes (/).
   HINSTANCE dllHandle = NULL;
   dllHandle = LoadLibrary((argv[1]));
   if(dllHandle==NULL){
      fprintf(stderr, "Error loading MEX file.\n");
      return -1;
   }

 /*
 Preferably use the mexFunctionInner, which skips registering mexAtExit which causes crashes.
 */
   mexfunction = (mexFunction_t)GetProcAddress(dllHandle,"mexFunctionInner");
   if(mexfunction==NULL){
      mexfunction = (mexFunction_t)GetProcAddress(dllHandle,"mexFunction");
      if(mexfunction==NULL){
         fprintf(stderr, "MEX file does not contain mexFunction\n");
         return -1;
      }
   }
#else
  void *handle = dlopen(argv[1], RTLD_NOW);
  if(!handle){
    fprintf(stderr, "Error loading MEX file: %s\n", dlerror());
    return -1;
  }

  mexfunction = (mexFunction_t)dlsym(handle, "mexFunctionInner");
  if(!mexfunction){
    mexfunction = (mexFunction_t)dlsym(handle, "mexFunction");
    if(!mexfunction){
       fprintf(stderr, "MEX file does not contain mexFunction\n");
       return -1;
    }
  }
#endif

   mxArray* res[1];
   res[0] = NULL;
   int argCount = getNoOfArgsFromMAT(argv[2]);
   const mxArray* prhs[argCount];
   getArgsFromMAT(argv[2],(mxArray**)prhs,(mxArray**)res, argCount);
   mxArray* plhs[1];


   if(res[0]!=NULL)
   {
       argCount--;
   }

   mexfunction(1,plhs,argCount,prhs);

   if(res!=NULL)
   {
     if(mxGetClassID(plhs[0])!=mxDOUBLE_CLASS)
     {
       fprintf(stderr, "Currently working with doubles only.\n");
       return -1;
     }
     double err=compMSE_d(res[0],plhs[0]);
     fprintf(stdout, "Error: %e\n", err);
     mxDestroyArray(res[0]);
   }


   for(int ii=0;ii<argCount;ii++)
      mxDestroyArray((mxArray*)prhs[ii]);

   mxDestroyArray(plhs[0]);

   //Free the library:
#if defined(_WIN32) || defined(__WIN32__)
   //  BOOL freeResult = FreeLibrary(dllHandle);
#else
   dlclose(handle);
#endif

}

double compMSE_d(const mxArray* a,const mxArray* b)
{
   if(a==NULL || b == NULL)
   {
       fprintf(stderr, "NULL pointer\n");
       return -1;
   }

   if(mxIsNumeric(a) && mxIsNumeric(b))
   {
     double* fr = mxGetPr(a);
     double* fhatr = mxGetPr(b);
     mwSize L = mxGetM(a)*mxGetN(a);
     double err = 0;
     for(mwIndex n=0;n<L;n++)
     {
        err += pow(fr[n]-fhatr[n],2);
     }
     err = sqrt(err);


     if(mxIsComplex(a) && mxIsComplex(b))
     {
        double* fi = mxGetPi(a);
        double* fhati = mxGetPi(b);
        for(mwIndex n=0;n<L;n++)
        {
           err += pow(fi[n]-fhati[n],2);
        }
     }
     return err;
   }
   else
   {
       fprintf(stderr, "Arrays are not numeric.\n");
       return -1;
   }


}


