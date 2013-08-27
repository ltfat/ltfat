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
#include "mat.h"
#include "matLoader.h"




/*
There is probably not a multiplatform way of hand loading of shared libraries, therefore 2 versions.
*/
#if defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#include <tchar.h>

typedef void (*mexFunction_t)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

int main(int argc,char* argv[])
{

   if(argc<2)
   {
      fprintf(stderr,"Correct parameters: NAMEOFMEX.MEX ARGS.mat\n");
      return -1;
   }

// NOTE: When specifying a path, be sure to use backslashes (\), not forward slashes (/).
   HINSTANCE dllHandle = NULL;
   dllHandle = LoadLibrary((argv[1]));
   if(dllHandle==NULL){
      fprintf(stderr, "Error loading MEX file.\n");
      return -1;
   }

   mexFunction_t mexfunction = NULL;
   mexfunction = (mexFunction_t)GetProcAddress(dllHandle,"mexFunction");
   if(mexfunction==NULL){
      fprintf(stderr, "MEX file does not contain mexFunction\n");
      return -1;
   }

   int argCount = getNoOfArgsFromMAT(argv[2]);
   const mxArray* prhs[argCount];
   getArgsFromMAT(argv[2],(mxArray**)prhs,argCount);
   mxArray* plhs[1];

   mexfunction(1,plhs,1,prhs);

   for(int ii=0;ii<argCount;ii++)
      mxDestroyArray((mxArray*)prhs[ii]);
   //Free the library:
   BOOL freeResult = FreeLibrary(dllHandle);

}
#else

#endif
