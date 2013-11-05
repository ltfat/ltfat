#include "block_interface.h"

void defaultSetter(biEntry* obj,mxArray* in)
{
   if(obj->var!=NULL)
      mxFree(obj->var);
	  
   obj->var=mxDuplicateArray(in);
   mxMakeArrayPersistent(obj->var);
}

mxArray* defaultGetter(biEntry* obj)
{
   return obj->var;
}


static biEntry vars[] = 
{
   // Default setters and getters
   {.name="Ls",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="Pos",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="BufCount",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="PlayChanList",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="RecChanList",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="PageList",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="PageNo",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="Skipped",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="BufLen",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="ClassId",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="AnaOverlap",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="SynOverlap",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="DispLoad",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   // Default setters and custom getters
   {.name="Source",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   {.name="",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
}


void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
   if(nlhs<1)
   {
     mexErrMsgText("BLOCK_INTERFACE: Not enough input arguments.");
   }

   mxArray* mxCmd = prhs[0];
   char command[COMMAND_LENGTH];
   size_t comStrLen = mxGetNumberOfElement(mxCmd);
   mxGetString(mxCmd,command,comStrLen+1);

   if(!strcmp(command,"set"))
   {
   }
   else if(!strcmp(command,"get"))
   {
   } 
 
}
