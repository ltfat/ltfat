#include "block_interface.h"
#include <string.h>


void defaultSetter(biEntry* obj,mxArray* in)
{
   if(obj->var!=NULL)
      mxDestroyArray(obj->var);
	  
   obj->var=mxDuplicateArray(in);
   mexMakeArrayPersistent(obj->var);
}

mxArray* defaultGetter(biEntry* obj)
{
   if(obj->var!=NULL)  	
	return obj->var;
   else
	return mxCreateDoubleMatrix(0,0,mxREAL); 
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
   {.name="IsLoop",.var=NULL,.setter=defaultSetter,.getter=defaultGetter},
   // Default setters and custom getters
   {.name="Source",.var=NULL,.setter=defaultSetter,.getter=getSource},
   {.name="EnqBufCount",.var=NULL,.setter=NULL,.getter=getEnqBufCount},
   {.name="ToPlay",.var=NULL,.setter=defaultSetter,.getter=getToPlay},
};


mxArray* getSource(biEntry* obj)
{
   if(obj->var==NULL || mxIsNumeric(obj->var))
      return mxCreateString("numeric");
   else
      return obj->var;
}

mxArray* getEnqBufCount(biEntry* obj)
{
   double retVal = 0.0;	
   biEntry* cmdStruct = lookupEntry("PageList",vars,sizeof(vars)/sizeof(*vars));
   if(cmdStruct->var!=NULL)
      retVal = mxGetNumberOfElements(cmdStruct->var);
  
   return mxCreateDoubleScalar(retVal);	
}

void incPageNo()
{
   biEntry* cmdStruct = lookupEntry("PageNo",vars,sizeof(vars)/sizeof(*vars));
   if(cmdStruct->var==NULL)
   {
      cmdStruct->var = mxCreateDoubleScalar(1.0);
      mexMakeArrayPersistent(cmdStruct->var);
   }
   else
   {
      double* pageNo = mxGetPr(cmdStruct->var);
      *pageNo = *pageNo + 1.0;
   }
}

mxArray* getToPlay(biEntry* obj)
{
   if(obj->var!=NULL)
   {  
      mxArray* duplicate = mxDuplicateArray(obj->var);	   
      mxDestroyArray(obj->var);
      obj->var = NULL;
      return duplicate;
   }
}



biEntry* lookupEntry(const char* name, biEntry* dict,size_t dictLen)
{	
   for(size_t ii=0;ii<dictLen;ii++)
   {
      if(!strncmp(name,dict[ii].name,strlen(dict[ii].name)))
      {
         return &dict[ii];
      }
   }
   return NULL;
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
   static int atExitRegistered = 0;
   if(!atExitRegistered)
   {
  //    mexAtExit();
   }   
   if(nrhs<1)
   {
     mexErrMsgTxt("BLOCK_INTERFACE: Not enough input arguments.");
   }

   mxArray* mxCmd = prhs[0];
   char command[COMMAND_LENGTH];
   size_t comStrLen = mxGetNumberOfElements(mxCmd);
   mxGetString(mxCmd,command,comStrLen+1);

   if(!strncmp(command,"set",3))
   {
      if(nrhs<2)
      {
         mexErrMsgTxt("BLOCK_INTERFACE: Not enough arguments for the set method.");
      }
      biEntry* cmdStruct = lookupEntry(command+3,vars,sizeof(vars)/sizeof(*vars));       
         
      if(cmdStruct==NULL)
      {
         mexErrMsgTxt("BLOCK_INTERFACE: Unrecognized set command.");
      }
      cmdStruct->setter(cmdStruct,prhs[1]);
      return;
   }
   else if(!strncmp(command,"get",3))
   {
      if(nrhs>1)
      {
         mexErrMsgTxt("BLOCK_INTERFACE: Too many input arguments for get method.");
      }
      biEntry* cmdStruct = lookupEntry(command+3,vars,sizeof(vars)/sizeof(*vars));       
      if(cmdStruct==NULL)
      {
         mexErrMsgTxt("BLOCK_INTERFACE: Unrecognized get command.");
      }
      plhs[0]=cmdStruct->getter(cmdStruct);
      return;
   }
   else if(!strcmp(command,"reset"))
   {
        


   } 
   else if(!strcmp(command,"incPageNo"))
   {
      incPageNo();
      return;
   } 

   mexErrMsgTxt("BLOCK_INTERFACE: Unrecognized command.");
 
}
