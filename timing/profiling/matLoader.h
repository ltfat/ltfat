#include "mat.h"
#include <string.h>

int getNoOfArgsFromMAT(const char* file)
{
  MATFile *pmat;
  const char **dir;
  int	  ndir;

  pmat = matOpen(file, "r");
  if (pmat == NULL) {
    printf("Error opening file %s\n", file);
    return(1);
  }

  dir = (const char **)matGetDir(pmat, &ndir);
  mxFree(dir);

  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(1);
  }
  return ndir;
}

int getArgsFromMAT(const char* file, mxArray* prhs[], int ndir)
{
  MATFile *pmat;
  const char *name;
  int	  ii;
  mxArray *pa;

  pmat = matOpen(file, "r");
  if (pmat == NULL) {
    printf("Error reopening file %s\n", file);
	return;
  }

  /* Read in each array. */
  //printf("\nReading in the actual array contents:\n");
  for (ii=0; ii<ndir; ii++) {
      pa = matGetNextVariable(pmat, &name);
	  if(strncmp(name,"arg",3))
	  {
	     printf("Variable name shoud be in the following format: arg## \n");
		 return -1;
	  }

	  char buf[4];

	  memcpy(buf,name+3,strlen(name)-3);
	  int idx = atoi(buf);

	  if(!(idx<ndir))
	  {
	     printf("Variable name shoud be in the following format: arg##, where ## is argument number less than total arg. count. \n");
		 return -1;
	  }

      if (pa == NULL) {
	  printf("Error reading in file %s\n", file);
	  return -1;
      }
	  prhs[idx] = pa;
  }

  if (matClose(pmat) != 0) {
      printf("Error closing file %s\n",file);
	  return -1;
  }
  return 0;
}
