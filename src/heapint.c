#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ltfat.h"

#define PI 3.1415926535897932384626433832795

typedef struct 
{
      int *h;
      int heapsize;
      int totalheapsize;
      const LTFAT_REAL *s;
}
LTFAT_NAME(heap);

void LTFAT_NAME(heap_insert)(LTFAT_NAME(heap) *h, const int key)
{
   int pos,pos2, swap;
   
   /*printf("Grow heap: heapsize %i, totalheapsize %i\n",h->heapsize,h->totalheapsize);*/

   /* Grow heap if necessary */
   if (h->totalheapsize==h->heapsize)
   {
      (h->totalheapsize)*=2;
      h->h = (int*)realloc(h->h,h->totalheapsize*sizeof(int));
   }  

   pos=h->heapsize;
   h->heapsize++;

   /*printf("heap: heapsize %i, pos %i\n",h->heapsize,pos);*/

   h->h[h->heapsize-1]=key;  
     
   /* printf("before heapint main loop\n"); */

   while (pos>0)
   {
      /* printf("pos %i\n",pos); */
      pos2=(pos-pos%2)/2;
      if (h->s[h->h[pos2]] < h->s[h->h[pos]] )
      {
	 swap=h->h[pos2];
	 h->h[pos2]=h->h[pos];
	 h->h[pos]=swap;
	 pos=pos2;
      }
      else
      {
	 break;
      }
   }
   
   /* { */
/*       printf("C INSERT heapsize: %i\n",h->heapsize); */
/*       for (ii=0; ii<h->heapsize; ii++) */
/* 	 printf("    %i\n",h->h[ii]); */
/*    } */   

}

int LTFAT_NAME(heap_delete)(LTFAT_NAME(heap) *h)
{
   
   int pos, maxchildpos, swap, key;
   LTFAT_REAL maxchildkey;
   
   /* Extract first element */
   key=h->h[0];
   
   /* Put last element on first elements place, and make the heap smaller. */
   h->h[0]=h->h[h->heapsize-1];
   h->heapsize--;
   
   /* Fix the just introduced problem. */
   pos=0;
   
   /*  %%%%%%%%%%%%%%
    * %
    * %  Is maxchildpos 0 or 1 indexed!
    * %
    * %
    * %%%%%%%%%%%%%
    */


   while (2*pos+1<h->heapsize)
   {
      if (2*pos+3>h->heapsize)
      {
	 maxchildkey=h->s[h->h[2*pos+1]];
	 maxchildpos=1;
      }
      else    
      {
	 if (h->s[h->h[2*pos+1]]>=h->s[h->h[2*pos+2]])
	 {
	    maxchildkey=h->s[h->h[2*pos+1]];
	    maxchildpos=1;
	 }
	 else
	 {
	    maxchildkey=h->s[h->h[2*pos+2]];
	    maxchildpos=2;
	 }
      }
      
      if (maxchildkey>h->s[h->h[pos]])
      {
	 swap=h->h[2*pos+maxchildpos];
	 h->h[2*pos+maxchildpos]=h->h[pos];
	 h->h[pos]=swap;
	 pos=2*pos+maxchildpos;
      }
      else
      {
	 break;
      }
   }

/*    { */
/*       printf("C DELETE heapsize: %i\n",h->heapsize); */
/*       for (ii=0; ii<h->heapsize; ii++) */
/* 	 printf("    %i %f\n",h->h[ii],h->s[h->h[ii]]); */
/*    }    */


   return key;
}



void LTFAT_NAME(heapint)(const LTFAT_REAL *s,
			 const LTFAT_REAL *tgrad,
			 const LTFAT_REAL *fgrad,
			 const int a, const int M, const int L, const int W,
			 LTFAT_REAL tol, LTFAT_REAL *phase)
{

  /* Declarations */
  int b, N, ii, jj, Imax, domainloop;
  int w, w_E, w_W, w_N, w_S;
  LTFAT_REAL maxs;
  int *donemask;
  LTFAT_NAME(heap) h;

  N=L/a;
  b=L/M;

  /* Main body */
  /* printf("Main body a: %i, M: %i, N, %i, L: %i, W: %i\n",a,M,N,L,W); */


  h.totalheapsize  =(int)(M*log((LTFAT_REAL)M));
  h.h              = (int*)ltfat_malloc(h.totalheapsize*sizeof(int));
  h.s              = s;
  h.heapsize       = 0;
    
  donemask = (int*)ltfat_malloc(M*N*W*sizeof(int));

  LTFAT_REAL *tgradw = (LTFAT_REAL*)ltfat_malloc(M*N*sizeof(LTFAT_REAL));
  LTFAT_REAL *fgradw = (LTFAT_REAL*)ltfat_malloc(M*N*sizeof(LTFAT_REAL));

  /* Set the phase to zero initially */
  for (ii=0; ii<M*N*W; ii++)
  {
     phase[ii]=0.0;
  }

  /* Rescale the partial derivatives to a torus of length L. We copy
     to work arrays because the input is const. */
  for (jj=0; jj<N; jj++)
  {
     for (ii=0; ii<M; ii++)
     {
       tgradw[ii+jj*M] = tgrad[ii+jj*M]*2*PI*a/L;
       fgradw[ii+jj*M] = fgrad[ii+jj*M]*2*PI*b/L+2*PI*jj*a*b/L;
     }
  }

  maxs=-1e99;
  Imax=0; /* This is just included to avoid a warning from the C++ compiler */
  for (ii=0; ii<M*N; ii++)
  {
     if (s[ii]>maxs)
     {
	maxs=s[ii];
	Imax=ii;
     }
  }

/*   printf("Mark\n"); */

  /* Mark all the small elements as done, they get a zero phase.
   * Code 5
   */
  for (ii=0; ii<M*N*W; ii++)
  {
     if (s[ii]<tol*maxs)
	donemask[ii]=5;
     else
	donemask[ii]=0;
  }
  

  domainloop=1;  


  while (domainloop)
  {
     /* printf("Main loop, Imax: %i\n",Imax); */

     h.heapsize=0;
  
     /* Put maximal element onto the heap and mark that it is done. It
      * will get zero phase.
      */
     LTFAT_NAME(heap_insert)(&h,Imax);
     
     /* printf("after heap_insert, Imax %i, h.heapize %i\n", Imax,h.heapsize); */

     donemask[Imax]=6;

/*      printf("after donemask\n",h.heapsize); */

     while (h.heapsize>0)
     {
	
	/* printf("Inner loop, heapsize: %i\n",h.heapsize); */
	/* Extract largest element from heap and delete it. */
	
       w = LTFAT_NAME(heap_delete)(&h);

	/* Try and put the four neighbours onto the heap.
	 * Integration by trapezoidal rule */

	/* North */
	w_N=(w-1+M)%M+w-w%M;
	if (!donemask[w_N])
	{
	 /*   printf("C North: %i %i\n",w,w_N); */
	   phase[w_N] = phase[w]+(fgradw[w]+fgradw[w_N])/2;
	   donemask[w_N]=1;
	   LTFAT_NAME(heap_insert)(&h,w_N);
	}

	/* South */
	w_S=((w+1)%M)+w-w%M;
	if (!donemask[w_S])
	{
/* 	   printf("C South: %i %i\n",w,w_S); */
	   phase[w_S] = phase[w]-(fgradw[w]+fgradw[w_S])/2;
	   donemask[w_S]=2;
	   LTFAT_NAME(heap_insert)(&h,w_S);
	}

	/* East */
	w_E = (w+M)%(M*N);
	if (!donemask[w_E])
	{
/* 	   printf("C East:%i %i\n",w,w_E); */
	   phase[w_E]=phase[w]+(tgradw[w]+tgradw[w_E])/2;
	   donemask[w_E]=3;
	   LTFAT_NAME(heap_insert)(&h,w_E);
	}

	/* West */
	w_W = (w-M+M*N)%(M*N);
	if (!donemask[w_W])
	{
/* 	   printf("C West: %i %i\n",w,w_W); */
	   phase[w_W]=phase[w]-(tgradw[w]+tgradw[w_W])/2;
	   donemask[w_W]=4;
	   LTFAT_NAME(heap_insert)(&h,w_W);
	}
            
     }

     /* Find the new maximal element 
      *
      * This code is currently missing  
      *
      */

     maxs=-1e99;
     domainloop=0;
     for (ii=0; ii<M*N; ii++)
     {
	if (!donemask[ii] && s[ii]>maxs)
	{
	   maxs=s[ii];
	   Imax=ii;
	   domainloop=1;
	}
     }


  }
  
  ltfat_free(tgradw);
  ltfat_free(fgradw);
  ltfat_free(donemask);
  ltfat_free(h.h);

}
      
