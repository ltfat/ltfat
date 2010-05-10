#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int int_max(const int a, const int b)
{
   return (a > b ? a : b);
}

int int_min(const int a, const int b)
{
   return (a < b ? a : b);
}

int makelarger(const int L, const int K)
{
   /* This is a floor operation */
   int o = (L/K)*K;

   /* Make it a ceil */
   if (L%K>0)
   {
      o += K;
   }

   return o;   
}

/* Extended Euclid algorithm. */
int gcd (const int a, const int b, int *r, int *s )
{
  int a1 = a;
  int b1 = b;
  int a2 = 1;
  int b2 = 0;
  int a3 = 0;
  int b3 = 1;
  int c, d;
  while ( b1 != 0 ) 
  {
      d=a1/b1;
      c = a1;
      a1 = b1;
      b1 = c-d*b1;

      c = a2;
      a2 = b2;
      b2 = c-d*b2;
      
      c = a3;
      a3 = b3;
      b3 = c-d*b3;
      
  }
   
  *r=a2;
  *s=a3;
  return a1;
}

int lcm(const int a, const int b)
{
  int junk_r, junk_s;

  int c = gcd(a, b, &junk_r, &junk_s);

  return (a*b/c);
}



void gabimagepars(const int Ls, const int x, const int y,
		  int *a, int *M, int *L, int *N, int *Ngood)
{


  *M = int_min(y,Ls);
  *N = int_max(x,Ls);
  
  /* Determine the minimum transform size. */
  int K = lcm(*M,*N);

  /* This L is good, but is it not the same as DGT will choose. */
  int Llong = makelarger(Ls,K);

  /* Fix a from the long L */
  *a=Llong/(*N);
  
  /* Now we have fixed a and M, so we can use the standard method of choosing L. */
  int Lsmallest=lcm(*a,*M);
  *L = makelarger(Ls, Lsmallest);

  /* We did not get N as desired. */
  *N=*L/(*a);
	 
  /* Number of columns to display */
  *Ngood=(Ls/(*a));
}
