#include "config.h"
#include "ltfat.h"



void fftindex(const ltfatInt N, ltfatInt *indexout)
{
    ltfatInt ii;

    if (N%2==0)
    {
        for (ii=0; ii<N/2+1; ii++)
        {
            indexout[ii]=ii;
        }
        for (ii=N/2; ii<N-1; ii++)
        {
            indexout[ii+1]=-N+ii+1;
        }
    }
    else
    {
        for (ii=0; ii<(N-1)/2+1; ii++)
        {
            indexout[ii]=ii;
        }
        for (ii=(N-1)/2; ii<N-1; ii++)
        {
            indexout[ii+1]=-N+ii+1;
        }
    }

}

ltfatInt imax(const ltfatInt a, const ltfatInt b)
{
   return (a > b ? a : b);
}

ltfatInt imin(const ltfatInt a, const ltfatInt b)
{
   return (a < b ? a : b);
}

/*
#define MAXFNC(T,prefix,suffix) \
T prefix##max##suffix(const T a, const T b) \
{                                    \
   return (a > b ? a : b);           \
}

MAXFNC(size_t,,_st)
MAXFNC(ptrdiff_t,,_pt)
MAXFNC(ltfatInt,int_,)

#undef MAXFNC

#define MINFNC(T,prefix,suffix) \
T prefix##min##suffix(const T a, const T b) \
{                                    \
   return (a < b ? a : b);           \
}

MINFNC(size_t,,_st)
MINFNC(ptrdiff_t,,_pt)
MINFNC(ltfatInt,int_,)

#undef MINFNC
*/

ltfatInt makelarger(const ltfatInt L, const ltfatInt K)
{
    /* This is a floor operation */
    ltfatInt o = (L/K)*K;

    /* Make it a ceil */
    if (L%K>0)
    {
        o += K;
    }

    return o;
}

/* Extended Euclid algorithm. */
ltfatInt gcd (const ltfatInt a, const ltfatInt b, ltfatInt *r, ltfatInt *s )
{
    ltfatInt a1 = a;
    ltfatInt b1 = b;
    ltfatInt a2 = 1;
    ltfatInt b2 = 0;
    ltfatInt a3 = 0;
    ltfatInt b3 = 1;
    ltfatInt c, d;
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

ltfatInt lcm(const ltfatInt a, const ltfatInt b)
{
    ltfatInt junk_r, junk_s;

    ltfatInt c = gcd(a, b, &junk_r, &junk_s);

    return (a*b/c);
}



void gabimagepars(const ltfatInt Ls, const ltfatInt x, const ltfatInt y,
                  ltfatInt *a, ltfatInt *M, ltfatInt *L, ltfatInt *N, ltfatInt *Ngood)
{


    *M = imin(y,Ls);
    *N = imax(x,Ls);

    /* Determine the minimum transform size. */
    ltfatInt K = lcm(*M,*N);

    /* This L is good, but is it not the same as DGT will choose. */
    ltfatInt Llong = makelarger(Ls,K);

    /* Fix a from the long L */
    *a=Llong/(*N);

    /* Now we have fixed a and M, so we can use the standard method of choosing L. */
    ltfatInt Lsmallest=lcm(*a,*M);
    *L = makelarger(Ls, Lsmallest);

    /* We did not get N as desired. */
    *N=*L/(*a);

    /* Number of columns to display */
    *Ngood=(Ls/(*a));
}

/* Determine the size of the output array of wfacreal and iwfacreal */
ltfatInt wfacreal_size(const ltfatInt L, const ltfatInt a, const ltfatInt M)
{

    ltfatInt h_a, h_m;

    const ltfatInt b=L/M;
    const ltfatInt c=gcd(a, M,&h_a, &h_m);
    const ltfatInt p=a/c;
    const ltfatInt d=b/p;

    /* This is a floor operation. */
    const ltfatInt d2= d/2+1;

    return d2*p*M;

}

ltfatInt nextPow2(ltfatInt x)
{
    ltfatInt bits = sizeof(x)*8;

    if(x==0)
        return 1;

    x--;
    (x) = ((x)>>1)  | (x);
    (x) = ((x)>>2)  | (x);
    (x) = ((x)>>4)  | (x);
    (x) = ((x)>>8)  | (x);
    (x) = ((x)>>16) | (x);
    if(bits>32)
        (x) = ((x)>>32) | (x);

    (x)++;
    return x;
}

ltfatInt nextfastfft(const ltfatInt x)
{
   ltfatInt xtmp = x;
    while (1)
    {
        ltfatInt m = xtmp;

        while ((m % 2) == 0)
            m /= 2;
        while ((m % 3) == 0)
            m /= 3;
        while ((m % 5) == 0)
            m /= 5;
        if (m <= 1)
            break;                    /* n is completely factorable by twos, threes, and fives */
        xtmp++;
    }
    return xtmp;
}


ltfatInt pow2(const ltfatInt x)
{
    return ((1)<<(x));
}

ltfatInt modPow2(const ltfatInt x, const ltfatInt pow2var)
{
    return ((x)&(pow2var-1));
}

int isPow2(ltfatInt x)
{
    return x==nextPow2(x);
}

ltfatInt ilog2(const ltfatInt x)
{
    ltfatInt tmp = 0;
    ltfatInt xtmp = x;
     while (xtmp >>= 1) ++tmp;
    return tmp;
}

// integer power by squaring
ltfatInt ipow(ltfatInt base, ltfatInt exp)
{
    ltfatInt result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}


ltfatInt filterbank_td_size(ltfatInt L, ltfatInt a, ltfatInt gl, ltfatInt offset, ltfatExtType ext)
{
    ltfatInt Lc = 0;
    if(ext==PER)
    {
        Lc = (ltfatInt) ceil( L/((double)a) );
    }
    else if(ext==VALID)
    {
        Lc = (ltfatInt) ceil( (L-(gl-1))/((double)a) );

    }
    else
    {
        Lc = (ltfatInt) ceil( (L + gl - 1 + offset )/((double)a) );
    }
    return Lc;
}


