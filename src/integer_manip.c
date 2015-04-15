#include "ltfat.h"


LTFAT_EXTERN_TOO void
fftindex(const ltfatInt N, ltfatInt *indexout)
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

LTFAT_EXTERN_TOO ltfatInt
imax(const ltfatInt a, const ltfatInt b)
{
   return (a > b ? a : b);
}

LTFAT_EXTERN_TOO ltfatInt
imin(const ltfatInt a, const ltfatInt b)
{
   return (a < b ? a : b);
}

LTFAT_EXTERN_TOO ltfatInt
makelarger(const ltfatInt L, const ltfatInt K)
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
LTFAT_EXTERN_TOO ltfatInt
gcd (const ltfatInt a, const ltfatInt b, ltfatInt *r, ltfatInt *s )
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

LTFAT_EXTERN_TOO ltfatInt
lcm(const ltfatInt a, const ltfatInt b)
{
    ltfatInt junk_r, junk_s;

    ltfatInt c = gcd(a, b, &junk_r, &junk_s);

    return (a*b/c);
}


LTFAT_EXTERN_TOO void
gabimagepars(const ltfatInt Ls, const ltfatInt x, const ltfatInt y,
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
LTFAT_EXTERN_TOO ltfatInt
wfacreal_size(const ltfatInt L, const ltfatInt a, const ltfatInt M)
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

LTFAT_EXTERN_TOO ltfatInt
nextPow2(const ltfatInt y)
{
    ltfatInt x = (ltfatInt) y;
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

LTFAT_EXTERN_TOO ltfatInt
nextfastfft(const ltfatInt x)
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

LTFAT_EXTERN_TOO ltfatInt
pow2(const ltfatInt x)
{
    return ((1)<<(x));
}

LTFAT_EXTERN_TOO ltfatInt
modPow2(const ltfatInt x, const ltfatInt pow2var)
{
    return ((x)&(pow2var-1));
}

LTFAT_EXTERN_TOO int
isPow2(const ltfatInt x)
{
    return x==nextPow2(x);
}

LTFAT_EXTERN_TOO int
ilog2(const ltfatInt x)
{
    ltfatInt tmp = 0;
    ltfatInt xtmp = x;
     while (xtmp >>= 1) ++tmp;
    return tmp;
}

// integer power by squaring
LTFAT_EXTERN_TOO ltfatInt
ipow(const ltfatInt base, const ltfatInt exp)
{
    ltfatInt baseTmp = (ltfatInt) base;
    ltfatInt expTmp = (ltfatInt) exp;
    ltfatInt result = 1;

    while (expTmp)
    {
        if (expTmp & 1)
            result *= baseTmp;
        expTmp >>= 1;
        baseTmp *= baseTmp;
    }

    return result;
}

LTFAT_EXTERN_TOO ltfatInt
filterbank_td_size(const ltfatInt L, const ltfatInt a, const ltfatInt gl,
                   const ltfatInt offset, const ltfatExtType ext)
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

LTFAT_EXTERN_TOO ltfatInt
ltfat_round(const double x)
{
    if (x < 0.0)
        return (ltfatInt)(x - 0.5);
    else
        return (ltfatInt)(x + 0.5);
}

LTFAT_EXTERN_TOO ltfatInt
positiverem(const ltfatInt a,const ltfatInt b)
{
    const ltfatInt c = a%b;
    return(c<0 ? c+b : c);
}

LTFAT_EXTERN_TOO ltfatInt
rangelimit(const ltfatInt a, const ltfatInt amin, const ltfatInt amax)
{
    ltfatInt c = a < amin? amin:a;
    c = c > amax? amax: c;
    return c;
}

