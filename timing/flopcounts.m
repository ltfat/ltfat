function [fcfac,fcfb]=flopcounts(a,M,L,Lg)
%
%   Compute the flop count for the factorization algorithm.
  
N=L/a;

[c,h_a,h_m]=gcd(a,M);
h_a=-h_a;
p=a/c;
q=M/c;
d=N/q;

  
fcfac = L*8*q+4*L*(1+q/p)*log2(d)+4*M*N*log2(M);

fcfb  = 8*L*Lg/a+4*M*N*log2(M);
