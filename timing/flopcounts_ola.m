function [fcfacola,fcola]=flopcounts(a,M,L,Lg,Lb)
%
%   Compute the flop count for the factorization algorithm.
  
N=L/a;

[c,h_a,h_m]=gcd(a,M);
h_a=-h_a;
p=a/c;
q=M/c;
d=N/q;

Nb=L/Lb;
rho=(Lb+Lg)/Lb;
  
fcfacola  = rho*(8*L*q+4*L*(1+q/p)*log2(d*rho)+4*M*N*log2(M));

fcola     = rho*(8*L*M+4*L*log2(rho*Lb)+4*M*N*log2(rho*N));


