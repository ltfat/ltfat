function [f]=ref_isfac(ff,L,a,M)
%REF_ISFAC  Reference inverse signal factorization
%   Usage: f=ref_sfac(ff,a,M);
%
%   Input parameters:
%         ff    : Factored signal
%         a     : Length of time shift.
%         b     : Length of frequency shift.
%   Output parameters:
%         f     : Output signal.
%

% Calculate the parameters that was not specified
W=prod(size(ff))/L;

N=L/a;
M=L/b;

% The four factorization parameters.
[c,h_a,h_m]=gcd(a,M);
p=a/c;
q=M/c;
d=N/q;

permutation=zeros(q*b,1);
P=stridep(p,b);

% Create permutation
for l=0:q-1
  for s=0:b-1
    permutation(l*b+1+s)=mod(P(s+1)-1-h_m*l,b)*M+l*c+1;
  end;
end;

f=ref_ifac(ff,W,c,d,p,q,permutation);


