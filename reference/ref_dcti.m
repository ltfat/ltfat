function c=ref_dcti(f)
%REF_DCTI  Reference Discrete Consine Transform type I
%   Usage:  c=ref_dcti(f);
%
%

L=size(f,1);
W=size(f,2);

if L==1
  % Doing the algorithm explicitly for L=1 does a division by
  % zero, so we exit here instead.
  c=f;
  return;
end;

% Create weights.
w=ones(L,1);
w(1)=1/sqrt(2);
w(L)=1/sqrt(2);

% Create transform matrix.
F=zeros(L);

for m=0:L-1
  for n=0:L-1
    F(m+1,n+1)=w(n+1)*w(m+1)*cos(pi*n*m/(L-1));
  end;
end;

F=F*sqrt(2/(L-1));
% Compute coefficients.
c=F'*f;


