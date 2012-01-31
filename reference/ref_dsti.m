function c=ref_dsti(f)
%REF_DSTI  Reference Discrete Sine Transform type I
%   Usage:  c=ref_dsti(f);
%
%

L=size(f,1);
W=size(f,2);

% Create transform matrix.
F=zeros(L);

for m=0:L-1
  for n=0:L-1
    F(m+1,n+1)=sin(pi*(n+1)*(m+1)/(L+1));
  end;
end;

F=F*sqrt(2/(L+1));
% Compute coefficients.
c=F'*f;


