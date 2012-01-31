function c=ref_edgt(f,g,a,M)
%REF_EDGT   Reference Even Discrete Gabor transform
%   Usage  c=ref_edgt(f,g,a,M);
%
%   The input window must be odd-centered.

L=size(f,1);
W=size(f,2);

N=L/a;
M=L/b;

F=zeros(L,M*N);

l=(0:L-1)';
for n=0:N-1
  for m=0:M-1
    F(:,1+m+n*M)=exp(2*pi*i*m.*(l+.5)*b/L).*circshift(g,n*a);
  end;
end;

c=F'*f;


