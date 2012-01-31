function f=ref_iedgtii(c,g,a,M)
%REF_IEDGTII   Inverse Reference Even DGT type II
%   Usage  f=ref_edgt(c,g,a,M);
%
%   If a is even, then the input window must be odd-centered of length 2L.
%   
%   If a is odd, then the input window must be even-centered of length 2L.
%

L2=size(g,1);
W=size(c,2);
L=L2/2;

b=L/M;
N=L/a;

F=zeros(L,M*N);

l=(0:L-1).';

lsecond=(2*L-1:-1:L).';

gnew=circshift(g,floor(a/2));
for n=0:N-1	   
  for m=0:M-1
    gshift=circshift(gnew,n*a);

    F(:,M*n+m+1)=exp(2*pi*i*m*(l+.5)/M).*gshift(l+1) + ...
	exp(-2*pi*i*m*(l+.5)/M).*gshift(lsecond+1);

  end;
end;

f=F*c;


