function c=ref_edgtii(f,g,a,M)
%REF_EDGTII   Reference Even Discrete Gabor transform type II
%   Usage  c=ref_edgt(f,g,a,M);
%
%   If a is even, then the input window must be odd-centered of length 2L.
%   
%   If a is odd, then the input window must be even-centered of length 2L.
%

L=size(f,1);
W=size(f,2);

N=L/a;
b=L/M;


%Fold=zeros(2*L,M*N);
F=zeros(L,M*N);

%l2=(0:2*L-1).';
l=(0:L-1).';

%lfirst=(0:L-1).';
lsecond=(2*L-1:-1:L).';

gnew=circshift(g,floor(a/2));
for n=0:N-1	   
  for m=0:M-1
    %Fold(:,M*n+m+1)=exp(2*pi*i*m*(l2+.5)/M).*circshift(gnew,n*a);

    gshift=circshift(gnew,n*a);
    %F(:,M*n+m+1)=exp(2*pi*i*m*(lfirst+.5)/M).*gshift(l+1) + ...
%	exp(2*pi*i*m*(lsecond+.5)/M).*gshift(lsecond+1);

    F(:,M*n+m+1)=exp(2*pi*i*m*(l+.5)/M).*gshift(l+1) + ...
	exp(-2*pi*i*m*(l+.5)/M).*gshift(lsecond+1);

  end;
end;

c=F'*f;


