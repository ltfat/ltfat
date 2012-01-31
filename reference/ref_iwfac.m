function [g]=ref_iwfac(gf,L,a,M)
%REF_IWFAC  Compute inverse window factorization
%   Usage: g=ref_iwfac(gf,L,a,M);
%
%   Input parameters:
%         gf    : Factored Window
%         L     : Length of transform.
%         a     : Length of time shift.
%         M     : Number of frequency bands.
%   Output parameters:
%         g     : Window function.
%

% Calculate the parameters that was not specified
R=prod(size(gf))/L;

N=L/a;
b=L/M;

% The four factorization parameters.
c=gcd(a,M);
p=a/c;
q=M/c;
d=N/q;

gf=reshape(gf,p,q*R,c,d);

% Scale by the sqrt(M) comming from Walnuts representation
gf=gf/sqrt(M);


% fft them
if d>1
  gf=ifft(gf,[],4);
end;

g=zeros(L,R);

% Set up the small matrices
for w=0:R-1
  for s=0:d-1
    for l=0:q-1
      for k=0:p-1
	g((1:c)+mod(k*M-l*a+s*p*M,L),w+1)=gf(k+1,l+1+q*w,:,s+1);
      end;
    end;
  end;
end;


