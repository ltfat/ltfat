function h=ref_tconv(f,g,a)
%REF_PCONV  Reference TCONV
%   Usage:  h=ref_tconv(f,g,a)
%
%   TCONV(f,g,a) computes the twisted convolution of f and g.

% AUTHOR: Peter L. Søndergaard

M=size(f,1);
N=size(f,2);

h=zeros(M,N);

theta=a/M;

for m=0:M-1
  for n=0:N-1
    for l=0:N-1
      for k=0:M-1
	h(m+1,n+1)=h(m+1,n+1)+f(k+1,l+1)*g(mod(m-k,M)+1,mod(n-l,N)+1)*...
	    exp(2*pi*i*theta*(m-k)*l);
      end;
    end;
  end;
end;


