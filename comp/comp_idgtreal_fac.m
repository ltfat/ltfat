function f=comp_idgtreal_fac(coef,gf,L,a,M)
%COMP_IDGTREAL_FAC  Full-window factorization of a Gabor matrix assuming.
%   Usage:  f=comp_idgtreal_fac(c,gf,L,a,M)
%
%   Input parameters:
%         c     : M x N array of coefficients.
%         gf    : Factorization of window (from facgabm).
%         a     : Length of time shift.
%         M     : Number of frequency shifts.
%   Output parameters:
%         f     : Reconstructed signal.
%
%   Do not call this function directly, use IDGT.
%   This function does not check input parameters!
%
%   If input is a matrix, the transformation is applied to
%   each column.
%
%   This function does not handle multidimensional data, take care before
%   you call it.
%
%   References: so07-2 st98-8

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

% Calculate the parameters that was not specified.
N=L/a;
b=L/M;
M2=floor(M/2)+1;

R=prod(size(gf))/L;

%W=size(coef,2)/(N*R);
W = size(coef,3);

N=L/a;
b=L/M;

[c,h_a,h_m]=gcd(a,M);
h_a=-h_a;
p=a/c;
q=M/c;
d=N/q;

ff=zeros(p,q*W,c,d,assert_classname(coef,gf));
C=zeros(q*R,q*W,c,d,assert_classname(coef,gf));
f=zeros(L,W,assert_classname(coef,gf));

% Apply ifft to the coefficients.
coef=ifftreal(coef,M)*sqrt(M);
  
% Set up the small matrices

coef=reshape(coef,M,N,R,W);

if p==1

  for rw=0:R-1
    for w=0:W-1
      for s=0:d-1
	for l=0:q-1
	  for u=0:q-1
	    C(u+1+rw*q,l+1+w*q,:,s+1)=coef((1:c)+l*c,mod(u+s*q+l,N)+1,rw+1,w+1);
	  end;
	end;
      end;
    end;
  end;
else
  % Rational oversampling
  for rw=0:R-1
    for w=0:W-1
      for s=0:d-1
	for l=0:q-1
	  for u=0:q-1
	    C(u+1+rw*q,l+1+w*q,:,s+1)=coef((1:c)+l*c,mod(u+s*q-l*h_a,N)+1,rw+1,w+1);
	  end;
	end;
      end;
    end;
  end;
end;

% FFT them
if d>1
  C=fft(C,[],4);
end;

% Multiply them
for r=0:c-1    
  for s=0:d-1
    CM=reshape(C(:,:,r+1,s+1),q*R,q*W);
    GM=reshape(gf(:,r+s*c+1),p,q*R);

    ff(:,:,r+1,s+1)=GM*CM;
  end;
end;

% Inverse FFT
if d>1
  ff=ifft(ff,[],4);
end;

% Place the result  
if p==1

  for s=0:d-1
    for w=0:W-1
      for l=0:q-1
	f((1:c)+mod(s*M+l*a,L),w+1)=reshape(ff(1,l+1+w*q,:,s+1),c,1);
      end;
    end;
  end;

else
  % Rational oversampling
  for s=0:d-1
    for w=0:W-1
      for l=0:q-1
	for k=0:p-1
	  f((1:c)+mod(k*M+s*p*M-l*h_a*a,L),w+1)=reshape(ff(k+1,l+1+w*q,:,s+1),c,1);
	end;
      end;
    end;
  end;

end;

f=real(f);








