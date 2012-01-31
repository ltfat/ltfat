function c=ref_dwilt(f,g,a,M)
%REF_DWILT  Reference Discrete Wilson Transform
%   Usage:  c=ref_dwilt(f,g,a,M);
%
%   M must be even.

% Setup transformation matrix.

L=size(f,1);
N=L/a;

F=zeros(L,M*N);

% Zero-extend g if necessary
g=fir2long(g,L);

l=(0:L-1).';

for n=0:N/2-1    
  % Do the unmodulated coefficient.
  F(:,2*M*n+1)=circshift(g,2*n*a);

  % m odd case
  for m=1:2:M-1
    F(:,m+2*M*n+1)   = sqrt(2)*sin(pi*m/M*l).*circshift(g,2*n*a);
    F(:,m+2*M*n+M+1) = sqrt(2)*cos(pi*m/M*l).*circshift(g,(2*n+1)*a);
  end;

  % m even case
  for m=2:2:M-1
    F(:,m+2*M*n+1)     = sqrt(2)*cos(pi*m/M*l).*circshift(g,2*n*a);
    F(:,m+2*M*n+M+1)   = sqrt(2)*sin(pi*m/M*l).*circshift(g,(2*n+1)*a);
  end;

  % Most modulated coefficient, Nyquest frequency.
  if mod(M,2)==0
    F(:,M+2*M*n+1)=(-1).^l.*circshift(g,2*n*a);
  else
    F(:,M+2*M*n+1)=(-1).^l.*circshift(g,(2*n+1)*a);
  end;
end;


c=F'*f;

