function f=ref_idwiltii(c,g,a,M)
%REF_DWILTII  Reference Inverse Discrete Wilson Transform type II
%   Usage:  f=ref_idwiltii(c,g,M);
%

% Setup transformation matrix.

L=size(g,1);
W=size(c,2);

F=zeros(L,L);

N=L/a;

l=(0:L-1).';

for n=0:N/2-1    
  % Do the unmodulated coefficient.
  F(:,2*M*n+1)=circshift(g,2*a*n); 
  
  % m odd case OK
  for m=1:2:M-1
    F(:,m+2*M*n+1) = sqrt(2)*sin(pi*m/M*(l+.5)).*circshift(g,2*n*a);
    F(:,m+2*M*n+M+1) = sqrt(2)*cos(pi*m/M*(l+.5)).*circshift(g,(2*n+1)*a);
  end;

  for m=2:2:M-1
    F(:,m+2*M*n+1) = sqrt(2)*cos(pi*m/M*(l+.5)).*circshift(g,2*n*a);
    F(:,m+2*M*n+M+1) = sqrt(2)*sin(pi*m/M*(l+.5)).*circshift(g,(2*n+1)*a);
  end;

  if rem(M,2)==0
    % Do the Nyquest wave.
    F(:,M+2*M*n+1)   = sin(pi*(l+.5)).*circshift(g,(2*n+1)*a); 
  else
    F(:,M+2*M*n+1)   = sin(pi*(l+.5)).*circshift(g,2*n*a); 
  end;     
  
end;

f=F*c;

