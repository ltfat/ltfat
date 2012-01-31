function f=ref_idwilt(c,g,a,M)
%REF_DWILT  Reference Inverse Discrete Wilson Transform
%   Usage:  f=ref_idwilt(c,g,a,M);
%

% Setup transformation matrix.

L=size(g,1);
W=size(c,2);
N=L/a;

F=zeros(L,M*N);



% Weight coefficients.

l=(0:L-1).';

pif=0;

if 1
  % This version uses sines and cosines to express the basis functions.

  for n=0:N/2-1    
    % Do the unmodulated coefficient.
    F(:,2*M*n+1)=circshift(g,2*a*n);
    
    % Setting this to -n*a should produce a time-invariant transform.
    timeinv=0; %-n*a;
    
    % m odd case
    for m=1:2:M-1
      F(:,m+2*M*n+1)   = sqrt(2)*sin(pi*m/M*(l+timeinv)+pif).*circshift(g,2*n*a);
      F(:,m+2*M*n+M+1) = sqrt(2)*cos(pi*m/M*(l+timeinv)+pif).*circshift(g,(2*n+1)*a);
    end;
    
    % m even case
    for m=2:2:M-1
      F(:,m+2*M*n+1)     = sqrt(2)*cos(pi*m/M*(l+timeinv)+pif).*circshift(g,2*n*a);
      F(:,m+2*M*n+M+1)   = sqrt(2)*sin(pi*m/M*(l+timeinv)+pif).*circshift(g,(2*n+1)*a);
    end;
    
    % Most modulated coefficient, Nyquest frequency.
    if mod(M,2)==0
      F(:,M+2*M*n+1)=(-1).^(l+timeinv).*circshift(g,2*n*a);
    else
      F(:,M+2*M*n+1)=(-1).^(l+timeinv).*circshift(g,(2*n+1)*a);
    end;
  end;

else

  % This version uses a cosine, 

  for n=0:N/2-1    
    % Do the unmodulated coefficient.
    F(:,2*M*n+1)=circshift(g,2*a*n);
    
    timeinv=-n*a;
    
    % m odd case
    for m=1:2:M-1
      F(:,m+2*M*n+1)   = sqrt(2)*cos(pi*m/M*(l+timeinv-M/2)+pif).*circshift(g,2*n*a);
      F(:,m+2*M*n+M+1) = sqrt(2)*sin(pi*m/M*(l+timeinv-M/2-a)+pif).*circshift(g,(2*n+1)*a);
    end;
    
    % m even case
    for m=2:2:M-1
      F(:,m+2*M*n+1)     = sqrt(2)*cos(pi*m/M*(l+timeinv-M/2)+pif).*circshift(g,2*n*a);
      F(:,m+2*M*n+M+1)   = sqrt(2)*sin(pi*m/M*(l+timeinv-M/2-a)+pif).*circshift(g,(2*n+1)*a);
    end;
    
    % Most modulated coefficient, Nyquest frequency.
    if mod(M,2)==0
      F(:,M+2*M*n+1)=(-1).^(l+timeinv).*circshift(g,2*n*a);
    else
      F(:,M+2*M*n+1)=(-1).^(l+timeinv-a).*circshift(g,(2*n+1)*a);
    end;
  end;



end;

f=F*c;

