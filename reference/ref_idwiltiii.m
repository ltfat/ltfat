function f=ref_idwiltiii(c,g,a,M)
%REF_IDWILTIII   Reference Inverse DWILT type III
%   Usage:  f=ref_idwiltiii(c,g,a,M);
%

L=size(g,1);
W=size(c,2);

N=L/a;

F=zeros(L,M*N);

l=(0:L-1)';

pif=pi/4;
for n=0:floor(N/2)-1
  for m=0:2:M-1
    F(:,1+m+2*n*M)=sqrt(2)*circshift(g,2*n*a).*cos((m+.5)*pi*l/M+pif);
    F(:,1+m+(2*n+1)*M)=sqrt(2)*circshift(g,(2*n+1)*a).*sin((m+.5)*pi*l/M+pif);
  end;
  for m=1:2:M-1
    F(:,1+m+2*n*M)=sqrt(2)*circshift(g,2*n*a).*sin((m+.5)*pi*l/M+pif);
    F(:,1+m+(2*n+1)*M)=sqrt(2)*circshift(g,(2*n+1)*a).*cos((m+.5)*pi*l/M+pif);
  end;
end;

f=F*c;

if 0
  pif=pi/4;
  for n=0:floor(N/2)-1
    for m=0:2:M-1
      %F(:,1+m+2*n*M)=sqrt(2)*circshift(g,2*n*a).*cos((m+.5)*(pi/M.*l-pi/2));
      %F(:,1+m+2*n*M)=sqrt(2)*circshift(g,2*n*a).*cos( ...
						      %       m*pi*l/M -m*pi/2+pi*l/(2*M)-pi/4);
      F(:,1+m+2*n*M)=sqrt(2)*circshift(g,2*n*a).*cos((m+.5)*pi*l/M+pif);
      
    end;
    for m=1:2:M-1
      %F(:,1+m+2*n*M)=sqrt(2)*circshift(g,2*n*a).*cos((m+.5)*(pi/M.*l-pi/2));
      %F(:,1+m+2*n*M)=sqrt(2)*circshift(g,2*n*a).*cos( ...
						      %       m*pi*l/M -m*pi/2+pi*l/(2*M)-pi/4);
      F(:,1+m+2*n*M)=sqrt(2)*circshift(g,2*n*a).*sin((m+.5)*pi*l/M+pif);
      
    end;
    %for m=0:M-1
      %  F(:,1+m+(2*n+1)*M)=sqrt(2)*circshift(g,(2*n+1)*a).*sin((m+.5)*(pi/M.*l-pi/2));
      %end;
    for m=0:2:M-1
      F(:,1+m+(2*n+1)*M)=sqrt(2)*circshift(g,(2*n+1)*a).*sin((m+.5)*pi*l/M+pif);
      
    end;
    for m=1:2:M-1
      F(:,1+m+(2*n+1)*M)=sqrt(2)*circshift(g,(2*n+1)*a).*cos((m+.5)*pi*l/M+pif);
      
    end;
    
  end;
end;


