function c=ref_dwiltiv(f,g,a,M)
%REF_DWILTIV   Reference DWILT type iv
%   Usage:  c=ref_dwiltiv(f,g,a,M);
%


L=size(g,1);

N=L/a;

F=zeros(L,M*N);

k=(0:L-1)';

pif=pi/4;
for n=0:floor(N/2)-1
  for m=0:2:M-1
    F(:,1+m+2*n*M)=sqrt(2)*circshift(g,2*n*a).*cos((m+.5)*pi*(k+.5)/M+pif);
    F(:,1+m+(2*n+1)*M)=sqrt(2)*circshift(g,(2*n+1)*a).*sin((m+.5)*pi*(k+.5)/M+pif);
  end;
  for m=1:2:M-1
    F(:,1+m+2*n*M)=sqrt(2)*circshift(g,2*n*a).*sin((m+.5)*pi*(k+.5)/M+pif);
    F(:,1+m+(2*n+1)*M)=sqrt(2)*circshift(g,(2*n+1)*a).*cos((m+.5)*pi*(k+.5)/M+pif);
  end;
end;

c=F.'*f;


