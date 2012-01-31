function f=ref_irdgt3(c,g,a,M)
%REF_IRDGT3  Reference Inverse Real DGT type 3
%   Usage:  c=ref_irdgt3(f,g,a,M);
%
%   Linear algebra version of the algorithm. Create big matrix
%   containing all the basis functions and multiply with the transpose.


L=size(g,1);

b=L/M;
N=L/a;

Mhalf=floor(M/2);

F=zeros(L,M*N);

l=(0:L-1).'/L;

for n=0:N-1
  
  for m=0:Mhalf-1
    F(:,M*n+2*m+1)=sqrt(2)*cos(2*pi*(m+.5)*b*l).*circshift(g,n*a);
    
    F(:,M*n+2*m+2)=sqrt(2)*sin(2*pi*(m+.5)*b*l).*circshift(g,n*a);
    
  end;

  if mod(M,2)==1
    F(:,M*(n+1))=cos(pi*L*l).*circshift(g,n*a);
  end;
  
end;

f=F*c;


