function f=ref_irdgt(c,g,a,M)
%REF_IRDGT  Reference Inverse Real DGT
%   Usage:  c=ref_rdgt(f,g,a,M);
%
%   Linear algebra version of the algorithm. Create big matrix
%   containing all the basis functions and multiply with the transpose.


L=size(g,1);

b=L/M;
N=L/a;

Mhalf=ceil(M/2);


F=zeros(L,M*N);

l=(0:L-1).';

for n=0:N-1

  % Do the unmodulated coefficient.
  F(:,M*n+1)=circshift(g,n*a);
  
  for m=1:Mhalf-1
    F(:,M*n+2*m)=sqrt(2)*cos(2*pi*m*l/M).*circshift(g,n*a);;
    
    F(:,M*n+2*m+1)=sqrt(2)*sin(2*pi*m*l/M).*circshift(g,n*a);;
    
  end;

  if mod(M,2)==0
    F(:,M*(n+1))=cos(pi*l).*circshift(g,n*a);;
  end;
  
end;

% dot-transpose will work because F is real.
f=F*c;


