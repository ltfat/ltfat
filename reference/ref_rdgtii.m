function c=ref_rdgtii(f,g,a,M)
%REF_RDGTII  Reference Real DGT type II
%   Usage:  c=ref_rdgt(f,g,a,M);
%
%   Linear algebra version of the algorithm. Create big matrix
%   containing all the basis functions and multiply with the transpose.


L=size(f,1);

b=L/M;
N=L/a;

Mhalf=ceil(M/2);


F=zeros(L,M*N);

l=(0:L-1).';

for n=0:N-1	 

  % Do the unmodulated coefficient.
  F(:,M*n+1)=circshift(g,n*a+floor(a/2));
  
  for m=1:Mhalf-1
    F(:,M*n+2*m)=sqrt(2)*cos(2*pi*m*(l+.5)/M).*circshift(g,n*a+floor(a/2));
    
    F(:,M*n+2*m+1)=sqrt(2)*sin(2*pi*m*(l+.5)/M).*circshift(g,n*a+floor(a/2));
    
  end;

  if mod(M,2)==0
    F(:,M*(n+1))=cos(pi*l).*circshift(g,n*a+floor(a/2));
  end;
  
end;

% dot-transpose will work because F is real.
c=F.'*f;


