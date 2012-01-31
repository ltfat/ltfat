function c=ref_gdgt(f,g,a,M,c_t,c_f,c_w)
%REF_GDGT  Reference generalized DGT
%   Usage:  c=ref_dgtiv(f,g,a,M,c_t,c_f,c_w);
%
%   Linear algebra version of the algorithm. Create big matrix
%   containing all the basis functions and multiply with the transpose.


L=size(f,1);

N=L/a;

F=zeros(L,M*N);

l=(0:L-1).';

for n=0:N-1	   
  for m=0:M-1
    F(:,M*n+m+1)=exp(2*pi*i*(m+c_f)*(l+c_t)/M).*circshift(g,n*a+c_w);
  end;
end;

c=F'*f;




