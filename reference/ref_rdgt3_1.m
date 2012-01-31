function c=ref_rdgt3_1(f,g,a,M)
%REF_RDGT3_1  Reference Real DGT type 3
%   Usage:  c=ref_rdgt3_1(f,g,a,M);
%
%   Compute a DGT3 and pick out the coefficients

L=size(f,1);
W=size(f,2);
R=size(g,2);

N=L/a;

Mhalf=floor(M/2);

cc=ref_gdgt(f,g,a,M,0,.5,0);
cc=reshape(cc,M,N,W);

c=zeros(M,N,W);

for m=0:Mhalf-1
  c(2*m+1,:,:)=sqrt(2)*real(cc(m+1,:,:));
  c(2*m+2,:,:)=-sqrt(2)*imag(cc(m+1,:,:));
end;

if mod(M,2)==1
  c(M,:,:)=cc((M+1)/2,:,:);
end;

c=reshape(c,M*N,W);


