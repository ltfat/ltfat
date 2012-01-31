function f=ref_irdgt3_1(c,g,a,M)
%REF_IRDGT3_1  Reference Inverse Real DGT type 3 by IDGT3
%   Usage:  f=ref_irdgt3_1(c,g,a,M);
%
%   Compute a complex coefficient layout for IDGT3

L=size(g,1);
W=size(c,2);
R=size(g,2);

b=L/M;
N=L/a;

Mhalf=floor(M/2);

c=reshape(c,M,N,W);
cc=zeros(M,N,W);

for m=0:Mhalf-1
  cc(m+1,:,:)=1/sqrt(2)*(c(2*m+1,:,:)-i*c(2*m+2,:,:));
  cc(M-m,:,:)=1/sqrt(2)*(c(2*m+1,:,:)+i*c(2*m+2,:,:));
end;

if mod(M,2)==1
  cc((M+1)/2,:,:)=c(M,:,:);
end;

cc=reshape(cc,M*N,W);

f=ref_igdgt(cc,g,a,M,0,.5,0);





