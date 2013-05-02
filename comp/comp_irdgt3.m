function cout=comp_irdgtiii(cin,a,M)
%COMP_IRDGTIII  Compute inverse real DGT type III.
% 
%   This is a computational routine. Do not call it
%   directly.

%   AUTHOR : Peter L. Søndergaard

N=size(cin,1)/M;
W=size(cin,2);
L=N*a;

cin=reshape(cin,M,N,W);

Mhalf=floor(M/2);

cout=zeros(M,N,W,assert_classname(cin));

for m=0:Mhalf-1
  cout(m+1,:,:)=1/sqrt(2)*(cin(2*m+1,:,:)-i*cin(2*m+2,:,:));
  cout(M-m,:,:)=1/sqrt(2)*(cin(2*m+1,:,:)+i*cin(2*m+2,:,:));
end;

if mod(M,2)==1
  cout((M+1)/2,:,:)=cin(M,:,:);
end;





