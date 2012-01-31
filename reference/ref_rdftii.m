function c=ref_rdftii(f)
%REF_RDFTII  Reference Real DFT type II
%   Usage:  c=ref_rdftii(f);
%
%   Compute RDFTII by explicit formulas.
%
%   The transform is orthonormal

L=size(f,1);
Lhalf=ceil(L/2);

F=zeros(L);

F(:,1)=ones(L,1);

l=(0:L-1).';
for m=1:Lhalf-1
  F(:,2*m)=sqrt(2)*cos(2*pi*m*(l+.5)/L);

  F(:,2*m+1)=sqrt(2)*sin(2*pi*m*(l+.5)/L);

end;

if mod(L,2)==0
  F(:,L)=cos(pi*l);
end;

F=F/sqrt(L);

% dot-transpose will work because F is real.
c=F.'*f;



