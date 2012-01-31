function c=ref_irdftiii(f)
%REF_IRDFTIII  Reference Inverse Real DFT type III
%   Usage:  c=ref_irdftiii(f);
%
%   Compute IRDFTIII by explicit formulas.
%
%   The transform is orthonormal

L=size(f,1);
Lhalf=floor(L/2);

F=zeros(L);

F(:,1)=ones(L,1);

l=(0:L-1).'/L;
for m=0:Lhalf-1
  F(:,2*m+1)=sqrt(2)*cos(2*pi*(m+.5)*l);

  F(:,2*m+2)=sqrt(2)*sin(2*pi*(m+.5)*l);
end;

if mod(L,2)==1
  F(:,L)=cos(pi*L*l);
end;

F=F/sqrt(L);

% dot-transpose will work because F is real.
c=F*f;



