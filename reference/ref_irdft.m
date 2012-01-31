function c=ref_irdft(f)
%REF_IRDFT  Reference Inverse Real DFT
%   Usage:  c=ref_irdft(f);
%
%   Compute IRDFT by explicit formulas.
%
%   The transform is orthonormal

L=size(f,1);
Lhalf=ceil(L/2);
Lend=Lhalf*2-1;


F=zeros(L);

F(:,1)=ones(L,1);

l=(0:L-1).'/L;
for m=1:Lhalf-1
  F(:,2*m)=sqrt(2)*cos(2*pi*m*l);

  F(:,2*m+1)=sqrt(2)*sin(2*pi*m*l);

end;

if mod(L,2)==0
  F(:,L)=cos(pi*L*l);
end;

F=F/sqrt(L);

% dot-transpose will work because F is real.
c=F*f;


