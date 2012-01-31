function c=ref_dcti_1(f)
%REF_DCTI_1  Reference Discrete Consine Transform type I
%   Usage:  c=ref_dcti_1(f);
%
%

L=size(f,1);
W=size(f,2);

if L==1
  c=f;
  return;
end;

R=1/sqrt(2)*[eye(L);...
	     [zeros(L-2,1),flipud(eye(L-2)),zeros(L-2,1)]];

R(1,1)=1;
R(L,L)=1;

c=R'*dft(R*f);



