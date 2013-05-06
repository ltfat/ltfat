function c=ref_dsti_1(f)
%REF_DSTI_1  Reference Discrete Sine Transform type I
%   Usage:  c=ref_dsti_1(f);
%
%

L=size(f,1);
W=size(f,2);

if L==1
  c=f;
  return;
end;

R=1/sqrt(2)*[zeros(1,L,assert_classname(f));...
	     eye(L);
	     zeros(1,L,assert_classname(f));...
	     -flipud(cast(eye(L),assert_classname(f)))];

c=i*R'*dft(R*f);



