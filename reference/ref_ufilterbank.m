function c=ref_ufilterbank(f,g,a);  
%REF_UFILTERBANK   Uniform filterbank by pconv
%   Usage:  c=ref_ufilterbank(f,g,a);
%
  
L=size(f,1);

N=L/a;

c=zeros(N,M,W,assert_classname(f));
  
for w=1:W
  for m=1:M
    c(:,m,w)=pfilt(f(:,w),g{m},a);
  end;
end;

  

