function c=ref_ufilterbank(f,g,a);  
%REF_UFILTERBANK   Uniform filterbank by pconv
%   Usage:  c=ref_ufilterbank(f,g,a);
%
  
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

[a,M,longestfilter,lcm_a]=assert_filterbankinput(g,a,1);

[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,'UFILTERBANK',0);

L=ceil(max(Ls,longestfilter)/lcm_a)*lcm_a;

N=L/a;

c=zeros(N,M,W,assert_classname(f));
  
for w=1:W
  for m=1:M
    c(:,m,w)=pfilt(f(:,w),g{m},a);
  end;
end;

  

