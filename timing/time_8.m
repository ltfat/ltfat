%TIME_8
%
%   The purpose of this test is to determine whether it pays of to do
%   special optimizations for p=1 q=2 and p=2 q=3%


Lr=[480000*sf^2,262144*sf^2,900*sf^2];
ar=[     600*sf,        512,       2];
Mr=[     800*sf,       1024,  600*sf];

for ii=1:length(Lr)

  L=Lr(ii);
  
  M=Mr(ii);
  a=ar(ii);
  
  g=rand(L,1);
  f=rand(L,1);
  
  [L, a, M]
  
  gf=comp_wfac(g,a,M);
  
  c1=mex_dgt_fac_1(f,gf,a,M,1);

  % This is an attempt of cacheflushing
  f=rand(L,1); gf=rand(size(gf));

  c2=mex_dgt_fac_7(f,gf,a,M,1);

end;


