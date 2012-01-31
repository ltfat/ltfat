%TIME_4
%
%   The purpose of this test is to determine whether it pays of to do
%   the r-loop as the outmost loop to make the internal buffers smaller
%   and hopefully to reduce the demand on the memory subsystem.
%

Lr=[480000*sf^2,262144*sf^2,900*sf^2];
ar=[     600*sf,        512,       2];
Mr=[     800*sf,       1024,  600*sf];

for ii=1:length(Lr)

  L=Lr(ii);
  
  M=Mr(ii);
  a=ar(ii);
  
  [L, a, M]
  
  N=L/a;
  c=gcd(a,M);
  p=a/c;
  q=M/c;
  d=N/q;

  gf=comp_wfac(g,a,M);

  f=rand(L,1);
  gf=rand(p*q,c*d);  
  c1=mex_dgt_fac_1(f,gf,a,M,1);
  
  f=rand(L,1);
  gf=rand(p*q,c*d);  
  c2=mex_dgt_fac_4(f,gf,a,M,1);

end;

