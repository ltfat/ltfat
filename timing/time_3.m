%TIME_3
%
%   The purpose of this test is to determine whether it pays of to do
%   contiquous FFTs in stead of strided FFTs
%

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

  N=L/a;
  c=gcd(a,M);
  p=a/c;
  q=M/c;
  d=N/q;

  f=rand(L,1);
  gf=rand(p*q,c*d);  
  c1=mex_dgt_fac_1(f,gf,a,M,1);
  
  f=rand(L,1);
  gf=rand(p*q,c*d);  
  c2=mex_dgt_fac_3(f,gf,a,M,1);

end;

