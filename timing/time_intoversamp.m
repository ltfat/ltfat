%TIME_INTOVERSAMP
%
%   The purpose of this test is to determine whether special code for
%   handling the integer oversampling case in the factorization routines
%   pays off in any way.


Lr=[1048576*2,2400];
ar=[      512,   2];
Mr=[     1024, 800];

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

  f=rand(L,1);
  gf=rand(p*q,c*d);  
  c1=mex_dgt_fac_1(f,gf,a,M);

  f=rand(L,1);
  gf=rand(p*q,c*d);
  c2=mex_dgt_fac_2(f,gf,a,M);

end;

