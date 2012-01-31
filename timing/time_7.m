%TIME_7
%
%   The purpose of this test is to compare the FB routines.
%
%   fb_1 is the original, using two passes through memory and no integer
%   optimizations.
%
%   fb_2 use only one pass and integer optimizations.


Lr=[480000*sf^2,480000*sf^2,262144*sf^2,262144*sf^2,900*sf^2];
ar=[     600*sf,     600*sf,        512,        512,       2];
Mr=[     800*sf,     800*sf,       1024,       1024,  600*sf];
gr=[     800*sf,  40*600*sf,       1024,     40*512,  600*sf];

for ii=1:length(Lr)

  L=Lr(ii);
  
  M=Mr(ii);
  a=ar(ii);

  gl=gr(ii);
  
  g=rand(L,1);
  f=rand(L,1);
  gfir=rand(gl,1);
  
  [L, a, M, gl]
  
  c1=mex_dgt_fb_1(f,gl,a,M,1);
  c2=mex_dgt_fb_2(f,gl,a,M,1);

end;


