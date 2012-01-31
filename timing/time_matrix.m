%TIME_MATRIX
%
%   The purpose of this test is to compare the fastest of the factorization
%   routines against the fastest of the filter bank routines for different
%   problems.
%
%   Some results: The FB computation seems to be completely memory bound:
%   increasing the length of the window does not impact the running time.
%
%   For the common 2xoversampling case with a long signal, the
%   factorization routine falls thorugh and is more than 10 times slower,
%   for the other case it is 3 times slower.


Lr=[480000*sf^2,480000*sf^2,262144*sf^2,262144*sf^2,900*sf^2,600, 600];
ar=[     600*sf,     600*sf,        512,        512,       2, 20,  20];
Mr=[     800*sf,     800*sf,       1024,       1024,  600*sf, 30,  40];
gr=[     800*sf,  40*600*sf,       1024,     40*512,  600*sf,400, 400];
Wr=[          2,          2,          2,          2,       1,600,   1];
nr=[          1,          1,          1,          1,       1,  1,1000];

for ii=1:length(Lr)

  L=Lr(ii);
  
  M=Mr(ii);
  a=ar(ii);

  gl=gr(ii);
  W=Wr(ii);
  nreps=nr(ii);
  
  g=rand(L,1);
  f=rand(L,W);
  gfir=rand(gl,1);
  
  gf=comp_wfac(g,a,M,1);
  
  tic;
  for jj=1:nreps
    c1=mex_dgt_fac_7(f,gf,a,M,0);
  end;
  t1=toc;
  
  % Flush the cache
  f=rand(L,W);
  
  tic;
    for jj=1:nreps
      c2=mex_dgt_fb_2(f,gl,a,M,0);
    end;
  t2=toc;

  %[L, W, a, M, gl]
  t1/t2
  
end;


