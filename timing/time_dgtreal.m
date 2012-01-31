%TIME_DGTREAL  Time the DGT versus DGTREAL
%

a=1000;
M=1500;

L=a*M;
W=4;

f=randn(L,W);

g=randn(L,1);

tic;
  c1 = dgt(f,g,a,M);
  
t1=toc;

disp(sprintf('Time to execute dgt:     %f',t1));

tic;
  c2 = dgtreal(f,g,a,M);
t2=toc;

disp(sprintf('Time to execute dgtreal: %f',t2));

