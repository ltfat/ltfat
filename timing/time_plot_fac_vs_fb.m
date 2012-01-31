

if 0
  
  % Long example
  reps=10;
  
  a=600;
  M=800;
  
  K=10*lcm(a,M);
  W=2;

  ndatapoints=10;
  
  glfac=3;

end;

if 1
  
  % Short example
  reps=1000;
  
  a=20;
  M=30;
  
  K=lcm(a,M);
  W=2;
  
  glfac=3;

end;
  
  data_fac = zeros(ndatapoints,1);
  data_fb  = zeros(ndatapoints,1);

  x=K*(1:ndatapoints);
  

for ii=1:ndatapoints
  L=x(ii);
  L
  
  f=randn(L,W);
  g=pgauss(L,1);
  gc=middlepad(g,glfac*M);
  
  tic;
  for jj=1:reps
    c=dgt(f,g,a,M);
  end;
  data_fac(ii)=toc/reps;
  
  tic;
  for jj=1:reps
    c=dgt(f,gc,a,M);
  end;
  data_fb(ii)=toc/reps;
    
end;
  
plot(x,data_fac,'x',...
     x,data_fb,'o');

legend;

