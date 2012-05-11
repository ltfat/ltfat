 
R=2;
L=24;

for wintype = 1:3
  switch wintype
   case 1
    a=6;
    M=4;
    g=pherm(24,0:R-1);
   case 2
    a=4;
    M=6;
    g=crand(L,R);
   case 3
    a=6;
    M=4;
    g=[pgauss(L,1),circshift(pgauss(L),a/2)];
  end;
  
  N=L/a;
  
  gd=gabdual(g,a,M);
  
  f=crand(24,1);
  r=zeros(L,1);
  c=zeros(M,N,R);
  for ii=1:R
    c=dgt(f,g(:,ii),a,M);
    r=r+idgt(c,gd(:,ii),a);
  end;
 
  norm(f-r)
  
end;