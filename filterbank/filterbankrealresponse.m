function gf=filterbankrealresponse(g,a,L)
  
  L2=floor(L/2)+1;

  gf=filterbankresponse(g,a,L);
  gf=gf+involute(gf);
  gf=gf(1:L2);
  
  

%OLDFORMAT
