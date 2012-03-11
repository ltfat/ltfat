function outsig=framet(insig,F);
  
switch(F.type)
 case 'dgt'
  outsig=dgt(insig,F.g,F.a,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
 case 'dgtreal'
  outsig=dgtreal(insig,F.g,F.a,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
 case 'dwilt'
  outsig=dwilt(insig,F.g,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
 case 'wmdct'
  outsig=wmdct(insig,F.g,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
end;

  