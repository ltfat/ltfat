function outsig=iframet(insig,F);
  
switch(F.type)
 case 'dgt'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[F.M,N,W]);
  outsig=idgt(insig,F.gd,F.a);
 case 'dgtreal'
  [MN,W]=size(insig);
  M2=floor(F.M/2)+1;
  N=MN/M2;
  insig=reshape(insig,[M2,N,W]);
  outsig=idgtreal(insig,F.gd,F.a,F.M);
 case 'dwilt'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[2*F.M,N/2,W]);
  outsig=idwilt(insig,F.gd);
 case 'wmdct'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[F.M,N,W]);
  outsig=iwmdct(insig,F.gd);
end;

  