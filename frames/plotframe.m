function outsig=plotframe(insig,F,varargin);
  
switch(F.type)
 case 'dgt'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[F.M,N,W]);

  plotdgt(insig,F.a,varargin);
 case 'dgtreal'
  [MN,W]=size(insig);
  M2=floor(F.M/2)+1;
  N=MN/M2;
  insig=reshape(insig,[M2,N,W]);
  plotdgtreal(insig,F.a,F.M,varargin);
 case 'dwilt'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[2*F.M,N/2,W]);
  plotdwilt(insig,varargin);
 case 'wmdct'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[2*F.M,N/2,W]);
  plotwmdct(insig,varargin);
end;

  