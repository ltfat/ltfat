function mwin=comp_nonsepwin2multi(g,a,M,s);
% Create multiwindow from non-sep win
  
L=size(g,1);

b=L/M;
mwin=zeros(L,s(2));
l=(0:L-1).'/L;
for ii=0:s(2)-1
  wavenum=mod(ii*s(1),s(2))*b/s(2);
  mwin(:,ii+1)=exp(2*pi*i*l*wavenum).*circshift(g,ii*a);
end;
