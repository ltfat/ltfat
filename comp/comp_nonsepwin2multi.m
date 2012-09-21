function mwin=comp_nonsepwin2multi(g,a,M,lt,L);
% Create multiwindow from non-sep win
  
Lw=size(g,1);

b=L/M;
mwin=zeros(L,lt(2));
l=(0:L-1).'/L;
for ii=0:lt(2)-1
  wavenum=mod(ii*lt(1),lt(2))*b/lt(2);
  mwin(:,ii+1)=exp(2*pi*i*l*wavenum).*circshift(g,ii*a);
end;
