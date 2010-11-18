function gd=comp_gabdual_long(g,a,M);
%COMP_GABDUAL_LONG  Compute dual window
%
%  This is a computational subroutine, do not call it directly, use
%  GABDUAL instead.
%
%  See also: gabdual
  
L=length(g);

gf=comp_wfac(g,a,M);

b=L/M;
N=L/a;

c=gcd(a,M);
d=gcd(b,N);
  
p=b/d;
q=N/d;

gdf=zeros(p*q,c*d);

G=zeros(p,q);
for ii=1:c*d
  % This essentially computes pinv of each block.

  G(:)=gf(:,ii);
  S=G*G';
  Gpinv=(S\G);

  gdf(:,ii)=Gpinv(:);
end;

gd=comp_iwfac(gdf,L,a,M);

if isreal(g)
  gd=real(gd);
end;
