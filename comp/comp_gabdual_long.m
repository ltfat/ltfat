function gd=comp_gabdual_long(g,a,M);
%COMP_GABDUAL_LONG  Compute dual window
%
%  This is a computational subroutine, do not call it directly, use
%  GABDUAL instead.
%
%  See also: gabdual
  
L=length(g);
R=size(g,2);

gf=comp_wfac(g,a,M);

b=L/M;
N=L/a;

c=gcd(a,M);
d=gcd(b,N);
  
p=b/d;
q=N/d;

gdf=zeros(p*q*R,c*d,assert_classname(g));

G=zeros(p,q*R,assert_classname(g));
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


