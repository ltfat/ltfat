function gt=comp_gabtight_long(g,a,M);
%COMP_GABTIGHT_LONG  Compute tight window
%
%  This is a computational subroutine, do not call it directly, use
%  GABTIGHT instead.
%
%  See also: gabtight
  
L=length(g);
R=size(g,2);

gf=comp_wfac(g,a,M);

b=L/M;
N=L/a;

c=gcd(a,M);
d=gcd(b,N);
  
p=b/d;
q=N/d;

G=zeros(p,q*R,assert_classname(g));



if p==1
  % Integer oversampling, including Wilson basis.
  for ii=1:c*d
    gf(:,ii)=gf(:,ii)/norm(gf(:,ii));
  end;
else
  for ii=1:c*d
    
    G(:)=gf(:,ii);
    
    % Compute thin SVD
    [U,sv,V] = svd(G,'econ');
    
    Gtight=U*V';  
    
    gf(:,ii)=Gtight(:);
  end;
end;

gt=comp_iwfac(gf,L,a,M);

if isreal(g)
  gt=real(gt);
end;


