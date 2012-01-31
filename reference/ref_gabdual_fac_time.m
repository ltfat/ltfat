function ff=ref_gabdual_fac_time(gf,L,a,M)
%REF_GABDUAL_FAC_TIME  Computes factorization of canonical dual window 
%

LR=prod(size(gf));
R=LR/L;

b=L/M;
N=L/a;

c=gcd(a,M);
d=gcd(b,N);
  
p=b/d;
q=N/d;

ff=zeros(p*q*R,c*d);

G=zeros(p,q*R);
for ii=1:c*d
  % This essentially computes pinv of each block.

  G(:)=gf(:,ii);
  S=G*G';
  Gpinv=(S\G);

  ff(:,ii)=Gpinv(:);
end;



