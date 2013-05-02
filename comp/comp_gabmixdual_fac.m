function gammaf=comp_gabmixdual_fac(gf1,gf2,L,a,M)
%COMP_GABMIXDUAL_FAC  Computes factorization of mix-dual.
%   Usage:  gammaf=comp_gabmixdual_fac(gf1,gf2,a,M)
%
%   Input parameters:
%      gf1    : Factorization of first window
%      gf2    : Factorization of second window
%      L      : Length of window.
%      a      : Length of time shift.
%      M      : Number of channels.
%
%   Output parameters:
%      gammaf : Factorization of mix-dual
%
%   GAMMAF is a factorization of a dual window of gf1
%
%   This function does not verify input parameters, call
%   GABMIXDUAL instead
%
%   See also:  gabmixdual, comp_fac, compute_ifac

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

LR=prod(size(gf1));
R=LR/L;

b=L/M;
N=L/a;

c=gcd(a,M);
d=gcd(b,N);
  
p=b/d;
q=N/d;

gammaf=zeros(p*q*R,c*d,assert_classname(gf1,gf2));

G1=zeros(p,q*R,assert_classname(gf1,gf2));
G2=zeros(p,q*R,assert_classname(gf1,gf2));
for ii=1:c*d

  G1(:)=gf1(:,ii);
  G2(:)=gf2(:,ii);
  S=G2*G1';
  Gpinv=M*S\G2;

  gammaf(:,ii)=Gpinv(:);
end;


