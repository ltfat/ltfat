function lambdas=comp_gfeigs(gf,L,a,M)
%COMP_GFEIGS_SEP
%   Usage:  lambdas=comp_gfeigs(gf,a,M);
%
%   Compute Eigenvalues of a Gabor frame operator in
%   the separable case.
%
%   This is a computational routine, do not call it directly.
%
%   See help on GFBOUNDS

%   AUTHOR : Peter L. Søndergaard.

LR=prod(size(gf));
R=LR/L;

b=L/M;
N=L/a;

c=gcd(a,M);
d=gcd(b,N);
p=b/d;
q=N/d;

% Initialize eigenvalues
AF=Inf;
BF=0;

% Holds subsubmatrix.
C=zeros(p,q*R,assert_classname(gf));

lambdas=zeros(p,c*d,assert_classname(gf));

% Iterate through all the subsubmatrices.
for k=0:c*d-1
  % Extract p x q*R matrix of array.
  C(:)=gf(:,k+1);
  
  % Get eigenvalues of 'squared' subsubmatrix.
  lambdas(:,1+k)=eig(C*C');
    
end;

% Clean eigenvalues, they are real, and
% scale them correctly.
lambdas=real(lambdas);

% Reshape and sort.
lambdas=sort(lambdas(:));









