function c=comp_nonsepdgt_multi(f,g,a,M,lt)
%COMP_NONSEPDGT_MULTI  Compute Non-separable Discrete Gabor transform
%   Usage:  c=comp_nonsepdgt_multi(f,g,a,M,lt);
%
%   This is a computational subroutine, do not call it directly.

%   AUTHOR : Nicki Holighaus and Peter L. Soendergaard
%   TESTING: TEST_NONSEPDGT
%   REFERENCE: REF_NONSEPDGT

% Assert correct input.

L=size(f,1);
W=size(f,2);
N=L/a;

% ----- algorithm starts here, split into sub-lattices ---------------

c=zeros(M,N,W);

mwin=comp_nonsepwin2multi(g,a,M,lt);

% simple algorithm: split into sublattices

for ii=0:lt(2)-1
    c(:,ii+1:lt(2):end,:)=comp_dgt(f,mwin(:,ii+1),lt(2)*a,M,L,0);
end;

% Phase factor correction 
E = zeros(1,N);
for win=0:lt(2)-1
    for n=0:N/lt(2)-1
        E(win+n*lt(2)+1) = exp(-2*pi*i*a*n*rem(win*lt(1),lt(2))/M);
    end;
end;
E

for w=1:W
    c(:,:,w) = c(:,:,w).*repmat(E,M,1);
end;
