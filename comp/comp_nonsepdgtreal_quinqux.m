function c=comp_nonsepdgt_multi(f,g,a,M)
%COMP_NONSEPDGT_MULTI  Compute Non-separable Discrete Gabor transform
%   Usage:  c=comp_nonsepdgt_multi(f,g,a,M);
%
%   This is a computational subroutine, do not call it directly.

%   AUTHOR : Nicki Holighaus and Peter L. Soendergaard
%   TESTING: TEST_NONSEPDGT
%   REFERENCE: REF_NONSEPDGT

lt=[1 2];

L=size(f,1);
W=size(f,2);
N=L/a;
M2=floor(M/2)+1;

% ----- algorithm starts here, split into sub-lattices ---------------

c=zeros(M2,N,W);

mwin=comp_nonsepwin2multi(g,a,M,lt,L);

% simple algorithm: split into sublattices

for ii=0:1
    c(:,ii+1:lt(2):end,:)=comp_dgtreal(f,mwin(:,ii+1),lt(2)*a,M,[0 1],0);
end;

% Phase factor correction 
E = zeros(1,N);
for win=0:lt(2)-1
    for n=0:N/lt(2)-1
        E(win+n*lt(2)+1) = exp(-2*pi*i*a*n*rem(win,2)/M);
    end;
end;

c=bsxfun(@times,c,E);
