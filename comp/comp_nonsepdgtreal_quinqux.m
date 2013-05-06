function c=comp_nonsepdgtreal_quinqux(f,g,a,M)
%COMP_NONSEPDGTREAL_QUINQUX  Compute Non-separable Discrete Gabor transform
%   Usage:  c=comp_nonsepdgtreal_quinqux(f,g,a,M);
%
%   This is a computational subroutine, do not call it directly.

%   AUTHOR : Nicki Holighaus and Peter L. Søndergaard
%   TESTING: TEST_NONSEPDGT
%   REFERENCE: REF_NONSEPDGT

lt=[1 2];

L=size(f,1);
W=size(f,2);
N=L/a;
M2=floor(M/2)+1;

% ----- algorithm starts here, split into sub-lattices ---------------

c=zeros(M,N,W,assert_classname(f,g));

mwin=comp_nonsepwin2multi(g,a,M,[1 2],L);

% simple algorithm: split into sublattices

for ii=0:1
    c(:,ii+1:2:end,:)=comp_dgt(f,mwin(:,ii+1),2*a,M,[0 1],0,0,0);
end;

% Phase factor correction 
E = zeros(1,N,assert_classname(f,g));
for win=0:1
    for n=0:N/2-1
        E(win+n*2+1) = exp(-2*pi*i*a*n*rem(win,2)/M);
    end;
end;

c=bsxfun(@times,c(1:M2,:,:),E);
