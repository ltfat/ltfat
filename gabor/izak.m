function f=izak(c);
%IZAK  Inverse Zak transform
%   Usage:  f=izak(c);
%
%   `izak(c)` computes the inverse Zak transform of *c*. The parameter of
%   the Zak transform is deduced from the size of *c*.
%
%   See also:  zak
%
%   References: ja94-4 bohl97-1

%   AUTHOR : Peter L. Søndergaard
%   TESTING: TEST_ZAK
%   REFERENCE: OK

complainif_argnonotinrange(nargin,1,1,mfilename);

a=size(c,1);
N=size(c,2);
W=size(c,3);

L=a*N;

% Create output matrix.
f=zeros(L,W,assert_classname(c));

for ii=1:W
  % Iterate through third dimension of c.
  % We use a normalized DFT, as this gives the correct normalization
  % of the Zak transform.
  f(:,ii)=reshape(idft(c(:,:,ii),[],2),L,1);
end;
