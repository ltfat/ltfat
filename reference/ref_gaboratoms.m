function G=ref_gaboratoms(g,spoints);
%REF_GABORATOMS  Create Gabor transformation matrix
%   Usage:  G=ref_gaboratoms(g,spoints);
%
%   Given a set a TF-sampling points, as returned
%   by REF_GABORLATTICE, returns the corresponding
%   Gabor matrix. Each column of the output matriorx is
%   a Gabor atom.
%

L=size(g,1);
W=size(g,2);
MN=size(spoints,2);

% Create Gabor matrix.
G=zeros(L,MN*W);
jj=(0:L-1).';

% Calculate atoms from sampling points.
for w=0:W-1
  for p=1:MN;
    G(:,p+w*MN)=exp(2*pi*i*spoints(2,p)*jj/L).*circshift(g(:,w+1),spoints(1,p));
  end;
end;




