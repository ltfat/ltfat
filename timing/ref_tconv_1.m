function h=ref_tconv_1(f,g,a)
%REF_TCONV_1  Reference TCONV
%   Usage:  h=ref_tconv_1(f,g)
%
%   REF_TCONV_1(f,g,a) computes the twisted convolution of f and g.
%
%   Version for sparse matrices without precomputation of the 
%   position of non-zeros coefficients.

% AUTHOR: Floret Jaillet

L=size(f,1);

h=zeros(L,L);


% precompute the Lth roots of unity
% Optimization note : the special properties and symetries of the 
% roots of unitycould be exploited to reduce this computation.
% Furthermore here we precompute every possible root if some are 
% unneeded. 
temp=exp((-i*2*pi/L)*(0:L-1)');
[rowf,colf,valf]=find(f);
[rowg,colg,valg]=find(g);

h=sparse(L,L);

for indf=1:length(valf)
  for indg=1:length(valg)
    m=mod(rowf(indf)+rowg(indg)-2, L);
    n=mod(colf(indf)+colg(indg)-2, L);
    h(m+1,n+1)=h(m+1,n+1)+valf(indf)*valg(indg)*temp(mod((m-(rowf(indf)-1))*(colf(indf)-1),L)+1);
  end
end



