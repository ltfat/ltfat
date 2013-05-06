function gd=ref_gabtight(g,a,M)
%REF_GABTIGHT   Reference GABTIGHT
%   Usage:  gd=ref_gabtight(g,a,M);
%
%   Calculate the canonical tight window by simple linear algebra
g = double(g);

G=tfmat('dgt',g,a,M);

gd=(G*G')^(-1/2)*g;

