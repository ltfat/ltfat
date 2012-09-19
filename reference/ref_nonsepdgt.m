function [c]=ref_nonsepdgt(f,g,a,M,lt)
%REF_NONSEPDGT  Reference non-sep Discrete Gabor transform.
%   Usage:  c=ref_dgt(f,g,a,M,lt);
%
%   Linear algebra version of the algorithm. Create big matrix
%   containing all the basis functions and multiply with the transpose.

% Calculate the parameters that was not specified.
L=size(f,1);
g=fir2long(g,L);

b=L/M;
N=L/a;
W=size(f,2);
R=size(g,2);

% Create 2x2 grid matrix..
V=[a,0;...
   b/lt(2)*lt(1),b];
  
% Create lattice and Gabor matrix.
lat=ref_lattice(V,L);
G=ref_gaboratoms(g,lat);
  
% Apply matrix to f.
c=G'*f;

% reshape to correct output format.

c=reshape(c,M,N,R*W);


