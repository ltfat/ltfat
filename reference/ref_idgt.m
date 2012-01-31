function f=ref_idgt(c,g,a,M)
%REF_DGT  Reference Inverse Discrete Gabor transform.
%   Usage:  c=ref_idgt(f,g,a,M);
%
%   Linear algebra version of the algorithm. Create big matrix
%   containing all the basis functions and multiply with the transpose.

% Calculate the parameters that was not specified.
L=size(g,1);

b=L/M;
N=L/a;
W=size(c,2);

% Create 2x2 grid matrix..
V=[a,0;
   0,b];
  
% Create lattice and Gabor matrix.
lat=ref_lattice(V,L);
G=ref_gaboratoms(g,lat);
  
% Apply matrix to c.
f=G*c;


