function [c]=ref_nonsepdgt(f,g,a,M)
%REF_NONSEPDGT  Reference non-sep Discrete Gabor transform.
%   Usage:  c=ref_dgt(f,g,a,M,s);
%
%   Linear algebra version of the algorithm. Create big matrix
%   containing all the basis functions and multiply with the transpose.

% Calculate the parameters that was not specified.
L=size(g,1);

b=L/M;
N=L/a;
W=size(f,2);
R=size(g,2);

% Create 2x2 grid matrix..
V=[a,0;
   b/s(2)*s(1),b];
  
% Create lattice and Gabor matrix.
lat=ref_lattice(V,L);
G=ref_gaboratoms(g,lat);
  
% Apply matrix to f.
c=G'*f;

% reshape to correct output format.

c=reshape(c,M,N,R*W);


