function [c]=ref_dgtns(f,gamma,V)
%REF_DGT  Reference Discrete Gabor transform for non-separable lattices.
%   Usage:  c=ref_dgtns(f,gamma,V);
%
%   Linear algebra version of the algorithm. Create big matrix
%   containing all the basis functions and multiply with the transfpose.

% Calculate the parameters that was not specified.
L=size(gamma,1);

% Create lattice and Gabor matrix.
lat=ref_lattice(V,L);
G=ref_gaboratoms(gamma,lat);
  
% Apply matrix to f.
c=G'*f;


