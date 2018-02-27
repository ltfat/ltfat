function s=MyOpTranspose_witharg(x,SubSet)
%
% Copyright (c) 2007 Thomas Blumensath
%
% The University of Edinburgh
% Email: thomas.blumensath@ed.ac.uk
% Comments and bug reports welcome
%
% This file is part of sparsity Version 0.1
% Created: April 2007
%
% Part of this toolbox was developed with the support of EPSRC Grant
% D000246/1
%
% Please read COPYRIGHT.m for terms and conditions.

lx=length(x);
if rem(lx/2,1);
   error('We require input to MyOpTranspose_witharg to be of even length.') 
end

s           = zeros(2*length(SubSet),1);
IndexSet    = logical([SubSet; zeros(size(SubSet))]);
s(IndexSet) = sqrt(sum(SubSet))*(x(1:lx/2) + i * x(lx/2+1:lx));
s           = 2 * 2*length(SubSet) / lx * real(ifft(s));

