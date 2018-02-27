function x = MyOp_witharg(s,SubSet)
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

ls=length(s);
if rem(ls/2,1);
   error('We require input to MyOp_witharg to be of even length.') 
end

S= 1/sqrt(sum(SubSet)) * fft(s);
x=[real(S(SubSet)) ; imag(S(SubSet))];
