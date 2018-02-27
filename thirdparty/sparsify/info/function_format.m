function function_format
% Help file to explain how to create function files that calculate a
% linear operator and its transpose.
%
% You need two functions, one from the sparse domain to the observation
% domain and one form the observation domain to the sparse domain.
%
% Function format for operator form sparse domain to observation domain:
% Example:
%       x = function MyOp(s)
%       x = fft(s);
%       x = x(SubSet);
%
% For simplicity we do not assume the operator to require additional
% arguments. If you need to specify different arguments, create a function
% that takes arguments
% Example:
%       function x = MyOp_witharg(s,SubSet)
%       % For example arguments contains SubSet
%       S= 1/sqrt(sum(SubSet)) * fft(s);
%       x=[real(S(SubSet)) ; imag(S(SubSet))];
%
%
% Then specify a new function in the code using 
% Example:
% MyOp =@(z) MyOp_witharg(z,arguments)
%
% Then hand MyOp as function handle to any of the functions in the toolbox
%
% Do the same for the transpose operator:
% Example:
%       function s = MyOpTranspose(x)
%       s=zeros(N)
%       s(SubSet)=x;
%       s = ifft(s);
%
% Again, specify additional arguments by creating  
% Example:
%       function x = MyOpTranspose_witharg(x,SubSet)
%       % For example arguments is SubSet
%         lx=length(x);
%         s           = zeros(2*length(SubSet),1);
%         IndexSet    = logical([SubSet; zeros(size(SubSet))]);
%         s(IndexSet) = sqrt(sum(SubSet))*(x(1:lx/2) + i * x(lx/2+1:lx));
%         s           = 2 * 2 * length(SubSet) / lx * real(ifft(s));
%
% Then specify a new function in the code using 
% Example:
% MyOpTranspose =@(z) MyOpTranspose_witharg(z,arguments)
%
% Other methods are possible to specify function properties, such as using
% global variables, ....
%
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




