function out = mtimes(P,x)
% This function specifies what happens if we write P*x
% in matlab where P is an object from class MyObjectName
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
% Please read COPYRIGHT.m for terms and conditions.
if P.adjoint == 0 %calculate P*x
  % Specify what shall happen here
  % for example choose a subset of the fft, where the subset is
  % specified by object_propert1
    s           = x;
    ls          = length(s);
    if rem(ls/2,1);
    error('We require input to MyOp_witharg to be of even length.') 
    end

    S           = 1/sqrt(sum(P.SubSet)) * fft(s);
    out         = [real(S(P.SubSet)) ; imag(S(P.SubSet))];

else %Calculate the transpose operator P'*x
  % Specify what shall happen here
  % for example choose the transpose of the above operation, 
  % where the subset is specified by object_propert1
    
    lx          = length(x);
    if rem(lx/2,1);
       error('We require input to MyOpTranspose_witharg to be of even length.') 
    end

    s           = zeros(2*length(P.SubSet),1);
    IndexSet    = logical([P.SubSet; zeros(size(P.SubSet))]);
    s(IndexSet) = sqrt(sum(P.SubSet))*(x(1:lx/2) + i * x(lx/2+1:lx));
    out         = 2 * 2 * length(P.SubSet) / lx * real(ifft(s));

end