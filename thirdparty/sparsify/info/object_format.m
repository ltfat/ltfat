function object_format
% Help file to explain how to create an object to calculate linear
% operator and their transpose.
% The same object style is used in [1].
%
% To create your own objects you need to:
%    a) Define a new object class
%    b) Overload two operations for this object, multiplication (mtimes) 
%       and transpose (ctranspose)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% STEP BY STEP GUIDE:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   1) decide on an object name (say MyObjectName)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   2) create a folder called @MyObjectName
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   3) create an m file for class constructor called MyObjectName.m 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example:
%
%       function  out = MyObjectName(operator_property1,operator_property2,....)
%       % MyObjectName linear operator class constructor
%       % operator_property1 = Some property (e.g. size of operator)
%       % operator_property2 = Some other property (e.g. random subset of a matrix)
%       out.adjoint = 0; 
%       out.perator_property1 = perator_property1;
%       out.perator_property2 = perator_property2;
%
%       %     Create the Object
%       out = class(out,'MyObjectName');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   4) create m file called mtimes.m within @MyObjectName
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example:
%
%   function out = mtimes(P,x)
%   % This function specifies which happens if we write P*x
%   % in matlab where P is an object from class MyObjectName
%   if P.adjoint == 0 %calculate P*x
%       % Specify what shall happen here
%       % for example choose a subset of the fft, where the subset is
%       % specified by object_propert1
%           S           = 1/sqrt(sum(P.SubSet)) * fft(x);
%           out         = [real(S(P.SubSet)) ; imag(S(P.SubSet))];
% 
%   else %Calculate the transpose operator P'*x
%       % Specify what shall happen here
%       % for example choose the transpose of the above operation, 
%       % where the subset is specified by object_propert1
%     
%       lx          = length(x);
%       s           = zeros(2*length(P.SubSet),1);
%       IndexSet    = logical([P.SubSet; zeros(size(P.SubSet))]);
%       s(IndexSet) = sqrt(sum(P.SubSet))*(x(1:lx/2) + i * x(lx/2+1:lx));
%       out         = 2 * 2 * length(P.SubSet) / lx * real(ifft(s));
% 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   5) create m file called ctranspose.m within @MyObjectName
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example:
%
%       function  out = ctranspose(P)
%       % Function to toggle the transpose property of object P of class
%       % MyObjectName
%       P.adjoint = xor(P.adjoint,1);
%       out = P;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   6) Create an object in your code giving it the desired properties
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example:
%
%       P  = MyObjectName(operator_property1,operator_property2,....);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   7) Now you can use P and P' in the following way
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example:
%       x=P*s;
%       s=P'*x;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% References 
%   [1] Kwangmoo Koh, Seung-Jean Kim, and Stephen Boyd, l1_ls: Simple 
%       Matlab Solver for l1-regularized Least Squares Problems
%       http://www.stanford.edu/~boyd/l1_ls/
%
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