function  out = ctranspose(P)
% Function to toggle the transpose property of object P of class
% MyObjectName
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
    P.adjoint = xor(P.adjoint,1);
    out = P;