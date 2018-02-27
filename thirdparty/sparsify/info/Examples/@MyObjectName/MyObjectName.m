function  out = MyObjectName(SubSet)
% MyObjectName linear operator class constructor
% operator_property1 = Some property (e.g. size of operator)
% operator_property2 = Some other property (e.g. random subset of a matrix)
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
      out.adjoint = 0; 
      out.SubSet  = SubSet;


  %     Create the Object
      out = class(out,'MyObjectName');