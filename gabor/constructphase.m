function c=constructphase(s,g,a)
%CONSTRUCTPHASE  Construct the phase of a DGT
%   Usage:  c=constructphase(s,g,a);
%           c=constructphase(s,g,a,tol);
%
%   CONSTRUCTPHASE(s,g,a) will construct a suitable phase for the postive
%   valued coefficients s. 
%  
%   If s is the absolute values of the Gabor coefficients of a signal
%   obtained using the window g and time-shift _a, i.e.:
%
%C     c=dgt(f,g,a,M);
%C     s=abs(c);
%
%   then CONSTUCTPHASE(s,g,a) will attempt to reconstruct _c.
%
%   The window g must be Gaussian, i.e. g must have the value 'gauss' or
%   be a cell array {'gauss',tfr}
%
%   CONSTRUCTPHASE(s,g,a,tol) does as above, but sets the phase of
%   coefficients less than tol to zero phase. This speeds up the
%   computation. By default, tol has the value 1e-10.
%
%   This function requires a computational subroutine that is only
%   available in C. Use LTFATMEX to compile it.
%
%   See also:  dgt, gabphasegrad, ltfatmex
%
%   Demos:  demo_constructphase
  
% AUTHOR: Peter Soendergaard
M=size(s,1);
N=size(s,2);
L=N*a;
  
% Obtain the vectors.
[tgrad,fgrad] = gabphasegrad('abs',s,g,a);

tol=1e-10;

newphase=comp_heapint(s,tgrad,fgrad,a,tol);
%newphase=heapint(s,tgrad,fgrad,a,tol);

c=s.*exp(i*newphase);


