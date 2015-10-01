function [c,newphase,tgrad,fgrad]=constructphase(s,g,a,tol)
%CONSTRUCTPHASE  Construct the phase of a DGT
%   Usage:  c=constructphase(s,g,a);
%           c=constructphase(s,g,a,tol);
%
%   `constructphase(s,g,a)` will construct a suitable phase for the postive
%   valued coefficients *s*. 
%  
%   If *s* is the absolute values of the Gabor coefficients of a signal
%   obtained using the window *g* and time-shift *a*, i.e.:
%
%     c=dgt(f,g,a,M);
%     s=abs(c);
%
%   then `constuctphase(s,g,a)` will attempt to reconstruct *c*.
%
%   The window *g* must be Gaussian, i.e. *g* must have the value `'gauss'`
%   or be a cell array `{'gauss',tfr}`. 
%
%   `constructphase(s,g,a,tol)` does as above, but sets the phase of
%   coefficients less than *tol* to random values. 
%   By default, *tol* has the value 1e-10.
%
%   This function requires a computational subroutine that is only
%   available in C. Use |ltfatmex| to compile it.
%
%   See also:  dgt, gabphasegrad, ltfatmex
%

% AUTHOR: Peter L. Søndergaard, Zdenek Prusa

if nargin<4
    tol=1e-10;
else
    if ~iscalar(tol)
        error('%s: tol must be scalar.',upper(mfilename));
    end
end

% Compute phase gradients, check parameteres
[tgrad,fgrad] = gabphasegrad('abs',s,g,a,2);
% Build the phase (calling a MEX file)
newphase=comp_heapint(s,tgrad,fgrad,a,tol);
% Set phase of the coefficients below tol to random values
absthr = max(s(:))*tol;
toosmallidx = s<absthr;
zerono = numel(find(toosmallidx));
newphase(toosmallidx) = rand(zerono,1)*2*pi;
% Combine the magnitude and phase
c=s.*exp(1i*newphase);

