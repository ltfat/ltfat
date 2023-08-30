function c=comp_dgtreal(f,g,a,M,lt,phasetype)
%COMP_DGTREAL  Compute a DGTREAL
%   Usage:  c=comp_dgtreal(f,g,a,M,lt,phasetype);
%
%   Input parameters:
%         f     : Input data
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of modulations.
%         L     : Length of transform to do.
%   Output parameters:
%         c     : M/2+1*N array of coefficients.
%

%   AUTHOR : Peter L. Søndergaard.

if lt(2)==1
    c = comp_sepdgtreal(f,g,a,M,phasetype);
else
    % Quinqunx lattice
    c = comp_nonsepdgtreal_quinqux(f,g,a,M);            
end;





