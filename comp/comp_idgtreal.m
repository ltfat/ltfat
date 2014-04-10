function f=comp_idgtreal(coef,g,a,M,lt,phasetype)
%COMP_IDGTREAL  Compute IDGTREAL
%   Usage:  f=comp_idgtreal(c,g,a,M,lt,phasetype);
%
%   Input parameters:
%         c     : Array of coefficients.
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of modulations.
%         lt    : lattice type
%   Output parameters:
%         f     : Signal.
%

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: TEST_DGT
%   REFERENCE: OK

N=size(coef,2);
L=N*a;

if lt(2)==1
    f = comp_isepdgtreal(coef,g,L,a,M,phasetype);
else
    % Quinqux lattice
    f=comp_inonsepdgtreal_quinqux(coef,g,a,M);
end;

