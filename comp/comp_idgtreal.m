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
b=L/M;



M2=floor(M/2)+1;
M2short=ceil(M/2);

if lt(2)==1
    if phasetype==1
        TimeInd = (0:(N-1))/N;
        FreqInd = (0:(M2-1))*b;
        
        phase = FreqInd'*TimeInd;
        phase = exp(-2*i*pi*phase);
        
        % Handle multisignals
        coef = bsxfun(@times,coef,phase);
        
    end;
    
    f = comp_isepdgtreal(coef,g,L,a,M);
    
else
    % Quinqux lattice
    f=comp_inonsepdgtreal_quinqux(coef,g,a,M);
    
end;

