function [coef]=comp_sepdgtreal(f,g,a,M,phasetype)
%COMP_SEPDGTREAL  Filter bank DGT
%   Usage:  c=comp_sepdgtreal(f,g,a,M);
%  
%   This is a computational routine. Do not call it directly.

%   See help on DGT.

%   AUTHOR : Peter L. Søndergaard.

L=size(f,1);
Lwindow=size(g,1);

if Lwindow<L
    % Do the filter bank algorithm
    % Periodic boundary conditions
    coef=comp_dgtreal_fb(f,g,a,M);

else
    % Do the factorization algorithm 
    coef=comp_dgtreal_long(f,g,a,M);
end;

% Change the phase convention from frequency-invariant to
% time-invariant
if phasetype==1
    N=L/a;
    M2=floor(M/2)+1;
    
    TimeInd = (0:(N-1))*a;
    FreqInd = (0:(M2-1))/M;
         
    phase = FreqInd'*TimeInd;
    phase = exp(2*1i*pi*phase);
    coef=bsxfun(@times,coef,phase);
end;

    
    
