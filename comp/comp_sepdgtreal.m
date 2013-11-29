function [coef]=comp_sepdgtreal(f,g,a,M)
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
        %c=reshape(c,M2,N,W);
        
    else
        % Do the factorization algorithm 
        coef=comp_dgtreal_long(f,g,a,M);
        %c=reshape(c,M2,N,W);

    end;