function [f]=comp_isepdgtreal(coef,g,L,a,M)
%COMP_ISEPDGTREAL  Separable IDGT.
%   Usage:  f=comp_isepdgtreal(c,g,L,a,M);
%       
%   This is a computational routine. Do not call it directly.
%
%   Input must be in the M x N x W format, so the N and W dimension is
%   combined.
%
%   See also: idgt

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

Lwindow=size(g,1);


    if L==Lwindow
        % Do full-window algorithm.
        
        % Get the factorization of the window.
        %gf = comp_wfac(g,a,M);      
        
        % Call the computational subroutine.
        f = comp_idgtreal_long(coef,g,L,a,M);
        
    else
        % Do filter bank algorithm.
        % Call the computational subroutine.
        f = comp_idgtreal_fb(coef,g,L,a,M);
    end;