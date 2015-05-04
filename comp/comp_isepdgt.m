function [f]=comp_isepdgt(coef,g,L,a,M,phasetype)
%COMP_ISEPDGT  Separable IDGT.
%   Usage:  f=comp_isepdgt(c,g,L,a,M);
%       
%   This is a computational routine. Do not call it directly.
%
%   Input must be in the M x N x W format.
%
%   See also: idgt

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

Lwindow=size(g,1);

% FIXME : Calls non-comp function 
if phasetype==1
    coef=phaseunlock(coef,a);
end;

if L==Lwindow
    % Do full-window algorithm.
    % coef=reshape(coef,M,prod(size(coef))/M);

    % Get the factorization of the window.
    %gf = comp_wfac(g,a,M);      

    % Call the computational subroutine.
    f  = comp_idgt_long(coef,g,L,a,M);
else
    %coef=reshape(coef,M,prod(size(coef))/M);
    % Do filter bank algorithm.
    % Call the computational subroutine.

    f=comp_idgt_fb(coef,g,L,a,M);
end;
