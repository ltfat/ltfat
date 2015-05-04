function [f]=comp_isepdgtreal(coef,g,L,a,M,phasetype)
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

if phasetype==1
    % Change from time-invariant phase convention to a
    % frequency-invariant one
    b=L/M;
    M2=floor(M/2)+1;
    N=size(coef,2);
    %M2short=ceil(M/2);

    TimeInd = (0:(N-1))/N;
    FreqInd = (0:(M2-1))*b;

    phase = FreqInd'*TimeInd;
    phase = exp(-2*1i*pi*phase);

    % Handle multisignals
    coef = bsxfun(@times,coef,phase);
end;


if L==Lwindow
    % Do full-window algorithm.

    % Call the computational subroutine.
    f = comp_idgtreal_long(coef,g,L,a,M);

else
    % Do filter bank algorithm.
    % Call the computational subroutine.
    f = comp_idgtreal_fb(coef,g,L,a,M);
end;
