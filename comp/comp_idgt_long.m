function f=comp_idgt_long(coef,g,L,a,M)
%COMP_IDGT_FAC  Full-window factorization of a Gabor matrix.
%   Usage:  f=comp_idgt_long(c,g,L,a,M)
%
%   Input parameters:
%         c     : M x N x W array of coefficients.
%         g     : Window.
%         a     : Length of time shift.
%         M     : Number of frequency shifts.
%   Output parameters:
%         f     : Reconstructed signal.
%
%   Do not call this function directly, use IDGT.
%   This function does not check input parameters!
%
%   If input is a matrix, the transformation is applied to
%   each column.
%
%   This function does not handle multidimensional data, take care before
%   you call it.
%
%   References: so07-2 st98-8

%   AUTHOR : Peter L. Søndergaard.

% Get the factorization of the window.
gf = comp_wfac(g,a,M);      

% Call the computational subroutine.
f  = comp_idgt_fac(coef,gf,L,a,M);
