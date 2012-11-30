function cout=comp_dgt_long(f,g,a,M)
%COMP_DGT_LONG  Gabor transform using long windows.
%   Usage:  c=comp_dgt_long(f,g,a,M);
%
%   Input parameters:
%         f      : Factored input data
%         g      : Window.
%         a      : Length of time shift.
%         M      : Number of channels.
%   Output parameters:
%         c      : M x N*W*R array of coefficients, where N=L/a
%
%   Do not call this function directly, use DGT instead.
%   This function does not check input parameters!
%
%   References: so07-2 st98-8

%   AUTHOR : Peter L. Søndergaard.

gf=comp_wfac(g,a,M);
  
% Compute the window application
cout=comp_dgt_walnut(f,gf,a,M);

% Apply dft modulation.
cout=fft(cout)/sqrt(M);

% Change to the right shape
L=size(f,1);
W=size(f,2);
N=L/a;
cout=reshape(cout,M,N,W);


