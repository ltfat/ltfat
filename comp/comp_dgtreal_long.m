function cout=comp_dgtreal_long(f,g,a,M)
%COMP_DGTREAL_LONG  Full-window factorization of a Gabor matrix.
%   Usage:  c=comp_dgtreal_long(f,g,a,M);
%
%   Input parameters:
%         f      : Factored input data
%         g      : Window
%         a      : Length of time shift.
%         M      : Number of channels.
%   Output parameters:
%         c      : M x N*W*R array of coefficients, where N=L/a
%
%   Do not call this function directly, use DGT instead.
%   This function does not check input parameters!
%
%   References: so07-2 st98-8

%   AUTHOR    : Peter L. Søndergaard.
%   TESTING   : TEST_DGT
%   REFERENCE : OK 
  
L=size(f,1);
W=size(f,2);
N=L/a;
M2=floor(M/2)+1;
  
gf=comp_wfac(g,a,M);

% Compute the window application
% We know the output is real, but comp_dgt_walnut cannot detect this, so
% we force the output to be real.
cout=real(comp_dgt_walnut(f,gf,a,M));

% FFT with only positive frequencies
cout=fftreal(cout)/sqrt(M);
cout=reshape(cout,M2,N,W);


