function z = comp_fftanalytic(f)
%COMP_FFTANALYTIC Compute analytic representation
%   Usage z = comp_fftanalytic(f, Ls);
%
% `comp_fftanalytic(f)` computes the analytic representation of f.  
% The analytic representation is computed through the FFT of f.
%

% AUTHOR: Jordy van Velthoven

Ls = size(f,1);
H = floor(Ls/2);

z = fft(f);
z(2:Ls-H,:) = 2*z(2:Ls-H,:);
z(H+2:Ls,:) = 0;
z = ifft(z);



