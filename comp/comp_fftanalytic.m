function z = comp_fftanalytic(s)
%COMP_FFTANALYTIC Compute analytic representation
%
%   Usage: z = comp_fftanalytic(s);
%
%   `comp_fftanalytic(s)` computes the analytic representation of s.  
%   The analytic representation is computed through the FFT of f.
%

% AUTHOR: Jordy van Velthoven

Ls = size(s,1);
H = floor(Ls/2);

z = fft(s);
z(2:Ls-H,:) = 2*z(2:Ls-H,:);
z(H+2:Ls,:) = 0;
z = ifft(z);



