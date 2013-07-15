function [h,g,a,info]=wfilt_lemarie(N)
%WFILT_LEMARIE  Battle and Lemarie filters
%   Usage: [h,g,a]=wfilt_lemarie(N)
%
%   Input parameters:
%         N     : Filter length.
%
%   `[h,g,a]=wfilt_lemarie(N)` calculates coeficients of orthonormal
%   Battle-Lemarie wavelets. Filter coefficients are obtainded by
%   frequency domain sampling and trunctating the impulse response.
%   
%
%   References: mallat89atheory

% Original copyright goes to:
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
% Author: Jose Martin Garcia
% e-mail: Uvi_Wave@tsc.uvigo.es

num_coefs = N;
L = 1024;
H = wfreq_lemarie(L);
hh=real(ifft(H{1},L));
hh=[ hh(L-floor(num_coefs/2)+1:L) hh(1:ceil(num_coefs/2))];
hh=hh/norm(hh);

g{1} = fliplr(hh);
g{2} = -(-1).^(1:length(hh)).*g{1}(end:-1:1);
 
h=g;

a= [2;2];
info.istight = 1;




