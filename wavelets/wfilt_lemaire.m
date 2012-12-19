function [h,g,a]=wfilt_lemaire(N)
%WFILT_LEMARIE  Battle and Lemarie filters.
%   Usage: [h,g,a]=wfilt_lemaire(N)
%
%   Input parameters:
%         N     : Filter length.
%
%   `[h,g,a]=wfilt_lemaire(N)` calculates coeficients of orthonormal
%   Battle-Lemarie wavelets. Filter coefficients are obtainded by
%   frequency domain sampling and trunctating the impulse response.
%
% References: mallat89atheory
%
% Original copyright goes to:
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
% Author: Jose Martin Garcia
% e-mail: Uvi_Wave@tsc.uvigo.es

num_coefs = N;
L = 1024;
H = wfreq_lemaire(L);
hh=real(ifft(H{1},L));
hh=[ hh(L-floor(num_coefs/2)+1:L) hh(1:ceil(num_coefs/2))];
hh=sqrt(2)/sum(hh)*hh;

g{1} = fliplr(hh);
g{2} = -(-1).^(1:length(hh)).*g{1}(end:-1:1);
 
h{1}=g{1}(length(g{1}):-1:1);
h{2}=g{2}(length(g{2}):-1:1);

a= [2;2];




