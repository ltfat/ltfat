function [H,G] = wfreq_lemarie(L)
%WFREQ_LEMARIE  Battle and Lemarie filters frequency resp. sampling
%   Usage: [H,G]=wfreq_lemarie(L)
%
%   Input parameters:
%         N     : Number of samples of the frequency response.
%
%   `[H,G]=wfreq_lemaire(L)` calculates $L$ samples of the Battle and
%   Lemarie filters frequency responses.
%
%   References: mallat89atheory
%
%

% Original copyright goes to:
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
% Author: Jose Martin Garcia
% e-mail: Uvi_Wave@tsc.uvigo.es


% frequency axis
w=[0:2*pi/L:2*pi*(1-1/L)];
w(1)=eps;
w(L/2+1)=w(L/2+1)+1e-15;

% calculation of frequency response of analysis lowpass filter 
num=0;den=0;
K=36;
for k=-K:K,
	num=1./((w+2*pi*k).^8)+num;
	den=1./((2*w+2*pi*k).^8)+den;
end
H = cell(2,1);
H{1}=sqrt(num./(2.^8*den));
H{1}(1)=1;

H{2} = fftshift(H{1});
G = cell(2,1);
G{1} = fliplr(H{1});
G{2} = fliplr(H{2});



