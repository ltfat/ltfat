%DEMO_DGT_PARAMETRIZE shows how to parametrize the DGT algorithm
%
%   This demo shows two ways of parametrizing the discrete Gabor
%   transform. In the first part, it is shown how to parametrize
%   the dgt to yield a specified signal length for a given window
%   length and overlap.
%   In the second part, the dgt algorithm is used to calculate a
%   classical short-time Fourier transform with `a = b = 1`.
%
%   .. figure::
%
%      Spectrogram of the 'gspi' signal.
%
%
%   See also: dgt, dgtreal, idgt, idgtreal, gabimagepars, demo_dgt

[insig, fs] = gspi;
L = length(insig);
%specify the the parameters
window = 'hanning';
windowlength = 2048;
overlap = 4;

%calculate the minimum hopsize for the given window length
a = floor(windowlength/overlap);

%to exactly maintain the length of the input signal throughout calculations, select M, 
%the number of fft bins (or frequency channels) such that lcm(a,M) divides length(insig)
%by an integer

%select a roughly appropriate starting value
M = 240;
%then successively increase it until length(insig)/lcm(a,M) yields an
%integer result
while length(insig)/lcm(a,M) ~= ceil(length(insig)/lcm(a,M))
    M = M + 1;
end

%specify the dual window used for reconstruction
dualwindow = [window,'dual'];

% use the routines for real valued signals
coef=dgtreal(insig, window, a, M, L);
plotdgtreal(coef, a, M, fs, 90);
outsig = idgtreal(coef,window,a,M);

fprintf('Reconstruction error:      %e\n',norm(insig-outsig(1:length(insig))));

%--------------------------------------------------------------------------
%calculate an STFT

[insig, fs] = greasy;
L = length(insig);

%specify a explicitly and b implicitly via M
a = 1;
M = L;

%the dgtlength is exactly the signal length
L = dgtlength(L,a,M);
window = 'gauss';

%calculate and plot the coefficients
coef = dgt(insig, window, a, M, L);
plotdgt(coef, a, fs, 90);



