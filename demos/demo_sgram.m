%DEMO_sgram  Give a simple spectrogram demo
%
%   This script displays two spectrogram of the "greasy" test signal,
%   demonstrating some common switches.
%
%   FIGURE 1 Spectrogram using default window
%
%     The figure shows a spectrogram of the 'greasy' test signal. The
%     magnitude of the Gagor coefficients is shown on a logarithmic
%     scale, and only the largest coefficients are shows, corresponding
%     to a dynamic range of 50 dB.
%
%   FIGURE 2 Spectrogram using longer window
%
%     Same spectrogram as Figure 1, but now using a longer window, with a
%     window length of 20 ms. This gives a sharper with of the partials,
%     at the expense of a more blurred wiew along the time dimension.
%
%   See also: sgram, greasy

% Test signal
f=greasy;

% Sampling frequency of the signal
fs=16000;

% Default window
figure(1)
sgram(f,fs,'dynrange',50);

% Longer window, 20 ms.
figure(2)
sgram(f,fs,'wlen',round(20/1000*fs),'dynrange',50);
