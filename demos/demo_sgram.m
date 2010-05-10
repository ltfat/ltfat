%DEMO_sgram  Give a simple spectrogram demo
%
%   This script loads the "linus" signal, and uses 'sgram' on it,
%
%   FIGURE 1 color coded image of the spectrogram modulus logarithm
%

f = linus;
figure()
sgram(linus,8000);
