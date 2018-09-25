%DEMO_FRSYNABS  Construction of a signal with a given spectrogram
%
%   This demo demonstrates iterative reconstruction of a spectrogram.
%
%   .. figure::
%
%      Original spectrogram
%
%      This figure shows the target spectrogram
%
%   .. figure::
%
%      Linear reconstruction
%
%      This figure shows a spectrogram of a linear reconstruction of the
%      target spectrogram.
%
%   .. figure::
%
%      Iterative reconstruction using the Griffin-Lim method.
%
%      This figure shows a spectrogram of an iterative reconstruction of the
%      target spectrogram using the Griffin-Lim projection method.
%
%   The BFGS method makes use of the minFunc software. To use the BFGS method, 
%   please install the minFunc software from:
%   `<http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html>`_.
%
%   See also:  isgramreal, isgram

s=ltfattext;

figure(1);
imagesc(s);
colormap(gray);
axis('xy');

figure(2);
F = frame('dgtreal','gauss',8,800); 
scoef = framenative2coef(F,s);
sig_lin = frsyn(F,sqrt(scoef));
sgram(sig_lin,'dynrange',100);

figure(3);
sig_griflim = frsynabs(F,scoef);
sgram(sig_griflim,'dynrange',100);

