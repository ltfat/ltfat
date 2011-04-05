%DEMO_ISGRAM  Contruction of a signal with a given spectrogram
%
%   This demo demonstrates iterative reconstruction of a spectrogram.
%
%   FIGURE 1 Original spectrogram
%
%      This figure shows the target spectrogram
%
%   FIGURE 2 Linear reconstruction
%
%      This figure shows a spectrogram of a linear reconstruction of the
%      target spectrogram.
%
%   FIGURE 3 Iterative reconstruction using the Griffin-Lim method.
%
%      This figure shows a spectrogram of an iterative reconstruction of the
%      target spectrogram using the Griffin-Lim projection method.
%
%   FIGURE 4 Iterative reconstruction using the BFGS method
%
%      This figure shows a spectrogram of an iterative reconstruction of the
%      target spectrogram using the BFGS method.

%   See also:  isgramreal, isgram, 

s=ltfattext;

figure(1);
imagesc(s);
colormap(gray);
axis('xy');

figure(2);
sig_lin = idgtreal(s,'gauss',8,800);
sgram(sig_lin,'dynrange',100);

figure(3);
sig_iter = isgramreal(s,'gauss',8,800);
sgram(sig_iter,'dynrange',100);

figure(4);
sig_iter = isgramreal(s,'gauss',8,800,'bfgs');
sgram(sig_iter,'dynrange',100);