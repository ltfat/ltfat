%DEMO_GABMIXDUAL  How to use GABMIXDUAL
%
%   This script illustrates how one can produce dual windows
%   using GABMIXDUAL
%
%   The demo constructs a dual window that is more concentrated in
%   the time domain by mixing the original Gabor window by one that is
%   extremely well concentrated. The result is somewhat in the middle
%   of these two.
% 
%   The lower framebound of the mixing Gabor system is horrible,
%   but this does not carry over to the gabmixdual.
%
%   .. figure::
%
%      Gabmixdual of two Gaussians.
%
%      The first row of the figure shows the canonical dual window
%      of the input window, which is a Gaussian function perfectly
%      localized in the time and frequency domains.
%
%      The second row shows the canonical dual window of the window we
%      will be mixing with: This is a Gaussian that is 10 times more
%      concentrated in the time domain than in the frequency domain.
%      The resulting canonical dual window has rapid decay in the time domain.
%
%      The last row shows the gabmixdual of these two. This is a non-canonical
%      dual window of the first Gaussian, with decay resembling that of the
%      second.
%
%   See also: gabmixdual, gabdual

disp('Type "help demo_gabmixdual" to see a description of how this demo works.');

L=120;
a=10;
M=12;

% Compute frequency shift.
b=L/M;

% Optimally centered Gaussian
g1=pgauss(L);

% Compute and print framebounds.
[A1,B1]=gabframebounds(g1,a,M);
disp('');
disp('Framebounds of initial Gabor system:');
A1, B1

% Narrow Gaussian
g2=pgauss(L,.1);

% Compute and print framebounds.
[A2,B2]=gabframebounds(g2,a,M);
disp('');
disp('Framebounds of mixing Gabor system:');
A2, B2

% Create a gabmixdual. The window gd is a dual window to g1
gd=gabmixdual(g1,g2,a,M);

% Compute and print framebounds.
[Am,Bm]=gabframebounds(gd,a,M);
disp('');
disp('Framebounds of gabmixdual Gabor system:');
Am, Bm

% Create canonical duals, for plotting.
gc1=gabdual(g1,a,M);
gc2=gabdual(g2,a,M);

% Standard note on plotting:
%
% - The windows are all centered around zero, but this
%   is not visually pleasing, so the window must be
%   shifted to the middle by an FFTSHIFT

figure(1);

subplot(3,2,1);
plot(fftshift(gc1));
title('Canonical dual window.');  

subplot(3,2,2);
plot(20*log10(abs(fftshift(gc1))));
title('Decay of canonical dual window.');  

subplot(3,2,3);
plot(fftshift(gc2));
title('Can. dual of mix. window.');

subplot(3,2,4);
plot(20*log10(abs(fftshift(gc2))));
title('Decay of can.dual of mix. window.')

subplot(3,2,5);
plot(fftshift(gd));
title('Gabmixdual');

subplot(3,2,6);
plot(20*log10(abs(fftshift(gd))));
title('Decay of gabmixdual');
