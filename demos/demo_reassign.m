%DEMO_REASSIGN  Give demos of Gabor reassignment
%
%   This script loads a sample signal and computes.


% Load bat sonar signal
batsignal = bat;
siglength = length(batsignal);
fs = 143;

% Compute Gabor transform and phase gradients, and reassigned Gabor
% transform
a = 1;
M = 200;
Mover2 = M/2;
[tgrad, fgrad, batgt] = gabphasegrad('dgt',batsignal,'gauss',a,M);
batreass = gabreassign(abs(batgt).^2,tgrad,fgrad,1);

% Displays
ttt = (1:siglength)*1000/fs;
fff = (0:(Mover2-1))*fs/M;

figure('name','Gabor transform modulus');
imagesc(ttt,fff,abs(batgt(1:100,:))); axis xy;
xlabel('Time (msec.)');
ylabel('Frequency (Hz)');
title('Gabor transform modulus');

figure('name','Reassigned Gabor transform modulus');
imagesc(ttt,fff,abs(batreass(1:100,:))); axis xy;
xlabel('Time (msec.)');
ylabel('Frequency (Hz)');
title('Reassigned Gabor transform modulus');


%OLDFORMAT
