%DEMO_AUDIOSHRINK  Decomposition into tonal and transient parts
%
%   This demos shows how to do audio coding and "tonal + transient"
%   decomposition using group lasso shrinkage of two |wmdct| transforms
%   with different time-frequency resolutions.
%
%   The signal is transformed using two orthonormal |wmdct| bases.
%   Then group lasso shrinkage is applied to the two transforms
%   in order to:
%
%     * select fixed frequency lines of large |wmdct| coefficients on the
%       wide window |wmdct| transform
%
%     * select fixed time lines of large |wmdct| coefficients on the
%       narrow window |wmdct| transform
% 
%   The corresponding approximated signals are computed with the
%   corresponding inverse, |iwmdct|.
%
%   .. figure:: 
%
%      Plots and time-frequency images
%
%      The upper plots in the figure show the tonal parts of the signal, the
%      lower plots show the transients. The TF-plots on the left are the
%      corresponding wmdct coefficients found by appropriate group lasso
%      shrinkage
%
%   Corresponding reconstructed tonal and transient sounds may be
%   listened from arrays rec1 and rec2 (sampling rate: 44.1 kHz)
%


% Load audio signal and add noise
% -------------------------------
% Use the 'glockenspiel' signal.
sig=gspi;
fs=44100;

% Shorten signal
siglen = 2^16;
sig = sig(1:siglen);

% Add Gaussian white noise
nsig = sig + 0.01*randn(size(sig));

% Tonal layer
% -----------

% Create a WMDCT basis with 256 channels
F1=frametight(frame('wmdct','gauss',256));

% Group lasso and invert
c1 = franagrouplasso(F1,nsig,0.8,'soft','freq');
rec1 = frsyn(F1,c1);

% Transient layer
% ---------------

% Create a WMDCT basis with 32 channels
F2=frametight(frame('wmdct','gauss',32));

c2 = franagrouplasso(F2,nsig,0.5,'soft','time');
rec2 = frsyn(F2,c2);

% Plots
% -----

% Dynamic range for plotting
dr=50;
xplot=(0:siglen-1)/fs;

figure(1);
subplot(2,2,1);
plot(xplot,rec1);
xlabel('Time (s)');
axis tight;

subplot(2,2,2);
plotframe(F1,c1,fs,dr);

subplot(2,2,3);
plot(xplot,rec2);
xlabel('Time (s)');
axis tight;

subplot(2,2,4);
plotframe(F2,c2,fs,dr);

% Count the number of non-zero coefficients
N1=sum(abs(c1)>0);
N2=sum(abs(c2)>0);

p1 = 100*N1/siglen;
p2 = 100*N2/siglen;
p=p1+p2;

fprintf('Percentage of retained coefficients: %f + %f = %f\n',p1,p2,p);

disp('To play the original, type "soundsc(sig,fs)"');
disp('To play the tonal part, type "soundsc(rec1,fs)"');
disp('To play the transient part, type "soundsc(rec2,fs)"');
