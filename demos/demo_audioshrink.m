%DEMO_AUDIOSHRINK
%
%   This demos shows how to do audio coding and "tonal + transient"
%   decomposition using group lasso shrinkage of two WMDCT transforms
%   with different time-frequency resolutions.
%
%   The signal is transformed using two orthonormal WMDCT bases.
%   Then group lasso shrinkage is applied to the two transforms
%   in order to:
%
%     * select fixed frequency lines of large wmdct coefficients on the
%       wide window wmdct transform
%
%     * select fixed time lines of large wmdct coefficients on the
%       narrow window wmdct transform
% 
%   The corresponding approximated signals are computed with corresponding
%   inverse WMDCT.
%
%   FIGURE 1 plots and time-frequency images
%
%     The upper plots in the figure show the tonal parts of the signal, the
%     lower plots show the transients. The TF-plots on the left are the
%     corresponding wmdct coefficients found by appropriate group lasso
%     shrinkage
%
%   Corresponding reconstructed tonal and transient sounds may be
%   listened from arrays rec1 and rec2 (sampling rate: 44.1 kHz)
%


% Load audio signal and add noise
% -------------------------------
% Use the 'glockenspiel' signal.
sig=gspi;

% Shorten signal
SigLength = 2^16;
sig = sig(1:SigLength);

% Add Gaussian white noise
nsig = sig + 0.01*randn(size(sig));

% Tonal layer
% -----------
% Initializations
M = 256;
N = SigLength/M;

% Generate window
gamma = wilorth(M,SigLength);

% Compute wmdct coefficients
c = wmdct(nsig,gamma,M);

% Group lasso and invert
[c1,N1] = wilgrouplasso(c,0.8,'soft','freq');
rec1 = iwmdct(c1,gamma);

% Transient layer
% ---------------
% Initializations
M = 32;
N = SigLength/M;

% Generate window
gamma = wilorth(M,SigLength);

% Compute wmdct coefficients
c = wmdct(nsig,gamma,M);

[c2,N2] = wilgrouplasso(c,0.5,'time','soft');
rec2 = iwmdct(c2,gamma);

% Plots
% -----

figure(1);
subplot(2,2,1);
plot(rec1);
axis tight;

subplot(2,2,2);
imagesc(log(abs(c1)+0.00001));
set(gca,'ydir','normal');

subplot(2,2,3);
plot(rec2);
axis tight;

subplot(2,2,4);
imagesc(log(abs(c2)+0.00001));
set(gca,'ydir','normal');

p1 = 100*N1/SigLength;
p2 = 100*N2/SigLength;
p=p1+p2;

fprintf('Percentage of retained coefficients: %f + %f = %f\n',p1,p2,p);

if ispc && ~isoctave
    disp('Playing sounds: Original');
    wavplay(sig);
    disp('Playing sounds: Tonal part');
    wavplay(rec1);
    disp('Playing sounds: Transient part');
    wavplay(rec2);
end

