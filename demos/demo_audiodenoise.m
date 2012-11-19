%DEMO_AUDIODENOISE  Audio denoising using thresholding
%
%   This demos shows how to do audio denoising using thresholding
%   of WMDCT transform.
%
%   The signal is transformed using an orthonormal WMDCT transform
%   followed by a thresholding. Then the signal is reconstructed
%   and compared with the original.
%
%   .. figure::
%
%      Denoising
%
%      The figure shows the original signal, the noisy signal and denoised
%      signals using hard and soft threshholding applied to the WMDCT of the
%      noise signal.
%

% Load audio signal
% Use the 'glockenspiel' signal.
sig=gspi;

SigLength = 2^16;
sig = sig(1:SigLength);

% Initializations
NbFreqBands = 1024;
sigma = 0.5;
Relative_Threshold = 0.1;
tau = Relative_Threshold*sigma;

% Generate window
gamma = wilorth(NbFreqBands,SigLength);

% add noise to the signal
sigma = sigma * std(sig);
nsig = sig + sigma * randn(size(sig));

% Compute wmdct coefficients
c = wmdct(nsig,gamma,NbFreqBands);

% Hard Thresholding
chard=thresh(c,tau);

% Reconstruct
hrec = real(iwmdct(chard,gamma));

% Soft thresholding
csoft=thresh(c,tau,'soft');

% Reconstruct
srec = real(iwmdct(csoft,gamma));

% Plot
figure(1);
subplot(4,1,1); plot(sig); legend('Original');
subplot(4,1,2); plot(nsig); legend('Noisy');
subplot(4,1,3); plot(hrec); legend('Hard threshold');
subplot(4,1,4); plot(srec); legend('Soft threshold');

% Results
InputSNR = 20 *log10(std(sig)/std(nsig-sig));
OutputSNR_h = 20 *log10(std(sig)/std(hrec-sig));
OutputSNR_s = 20 *log10(std(sig)/std(srec-sig));

fprintf(' RESULTS:\n');
fprintf('      Input SNR: %f dB.\n',InputSNR);
fprintf('      Output SNR (hard): %f dB.\n',OutputSNR_h);
fprintf('      Output SNR (soft): %f dB.\n',OutputSNR_s);
fprintf(' Signals are stored in variables sig, nsig, hrec, srec\n');
