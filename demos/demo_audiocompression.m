%DEMO_AUDIOCOMPRESSION  Audio compression using N-term approx
%
%   This demos shows how to do audio compression using best N-term
%   approximation of an |wmdct| transform.
%
%   The signal is transformed using an orthonormal |wmdct| transform.
%   Then approximations with a fixed number *N* of coefficients are obtained
%   by:
%
%     * Linear approximation: The *N* coefficients with lowest frequency
%       index are kept.
%
%     * Non-linear approximation: The *N* largest coefficients (in
%       magnitude) are kept.
%
%   The corresponding approximated signal can be computed using |iwmdct|.
%
%   .. figure::
%
%      Rate-distorition plot
%
%      The figure shows the output Signal to Noise Ratio (SNR) as a function
%      of the number of retained coefficients.
%
%   Note: The inverse WMDCT is not needed for computing computing
%   SNRs. Instead Parseval theorem states that the norm of a signal equals
%   the norm of the sequence of its |wmdct| coefficients.

% Load audio signal
% Use the 'glockenspiel' signal.
sig=gspi;

% Shorten signal
L = 2^16;
sig = sig(1:L);

% Number of frequency channels
M = 1024;

% Number of time steps
N = L/M;

% Generate window
gamma = wilorth(M,L);

% Compute wmdct coefficients
c = wmdct(sig,gamma,M);


% L2 norm of signal
InputL2Norm = norm(c,'fro');

% Approximate, and compute SNR values
kmax = M;
kmin = kmax/32;             % 32 is an arbitrary choice
krange = kmin:32:(kmax-1);  % same remark

for k = krange,
    ResL2Norm_NL = norm(c-largestn(c,k*N),'fro');
    SNR_NL(k) = 20*log10(InputL2Norm/ResL2Norm_NL);
    ResL2Norm_L = norm(c(k:kmax,:),'fro');
    SNR_L(k) = 20*log10(InputL2Norm/ResL2Norm_L);
end


% Plot
figure(1);

set(gca,'fontsize',14);
plot(krange*N,SNR_NL(krange),'x-b',...
     krange*N,SNR_L(krange),'o-r');
axis tight; grid;
legend('Best N-term','Linear');
xlabel('Number of Samples', 'fontsize',14);
ylabel('SNR (dB)','fontsize',14);
