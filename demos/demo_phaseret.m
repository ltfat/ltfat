%DEMO_PHASERET Phase retrieval and phase difference
%
%   This demo demonstrates iterative reconstruction of a spectrogram and
%   the phase difference.
%
%   .. figure::
%
%      Original spectrogram
%
%      This figure shows the target spectrogram of an excerpt of the gspi
%      signal
%
%   .. figure::
%
%      Phase difference
%
%      This figure shows a difference between the original phase and the 
%      reconstructed using 100 iterations of a Fast Griffin-Lim algorithm.
%      Note: The figure in the LTFAT 2.0 paper differs slightly because it
%      is genarated using 1000 iterations.
%
%   See also:  frsynabs, plotframe

% Here the number of iterations is set to 100 to speedup the execution
maxit = 100;

[f,fs]=gspi;
f=f(10001:100000);
g='gauss';
a=100; M=1000;
F = frame('dgtreal',g,a,M);
c = frana(F,f);

% Spectrogram values
s=abs(c).^2;
% Original phase
theta=angle(c);

% Do the reconstruction using the magnitude only
r = frsynabs(F,sqrt(s),'fgriflim','maxit',maxit);

% Re-analyse using the original frame
c_r = frana(F,r);

% Obtain phase of the re-analysis
s_r = abs(c_r).^2;

theta_r=angle(c_r);

d1=abs(theta-theta_r);
d2=2*pi-d1;
anglediff=min(d1,d2);

% Plot the original spectrogram
figure(1);
plotframe(F,c,fs,'dynrange',50);

% Plot the phase difference
figure(2);
anglediff(abs(c)<10^(-50/20)) = 0;
plotframe(F,anglediff,fs,'lin');


