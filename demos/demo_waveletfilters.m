%DEMO_WAVELETFILTERS  Introduction to grid-like wavelet sampling
%
%   This demo shows how to perform grid-like wavelet sampling.
%
%   .. figure::
%
%      Frequency response of the filterbank and its dual.
%
%      The figure shows the filter bank response with linear scales.
%
%   .. figure::
%
%      3 types of frequency spacing achievable with waveletfilters.
%
%      The figure shows the three types of frequency spacing that can be
%      achieved with waveletfilters.
%
%   See also: waveletfilters, freqwavelet, filterbank, filterbankrealdual

%the input signal
[f,fs]=gspi;
Ls = length(f);

%filter bank settings
M0 = 511; %Desired number of channels (without 0 Hz-lowpass channel)
max_freq = 10;  % 10 corresponds to waveltfilters' internal Nyquist frequency
freq_step = max_freq/M0;
min_freqHz = fs/10*freq_step;

start_index = 1;
min_scale_freq = min_freqHz*start_index;
min_freqDiv10 = freq_step*start_index; %1/25; % By default, the reference scale for freqwavelet has center frequency 0.1
scales = 1./linspace(min_freqDiv10,max_freq,M0);
alpha = 1-2/(1+sqrt(5)); % 1-1/(goldenratio) delay sequence
delays = @(n,a) a*(mod(n*alpha+.5,1)-.5);
CauchyAlpha = 600;
[g, a,fc,L,info] = waveletfilters(Ls,scales,{'cauchy',CauchyAlpha},'uniform','single','energy', 'delay',delays, 'redtar', 8);


c=filterbank(f,{'realdual',g},a);
r=2*real(ifilterbank(c,g,a));
if length(r) > length(f)
    norm(r(1:length(f))-f)
else
    norm(r-f(1:length(r)))
end
 
% Plot frequency responses of individual filters
gd=filterbankrealdual(g,a,L);
figure(1);
subplot(2,1,1);
filterbankfreqz(gd,a,L,fs,'plot','linabs','posfreq');

subplot(2,1,2);
filterbankfreqz(g,a,L,fs,'plot','linabs','posfreq');