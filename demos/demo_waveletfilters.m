%DEMO_WAVELETFILTERS  Introduction to grid-like wavelet sampling
%
%   This demo shows three options to generate a wavelet filter bank.
%   Firstly, a filterbank with linearly spaced center frequencies is generated, 
%   by passing the desired frequency range directly. The filter bank
%   coefficients are visualized.
%   Secondly, waveletfilters is conventionally parametrized by direct
%   passing of the desired scales and the frequency response of the
%   filter bank and its dual are plotted.
%   Finally, a constant-Q style wavelet filter bank is generated.
%
%   .. figure::
%
%      Coefficients of an invertible wavelet filterbank with linearly 
%      spaced center frequencies.
%
%   .. figure::
%
%      Frequency response of the filterbank and its dual.
%
%      The figure shows the filter bank response with linear scales.
%
%
%   See also: waveletfilters, freqwavelet, lowdiscrepancy, filterbankrealdual

%the input signal
[f,~]=gspi;
Ls = length(f);

%wavelet types
wvlt = {{'cauchy', 600}, {'morlet', 3}};

cauchyalpha = [300, 700, 1200, 1700, 2100];

%design an invertible wavelet filter bank with linearly spaced center
%frequencies
fs = 2;
M = 64;
MC = 3;
fmin = MC/M;
fmax = fs/2; %up to Nyquist frequency
numscales = M-MC+1;
redundancy = 4;
delays = lowdiscrepancy('digital');

%calling waveletfilters by passing the frequency range
[g,a,fc, L]=waveletfilters(Ls,'linear',fs,fmin,fmax,numscales,...
    {'cauchy',cauchyalpha(1)},'uniform','redtar',redundancy,'repeat','delay',delays);

c = filterbank(f, g, a);
    
figure(1)
plotfilterbank(c, a)


%call waveletfilters by passing the scales directly
M = 127; %desired number of channels (without low frequency compensation channels)

%derive the associated scales
max_freq = 10;  % 10 corresponds to waveltfilters' internal Nyquist frequency
min_freq = max_freq/M;
scales = 1./linspace(min_freq,max_freq,M);

[g, a,fc_scales,L,info] = waveletfilters(Ls,scales,wvlt{2},'uniform',...
    'single','energy', 'delay',delays, 'redtar', 8);

gd=filterbankrealdual(g,a,L);
% Plot frequency responses of individual filters
figure(2);
subplot(2,1,1);
filterbankfreqz(gd,a,L,fs,'plot','linabs','posfreq');

subplot(2,1,2);
filterbankfreqz(g,a,L,fs,'plot','linabs','posfreq');
    

%calling waveletfilters in a constant-Q style
fmin = 20;
fmax = 8000;
bins = 16;
fs = 16000;
[g, a,fc_cqt, L] = waveletfilters(Ls,'bins', fs, fmin, fmax, bins, 'delay', lowdiscrepancy('kronecker'));


c=filterbank(f,g,a);
[fpcg,~,iterpcg] = ifilterbankiter(c,g,a,'pcg');
if length(fpcg) > length(f)
    err=norm(fpcg(1:length(f))-f);
else
    err=norm(fpcg-f(1:length(fpcg)));
end
fprintf('Reconstruction error:      %e\n',err);
