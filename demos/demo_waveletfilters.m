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
%
%   See also: waveletfilters, freqwavelet, filterbank, filterbankrealdual

%the input signal
[f,fs]=gspi;
Ls = length(f);

%wavelet types
wvlt = {{'cauchy', 300}, {'morlet', 3}};

%the frequency range and spacing to be covered
M0 = 511; %desired number of channels (without low frequency compensation channels)
fmax = fs/2; %Nyquist frequency
fmin = 20;

%derive the associated scales
max_freq = 10;  % 10 corresponds to waveltfilters' internal Nyquist frequency
freq_step = max_freq/M0;
min_freqHz = fs/10*freq_step;
start_index = 6; %set the desired number of compensation filters
min_scale_freq = min_freqHz*start_index;
min_freqDiv10 = freq_step*start_index;
scales = 1./linspace(min_freqDiv10,max_freq,M0);

%define the delay function
alpha = 1-2/(1+sqrt(5)); % 1-1/(goldenratio) delay sequence
delays = @(n,a) a*(mod(n*alpha+.5,1)-.5);

%calling waveletfilters by passing the frequency range
[~, ~,fc_freqlin,~,~] = waveletfilters(Ls,'linear', fs, fmin, fmax, M0);

%calling waveletfilters in a constant-Q style
[~, ~,fc_cqt,~,~] = waveletfilters(Ls,'bins', fs, fmin, fmax, 50);

%calling waveletfilters via scales
[g, a,fc_scales,L,info] = waveletfilters(Ls,scales,wvlt{1},'uniform',...
    'repeat','energy', 'delay',delays, 'redtar', 4);

fbbounds = filterbankrealbounds(g,a,L);

if ~isinf(fbbounds)
    gd=filterbankrealdual(g,a,L);
    % Plot frequency responses of individual filters
    figure(1);
    subplot(2,1,1);
    filterbankfreqz(gd,a,L,fs,'plot','linabs','posfreq');

    subplot(2,1,2);
    filterbankfreqz(g,a,L,fs,'plot','linabs','posfreq');
    
    c=filterbank(f,gd,a);
    fcons=2*real(ifilterbank(c,g,a));
    if length(fcons) > length(f)
        err=norm(fcons(1:length(f))-f)
    else
        err=norm(fcons-f(1:length(fcons)))
    end
    fprintf('Reconstruction error:      %e\n',err);
else
    c=filterbank(f,g,a);
    [fpcg,~,iterpcg] = ifilterbankiter(c,g,a,'pcg');
    if length(fpcg) > length(f)
        err=norm(fpcg(1:length(f))-f);
    else
        err=norm(fpcg-f(1:length(fpcg)));
    end
    fprintf('Reconstruction error:      %e\n',err);
end