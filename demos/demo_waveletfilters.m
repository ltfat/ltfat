%DEMO_WAVELETFILTERS  Introduction to grid-like wavelet sampling
%
%   This demo shows how to generate invertible wavelet filter banks with linear
%   frequency spacing.
%   Firstly, a conventional wavelet filter bank is generated.
%   In a second iteration, the application of a small delay allows for the
%   design of a filter bank with a linearly spaced frequency axis that
%   yields a similar reconstruction error as the first filter bank. The
%   reconstruction can, however, be further improved by adding low
%   frequency compensation filters as shown for the third filter bank.
%   Alternative choices of the delay generating function are possible.
%
%   .. figure::
%
%      Coefficients of the three filter banks.
%
%
%   See also: waveletfilters, freqwavelet, lowdiscrepancy, filterbankrealdual

%the input signal
[f,fs]=greasy;
Ls = length(f);

%retrieve the filters and downsampling factors for a wavelet filterbank
%with the following parameters:
fmin = 16;
fmax = 8000;
bins = 16;
[g_cq, a_cq, fc_cq] = waveletfilters(Ls,'bins', fs, fmin, fmax, bins);

%apply the filter bank thus generated
c_cq=filterbank(f,g_cq,a_cq);

%perform iterative reconstruction and calculate the error
fpcg = ifilterbankiter(c_cq,g_cq,a_cq,'pcg');
if length(fpcg) > length(f)
    err=norm(fpcg(1:length(f))-f);
else
    err=norm(fpcg-f(1:length(fpcg)));
end
fprintf('Reconstruction error:      %e\n',err);


%now, design an invertible wavelet filter bank covering the same frequency
%range with linearly spaced center frequencies

%use the same number of frequency channels as before
M = numel(fc_cq);
%delays is an anonymous function that corresponds to the desired low 
%discrepancy sequence; it can be passed directly to waveletfilters
delays = lowdiscrepancy('kronecker');

%to each sampled wavelet, a small delay 'delays' is applied to achieve even
%coverage of the time-frequency plane;
%a single compensation filter is added per default
[g_del, a_del, fc_del]=waveletfilters(Ls,'linear',fs,fmin,fmax,M, 'delay',delays);

c_del = filterbank(f, g_del, a_del);

%perform iterative reconstruction and calculate the error
fpcg = ifilterbankiter(c_del,g_del,a_del,'pcg');
if length(fpcg) > length(f)
    err=norm(fpcg(1:length(f))-f);
else
    err=norm(fpcg-f(1:length(fpcg)));
end
fprintf('Reconstruction error (with delay):      %e\n',err);

%now, further improve reconstruction by adding 3 compensation filters in the
%lower frequency range
MC = 3;
fmin = MC/M * fs;
numscales = M-MC+1;

[g_comp,a_comp]=waveletfilters(Ls,'linear',fs,fmin,fmax,numscales,'repeat','delay',delays);

c_comp = filterbank(f, g_comp, a_comp);

%perform iterative reconstruction and calculate the error
fpcg = ifilterbankiter(c_comp,g_comp,a_comp,'pcg');
if length(fpcg) > length(f)
    err=norm(fpcg(1:length(f))-f);
else
    err=norm(fpcg-f(1:length(fpcg)));
end
fprintf('Reconstruction error (kronecker sequence delay and low-f compensation):      %e\n',err);

figure
subplot(1,3,1)
plotfilterbank(c_cq, a_cq)
xlim([0 Ls])
subplot(1,3,2)
plotfilterbank(c_del, a_del)
xlim([0 Ls])
subplot(1,3,3)
plotfilterbank(c_comp, a_comp)
xlim([0 Ls])


delays = lowdiscrepancy('digital');
[g,a]=waveletfilters(Ls,'linear',fs,fmin,fmax,numscales,'repeat','delay',delays);

c = filterbank(f, g, a);

%perform iterative reconstruction and calculate the error
fpcg = ifilterbankiter(c,g,a,'pcg');
if length(fpcg) > length(f)
    err=norm(fpcg(1:length(f))-f);
else
    err=norm(fpcg-f(1:length(fpcg)));
end
fprintf('Reconstruction error (digital net delay and low-f compensation):      %e\n',err);


% M=64;
% MC = 3;
% fmin = MC/M;
% fmax = fs/2; %up to Nyquist frequency
% numscales = M-MC+1;
% redundancy = 4;
% delays = lowdiscrepancy('digital');
% 
% %calling waveletfilters by passing the frequency range
% [g,a,fc, L]=waveletfilters(Ls,'linear',fs,fmin,fmax,numscales,...
%     {'cauchy',cauchyalpha(1)},'uniform','redtar',redundancy,'repeat','delay',delays);
% 
% c = filterbank(f, g, a);
%     
% figure(1)
% plotfilterbank(c, a)
% 
% 
% %call waveletfilters by passing the scales directly
% M = 127; %desired number of channels (without low frequency compensation channels)
% 
% %derive the associated scales
% max_freq = 10;  % 10 corresponds to waveltfilters' internal Nyquist frequency
% min_freq = max_freq/M;
% scales = 1./linspace(min_freq,max_freq,M);
% 
% [g, a,fc_scales,L,info] = waveletfilters(Ls,scales,wvlt{2},'uniform',...
%     'single','energy', 'delay',delays, 'redtar', 8);
% 
% gd=filterbankrealdual(g,a,L);
% % Plot frequency responses of individual filters
% figure(2);
% subplot(2,1,1);
% filterbankfreqz(gd,a,L,fs,'plot','linabs','posfreq');
% 
% subplot(2,1,2);
% filterbankfreqz(g,a,L,fs,'plot','linabs','posfreq');
%     
% 
% %calling waveletfilters in a constant-Q style
% fmin = 20;
% fmax = 8000;
% bins = 16;
% fs = 16000;
% [g, a,fc_cqt, L] = waveletfilters(Ls,'bins', fs, fmin, fmax, bins);
% 
% 
% c=filterbank(f,g,a);
% [fpcg,~,iterpcg] = ifilterbankiter(c,g,a,'pcg');
% if length(fpcg) > length(f)
%     err=norm(fpcg(1:length(f))-f);
% else
%     err=norm(fpcg-f(1:length(fpcg)));
% end
% fprintf('Reconstruction error:      %e\n',err);
