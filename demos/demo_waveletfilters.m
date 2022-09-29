%DEMO_WAVELETFILTERS  Introduction to grid-like wavelet sampling
%
%   This demo shows how to generate invertible wavelet filter banks with linear
%   frequency spacing.
%
%   The following wavelet filter banks are produced:
%
%   * a conventional wavelet filter bank with logarithmic frequency spacing
%   * a wavelet filter bank with linear frequency spacing that achieves
%   similar reconstruction quality
%   * a wavelet filter bank that uses more than one compensation filter to
%   improve the time-frequency coverage of the transform in the
%   low-frequency range
%   * two wavelet filter banks using different types of wavelets with 
%   different Q-factors
%
%   .. figure::
%
%      Coefficients, frequency response, and filter bank response of a
%      conventional wavelet filter bank.
%
%   .. figure::
%
%      Comparison of the spacing of the center frequencies.
%
%   .. figure::
%
%      Comparison of the low-frequency range using a single and several 
%      compensation filters. Energy tends to be more evenly distributed 
%      when more than one compensation filter is used.
%
%   .. figure::
%
%      A Cauchy and a B-spline wavelet with similar Q-factors.
%
%   .. figure::
%
%      Frequency response of the Cauchy and B-spline wavelet filter banks.
%
%   .. figure::
%
%      Coefficients of the Cauchy and B-spline wavelet filter banks.
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
[g_cq, a_cq, fc_cq] = waveletfilters(Ls,'bins', fs, fmin, fmax, bins, 'uniform');

%calculate the dual windows and apply the filterbank
gd = filterbankrealdual(g_cq, a_cq, Ls);
c_cq=filterbank(f,gd,a_cq);

%invert the filterbank
frec = ifilterbank(c_cq, g_cq, a_cq);
%calculate the reconstruction error
if length(frec) > length(f)
    err=norm(frec(1:length(f))-f);
else
    err=norm(frec-f(1:length(frec)));
end
fprintf('Reconstruction error (log. f-spacing):      %e\n',err);

%plot the filter bank coefficients, the frequency response, and the
%filter bank response of the wavelet filter bank
figure(1)
subplot(3,1,1)
plotfilterbank(c_cq, a_cq)
subplot(3,1,2)
filterbankfreqz(g_cq,a_cq,Ls, 'plot', 'posfreq', 'dynrange', 60);
subplot(3,1,3)
filterbankresponse(g_cq,a_cq,Ls, 'plot', 'real', 'fs', fs);

%now, design an invertible wavelet filter bank covering the same frequency
%range with linearly spaced center frequencies

%use the same number of frequency channels as before
M = numel(fc_cq);
%delays is an anonymous function specifying the desired low 
%discrepancy sequence; it can be passed directly to waveletfilters
delays = lowdiscrepancy('kronecker');

%to each sampled wavelet, a small delay 'delays' is applied to achieve even
%coverage of the time-frequency plane;
%a single compensation filter is added per default
[g_del, a_del, fc_del]=waveletfilters(Ls,'linear',fs,fmin,fmax,M,...
    'delay',delays, 'uniform');

c_del = filterbank(f, g_del, a_del);

%invert the filterbank using a short form for calculating the dual window
frec = ifilterbank(c_del, {'realdual',g_del}, a_del);

%calculate the error
if length(frec) > length(f)
    err=norm(frec(1:length(f))-f);
else
    err=norm(frec-f(1:length(frec)));
end
fprintf('Reconstruction error (lin. f-spacing):      %e\n',err);

%compare the frequency spacing of the two filter banks
figure(2)
plot(fc_cq, 'x')
hold on
plot(fc_del, 'o')
xlabel('Number of wavelet filters')
ylabel('Center frequencies [Hz]')
xlim([0 numel(fc_del)])
grid on
legend({'Conventional f-spacing', 'Linear f-spacing'}, 'location', 'northwest')

%% adding 5 compensation filters in the lower frequency range
MC = 5;
fmin = MC/M * fs;
numscales = M-MC+1;

[g_comp,a_comp, fc_comp]=waveletfilters(Ls,'linear',fs,fmin,fmax,numscales,...
    'repeat','delay',delays, 'uniform');

c_comp = filterbank(f, g_comp, a_comp);


lowf_del = numel(fc_del(fc_del<500));
for ii = 1:lowf_del
    c_dlow{ii} = c_del{ii};
end

lowf_comp = numel(fc_comp(fc_comp<500));
for ii = 1:lowf_comp
    c_clow{ii} = c_comp{ii};
end

%investigate the differences in the low frequency range
figure(3)
subplot(2,1,1)
plotfilterbank(c_dlow,a_del(1:lowf_del), fc_del(1:lowf_del))
title('Low-f FB coefficients, single compensation filter')
subplot(2,1,2)
plotfilterbank(c_clow,a_comp(1:lowf_comp), fc_comp(1:lowf_comp))
title('Low-f FB coefficients, 5 compensation filters')

%% now, select the wavelet

wvlt1 = {'cauchy', 257};
[H1, info1]=freqwavelet(wvlt1,1000);
qest1 = info1.fc/info1.bw;

wvlt2 = {'fbsp', 4, 3};
[H2, info2]=freqwavelet(wvlt2,1000);

%estimate their Q-factor, should be roughly the same
qest1 = info1.fc/info1.bw;
qest2 = info2.fc/info2.bw;

figure(4)
plot(H1)
hold on
plot(H2)
xlim([0 100])
ylim([0 1])
grid on
legend1 = sprintf('Cauchy wavelet with Q_est= %0.2f', qest1);
legend2 = sprintf('B-spline wavelet with Q_est= %0.2f', qest2);
legend({legend1, legend2}, 'location', 'northeast')
title('Two wavelets with a similar Q-factor')

%now, increase the Q-factor for wavelet 2
wvlt2 = {'fbsp', 4, 10};
H2=freqwavelet(wvlt2,1000);

%specify a desired target redundancy and delay function
redundancy = 4;
delays = lowdiscrepancy('digital');

%pass the scales directly
%determine the frequency range to be covered
fn_wl = 10; %10 is waveletfilters' internal nyquist frequency
fmax = fn_wl;
freq_step = fn_wl/numscales;
%set fmin depending on the desired start index for the compensation filters
start_index = 3;
fmin = freq_step*start_index; 
scales = 1./linspace(fmin,fmax,numscales);

[g_w1, a_w1,fc_w1,Ls1,info1] = waveletfilters(Ls,scales,...
    wvlt1,'uniform','repeat','energy', 'delay',delays, 'redtar', redundancy);

[g_w2, a_w2,fc_w2,Ls2,info2] = waveletfilters(Ls,scales,...
    wvlt2,'uniform','repeat','energy', 'delay',delays, 'redtar', redundancy);

%compare the frequency response of the two filter banks...
figure(5)
subplot(2,1,1)
filterbankfreqz(g_w1,a_w1,Ls1, 'plot', 'posfreq', 'dynrange', 70);
title('FB frequency response, Cauchy wavelet with small Q-factor')
subplot(2,1,2)
filterbankfreqz(g_w2,a_w2,Ls2, 'plot', 'posfreq', 'dynrange', 70);
title('FB frequency response, B-spline wavelet with large Q-factor')

%...and their coefficients
c_w1 = filterbank(f, g_w1, a_w1);
c_w2 = filterbank(f, g_w2, a_w2);

figure(6)
subplot(2,1,1)
plotfilterbank(c_w1,a_w1, fc_w1)
title('FB coefficients, Cauchy wavelet with small Q-factor')
subplot(2,1,2)
plotfilterbank(c_w2,a_w2, fc_w2)
title('FB coefficients, B-spline wavelet with large Q-factor')
