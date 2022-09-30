%DEMO_WAVELETFILTERS  Introduction to grid-like wavelet sampling
%
%   This demo shows how to generate uniformly sampled invertible wavelet filter 
%   banks.
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
%   * the example from the paper "Grid-based Decimation for Invertible Wavelet
%   Transforms in Audio Processing", Section 4-C.
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
fmin = 40;
fmax = fs/2;
bins = 16;
[g_geo, a_geo, fc_geo, L_geo, info_geo] = waveletfilters(Ls,'bins', fs,...
    fmin, fmax, bins);

% Compute redundancy
if size(a_geo,2) == 2
    a_temp = a_geo(:,1)./a_geo(:,2);
    red_geo = 1/a_temp(1) + 2*sum(1./a_temp(2:end));
else
    red_geo = 1/a_geo(1) + 2*sum(1./a_geo(2:end));
end

% Compute coefficients
c_geo=filterbank(f,g_geo,a_geo);

% Plot filter bank coefficients, frequency responses and total filter bank
% response
figure(1)
subplot(3,1,1)
plotfilterbank(c_geo, a_geo)
subplot(3,1,2)
filterbankfreqz(g_geo,a_geo,Ls, 'plot', 'posfreq', 'dynrange', 60);
subplot(3,1,3)
filterbankresponse(g_geo,a_geo,Ls, 'plot', 'real', 'fs', fs);

% Compute frame bounds, the filter bank has nonuniform decimation and the
% wavelet is not bandlimited, such that the bounds are estimated by an
% iteration. Iterative computation is not available directly in the filterbank 
% module, but in the frames module

F = frame('filterbankreal',g_geo, a_geo, numel(g_geo));
[A,B]=framebounds(F,'iter');
FB_ratio_geo = B/A

%% now, construct a more complex filter bank...
fmax = fs/2;
fmin = 625;
bins = 8;
redundancy = 5;

[g_geored, a_geored, fc_geored, L_geored, info_geored] = waveletfilters(Ls,'bins',...
    fs, fmin, fmax, bins, 'redtar', redundancy, 'energy');

c_geored=filterbank(f,g_geored,a_geored);

if size(a_geo,2) == 2
    a_temp = a_geored(:,1)./a_geored(:,2);
    red_geored = 1/a_temp(1) + 2*sum(1./a_temp(2:end));
else
    red_geored = 1/a_geored(1) + 2*sum(1./a_geored(2:end));
end

figure(2)
subplot(3,1,1)
plotfilterbank(c_geored, a_geored)
subplot(3,1,2)
filterbankfreqz(g_geored,a_geored,Ls, 'plot', 'posfreq', 'dynrange', 60);
subplot(3,1,3)
filterbankresponse(g_geored,a_geored,Ls, 'plot', 'real', 'fs', fs);

F = frame('filterbankreal',g_geored, a_geored, numel(g_geored));
[A,B]=framebounds(F,'iter');
FB_ratio_geored = B/A

% the alternative scaling of the filters may yield a favourable filter bank
% response for some applications. also, the redundancy can be adjusted, but
% this may adversely affect the framebounds, and thus, invertibility.

%% for yet more freedom in the filter bank design,...
channels = numel(fc_geo);
delay = lowdiscrepancy('digital');

[g_lin, a_lin, fc_lin, L_lin, info_lin] = waveletfilters(Ls,'linear',...
    fs, fmin, fmax, channels, {'cauchy', 100}, 'uniform', 'delay', delay,...
    'redtar', 6, 'energy');

F = frame('filterbankreal',g_lin, a_lin, numel(g_lin));
[A,B]=framebounds(F,'iter');
FB_ratio_lin = B/A;

if size(a_lin,2) == 2
    a_temp = a_lin(:,1)./a_lin(:,2);
    red_lin = 1/a_temp(1) + 2*sum(1./a_temp(2:end));
else
    red_lin = 1/a_lin(1) + 2*sum(1./a_lin(2:end));
end

c_lin=filterbank(f,g_lin,a_lin);
figure(4)
subplot(3,1,1)
plotfilterbank(c_lin, a_lin)
subplot(3,1,2)
filterbankfreqz(g_lin,a_lin,Ls, 'plot', 'posfreq', 'dynrange', 60);
subplot(3,1,3)
filterbankresponse(g_lin,a_lin,Ls, 'plot', 'real', 'fs', fs);
% compare the frequency spacing of the two filter banks
figure(2)
plot(fc_geo, 'x')
hold on
plot(fc_lin, 'o')
xlabel('Number of wavelet filters')
ylabel('Center frequencies [Hz]')
xlim([0 numel(fc_lin)])
grid on
legend({'Conventional f-spacing', 'Linear f-spacing'}, 'location', 'northwest')

%% 
fs_intern = 20;

channels = 256;
fmax = fs/2;
fmax_intern = fmax*fs_intern/fs;
freq_step = fmax_intern/channels;
fmin_intern = freq_step;
fmin = fs/fs_intern * fmin_intern;
%adjust fmin with the start index
start_index = 3;
fmin = fmin*start_index;
fmin_intern = fmin_intern*start_index;
scales = 1./linspace(fmin_intern,fmax_intern,channels);
CauchyAlpha = 700;


[g, a,fc,L,info] = waveletfilters(Ls,scales,{'cauchy',CauchyAlpha},...
    'uniform','energy', 'delay',delay, 'redtar', 8);

c=filterbank(f,g,a);
figure(4)
subplot(3,1,1)
plotfilterbank(c_lin, a_lin)
subplot(3,1,2)
filterbankfreqz(g_lin,a_lin,Ls, 'plot', 'posfreq', 'dynrange', 60);
subplot(3,1,3)
filterbankresponse(g_lin,a_lin,Ls, 'plot', 'real', 'fs', fs);

%F = frame('filterbankreal',g, a, numel(g));
%[A,B]=framebounds(F,'iter');
%FB_ratio = B/A

if size(a_lin,2) == 2
    a_temp = a(:,1)./a(:,2);
    red = 1/a_temp(1) + 2*sum(1./a_temp(2:end));
else
    red = 1/a(1) + 2*sum(1./a(2:end));
end