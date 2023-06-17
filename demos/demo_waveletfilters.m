%DEMO_WAVELETFILTERS  Invertible wavelet filter banks via grid-like sampling
%
%   This demo shows how to generate invertible wavelet filter 
%   banks.
%
%   The following wavelet filter banks are produced:
%
%
%   .. figure::
%
%     The figure shows a non-uniformly sampled invertible filter bank with 
%     conventional, geometric frequency spacing.
%
%
%   .. figure::
%
%     The figure shows a non-uniformly sampled invertible filter bank with 
%     conventional, geometric frequency spacing and compensation 
%     filters that have a constant bandwidth in the low-frequency region.
%
%
%   .. figure::
%
%     The figure shows a uniformly sampled invertible filter bank with linear 
%     frequency spacing where the compensation channel index is passed directly.
%     invertibility is achieved by passing a small delay to the filter
%     generator.
%
%
%   .. figure::
%
%      The figure shows a comparison of the filter center frequency spacing.
%
%
%   .. figure::
%
%     The figure shows a uniformly sampled invertible filter bank with linear 
%     frequency spacing, where the scales are passed directly.
%
%
%   See also: waveletfilters, freqwavelet, lowdiscrepancy, filterbankrealdual

%the input signal
[f,fs]=greasy;
Ls = length(f);
fbounds = 1;

%% filter bank 1: non-uniformly sampled wavelet filterbank, 
%  geometric frequency spacing

%  Retrieve the filters and downsampling factors for a 
%  non-uniformly sampled wavelet filterbank with geometric frequency
%  spacing, using the following parameters:
fmin = 400;
fmax = fs/2;
bins = 12;
redundancy = 128;

% waveletfilters supports various parametrizations of several wavelets,
% with the Cauchy wavelet as the default. Its alpha parameter has a
% correspondence to its Q-factor, allowing for a relatively intuitive
% parametrization. For details on the supported wavelets, check |freqwavelet|
cauchyalpha = 100;

[g_geo, a_geo, fc_geo, L_geo, info_geo] = waveletfilters(Ls,'bins', fs,...
    fmin, fmax, bins, {'cauchy',cauchyalpha}, 'redtar', redundancy);
% Per default, one compensation (lowpass) filter covering the region from 0
% Hz to the wavelet region is added to waveletfilters' output.

% compute redundancy
red_geo = 1/a_geo(1) + 2*sum(1./a_geo(2:end));

% compute coefficients
c_geo=filterbank(f,g_geo,a_geo);

% plot filter bank coefficients, frequency responses and total filter bank
% response
figure(1)
subplot(3,1,1)
plotfilterbank(c_geo, a_geo)
title('Filter bank coefficients')
subplot(3,1,2)
filterbankfreqz(g_geo,a_geo,Ls, 'plot', 'posfreq', 'dynrange', 60);
title('Frequency response single filters')
subplot(3,1,3)
filterbankresponse(g_geo,a_geo,Ls, 'plot', 'real', 'fs', fs);
title('Total filter bank response')

% Compute the frame bounds. The filter bank has nonuniform decimation and the
% wavelet is not bandlimited. Thus, the bounds are estimated by an
% iterative procedure. Iterative computation of the frame bounds is not 
% available directly in the filterbank module, but in the frames module.
% For iteratively reconstructing the time-domain signal from the filter
% bank coefficients, check |ifilterbankiter|.

F = frame('filterbankreal',g_geo, a_geo, numel(g_geo));
[A,B]=framebounds(F,'iter');
FB_ratio_geo = B/A



%% filter bank 2: non-uniformly sampled invertible wavelet filter bank,
%  geometric frequency spacing

%  The flag 'repeat' specifies the usage of more than one compensation 
%  filter in the low frequency range and can improve the frame properties.

fmin = 120;

[g_geored, a_geored, fc_geored, L_geored, info_geored] = waveletfilters(Ls,'bins',...
    fs, fmin, fmax, bins,'repeat');

red_geored = 1/a_geored(1) + 2*sum(1./a_geored(2:end));


fprintf('Redundancy (geom. f-spacing fb 1):                %.2f\n',red_geo);
fprintf('Redundancy (geom. f-spacing fb 2):                %.2f\n',red_geored);


F = frame('filterbankreal',g_geored, a_geored, numel(g_geored));
[A,B]=framebounds(F,'iter');
FB_ratio_geored = B/A

c_geored=filterbank(f,g_geored,a_geored);

figure(2)
subplot(3,1,1)
plotfilterbank(c_geored, a_geored)
title('Filter bank coefficients')
subplot(3,1,2)
filterbankfreqz(g_geored,a_geored,L_geored, 'plot', 'posfreq', 'dynrange', 60);
title('Frequency response single filters')
subplot(3,1,3)
filterbankresponse(g_geored,a_geored,L_geored, 'plot', 'real', 'fs', fs);
title('Total filter bank response')

%% filter bank 3: invertible filter bank with linear frequency spacing

%  The frequency range to be covered is passed as input parameters.
%  Invertibility is achieved via the addition of a small delay to each 
%  wavelet sampling point. The delay is generated
%  from a kronecker low discrepancy sequence and passed to waveletfilters
%  as an anonymous function.
delay = lowdiscrepancy('kronecker');

% specifying a new sampling frequency and redundancy
fs = 2;
redundancy = 4;
% the number of channels overall
M = 1024;
% number of compensation channels
MC = 16;
% number of wavelet channels
wlchannels = M-MC+1;

% the frequency range covered by the wavelets
fmax = fs/2;
fmin = MC/M;

[g_lin,a_lin,fc_lin, L_lin, info_lin]=waveletfilters(Ls,'linear',fs,fmin,fmax,wlchannels,...
    {'cauchy',cauchyalpha},'uniform','redtar', redundancy, 'repeat','delay',delay, 'energy');
% The 'energy' flag is used to scale each filter such that its energy is
% *1*. For further scaling options, see the help of the function |setnorm|.


F = frame('filterbankreal',g_lin, a_lin, numel(g_lin));
[A,B]=framebounds(F,'iter');
FB_ratio_lin = B/A
c_lin=filterbank(f,g_lin,a_lin);

figure(3)
subplot(3,1,1)
plotfilterbank(c_lin, a_lin)
title('Filter bank coefficients')
subplot(3,1,2)
filterbankfreqz(g_lin,a_lin,Ls, 'plot', 'posfreq', 'dynrange', 60);
title('Frequency response single filters')
subplot(3,1,3)
filterbankresponse(g_lin,a_lin,Ls, 'plot', 'real', 'fs', fs);
title('Total filter bank response')


red_lin = 1/a_lin(1) + 2*sum(1./a_lin(2:end));

fprintf('Redundancy (lin. f-spacing fb 3):                %.2f\n',red_lin);


% inversion via the dual system
gd_lin = filterbankrealdual(g_lin, a_lin, L_lin);

% The |ifilterbank| function targets complex filter bank coefficients. For
% energy preservation when inverting real coefficients, taking twice the 
% real part of ifilterbank's result is necessary.
frec_lin = 2*real(ifilterbank(c_lin, gd_lin, a_lin));

% approximating the Q-factor of the wavelet
fprintf('Approximated Q-factor of the Cauchy wavelet with alpha parameter %i: %.2f \n',...
    cauchyalpha, info_lin.fc/info_lin.bw)


% compare the frequency spacing of the filter banks
figure(4)
plot(fc_geo/8000, 'x') %adjust for display purposes with the previous Nyquist f
hold on
plot(fc_lin, 'o')
xlabel('Number of wavelet filters')
ylabel('Center frequencies [Hz]')
xlim([0 numel(fc_geo)])
grid on
legend({'Conventional f-spacing', 'Linear f-spacing'}, 'location', 'northwest')


%% filter bank 4: invertible filter bank with linear frequency spacing

%  Here, the scales are passed directly.
%  waveletfilters' internal sampling frequency:
fs_intern = 20;
fs = 16000;
% The number of desired channels needs to be specified without the desired
% number of compensation channels; here, the same number of compensation
% channels as in filter bank 2 are used
channels = 512-info_geored.waveletstart-1;
fmax = fs/2; %maximum frequency to be covered in Hz
fmax_intern = fmax*fs_intern/fs;
freq_step = fmax_intern/channels;
fmin_intern = freq_step;
fmin = fs/fs_intern * fmin_intern; %minimum frequency to be covered in Hz

% adjust fmin with the start index from filter bank 2
start_index = info_geored.waveletstart;
fmin = fmin*start_index;
fmin_intern = fmin_intern*start_index;
scales = 1./linspace(fmin_intern,fmax_intern,channels);

cauchyalpha = 400;
redundancy = 4;

[g, a,fc,L,info] = waveletfilters(Ls,scales,{'cauchy',cauchyalpha},...
    'uniform','delay',delay, 'repeat', 'redtar', redundancy, 'energy');

% waveletfilters internal fmin = fc(info.waveletstart)*fs/2

fprintf('Approximated Q-factor of the Cauchy wavelet with alpha parameter %i: %.2f \n',...
    cauchyalpha, info.fc/info.bw)

c=filterbank(f,g,a);

figure(5)
subplot(3,1,1)
plotfilterbank(c, a)
title('Filter bank coefficients')
subplot(3,1,2)
filterbankfreqz(g,a,Ls, 'plot', 'posfreq', 'dynrange', 60);
title('Frequency response single filters')
subplot(3,1,3)
filterbankresponse(g,a,Ls, 'plot', 'real', 'fs', fs);
title('Total filter bank response')

F = frame('filterbankreal',g, a, numel(g));
[A,B]=framebounds(F,'iter');
FB_ratio = B/A

if size(a,2) == 2
    a_temp = a(:,1)./a(:,2);
    red = 1/a_temp(1) + 2*sum(1./a_temp(2:end));
else
    red = 1/a(1) + 2*sum(1./a(2:end));
end

% inversion via the dual system
gd = filterbankrealdual(g, a, L);
frec = 2*real(ifilterbank(c, gd, a));

% calculate the reconstruction error
if length(frec) > length(f)
   err=norm(frec(1:length(f))-f);
else
   err=norm(frec-f(1:length(frec)));
end
fprintf('Reconstruction error (lin. f-spacing fb 4):       %e\n',err);