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
%      Non-uniformly sampled wavelet filter bank
%
%      The figure shows a non-uniformly sampled invertible filter bank with 
%      conventional, geometric frequency spacing.
%
%
%   .. figure::
%
%      Non-uniformly sampled wavelet filter bank
%
%      The figure shows a non-uniformly sampled invertible filter bank with 
%      geometric frequency spacing, a higher redundancy and using
%      a wavelet with different Q-factor.
%
%
%   .. figure::
%
%      Uniformly sampled wavelet filter bank
%
%      The figure shows a uniformly sampled invertible filter bank with geometric 
%      frequency spacing, where the scales are passed directly as input
%      arguments.
%
%
%   .. figure::
%
%      Uniformly sampled wavelet filter bank
%
%      The figure shows a uniformly sampled invertible filter bank with linear 
%      frequency spacing where the compensation channel index is passed directly.
%      invertibility is achieved by passing a small delay to the filter
%      generator.
%
%
%   .. figure::
%
%      The figure shows a comparison of the filter center frequency spacing.
%
%
%
%
%   See also: waveletfilters, freqwavelet, lowdiscrepancy, filterbankrealdual

% The input signal
[f,fs]=greasy;
Ls = length(f);

%% filter bank 1: Non-uniformly sampled invertible wavelet filterbank, 
%  geometric frequency spacing.

% Retrieve the filters and downsampling factors for a 
% non-uniformly sampled wavelet filterbank with geometric frequency
% spacing, using the following parameters:
fmin = 40;
fmax = fs/2;
bins = 12;

[g_fb1, a_fb1, fc_fb1, L_fb1, info_fb1] = waveletfilters(Ls,'bins', fs,...
    fmin, fmax, bins);
% Per default, one compensation (lowpass) filter covering the region from 0
% Hz to the wavelet region is added to waveletfilters' output.

% Compute the redundancy. Count the DC channel only once, since it produces
% real-valued coefficients
red_fb1 = 1/a_fb1(1) + 2*sum(1./a_fb1(2:end));
fprintf('Redundancy (geom. f-spacing fb 1):                %.2f\n\n',red_fb1);

% Compute the filter bank coefficients, and plot them, along with the
% frequency responses of the single filters and the total filter bank
% response.
c_fb1=filterbank(f,g_fb1,a_fb1);

figure(1), clf
subplot(3,1,1)
plotfilterbank(c_fb1, a_fb1)
title('Filter bank coefficients (fb 1)')
subplot(3,1,2)
filterbankfreqz(g_fb1,a_fb1,L_fb1, 'plot', 'posfreq', 'dynrange', 60);
title('Frequency response single filters')
subplot(3,1,3)
filterbankresponse(g_fb1,a_fb1,L_fb1, 'plot', 'real', 'fs', fs);
title('Total filter bank response')

% Compute the frame bounds. The filter bank has nonuniform decimation and the
% wavelet is not bandlimited. Thus, the bounds are estimated by an
% iterative procedure. Iterative computation of the frame bounds is not 
% available directly in the "filterbank" module, but in the "frames" module.
% For iteratively reconstructing the time-domain signal from the filter
% bank coefficients, check |ifilterbankiter|.

F = frame('filterbankreal',g_fb1, a_fb1, numel(g_fb1));
[A,B]=framebounds(F,L_fb1,'iter','tol',1e-4,'maxit',100,'pcgtol',1e-4,'pcgmaxit',150);
FB_ratio_fb1 = B/A;
fprintf('Frame bound ratio (geom. f-spacing fb 1):                %.2f\n\n',FB_ratio_fb1);


%% filter bank 2: Non-uniformly sampled invertible wavelet filter bank,
%  geometric frequency spacing, with an extended set of input parameters.

% Various parametrizations of several wavelets are possible, using a Cauchy
% wavelet with alpha parameter 300 as the default. The filter banks constructed
% via `waveletfilters` exhibit the constant-Q property, and for Cauchy wavelets,
% the alpha parameter has a correspondence to its Q-factor, allowing for a 
% relatively intuitive parametrization. The properties of the wavelet used
% can be retrieved from the |freqwavelet| function.
 
% Approximate the Q-factor of the wavelet used in filter bank 1.
cauchyalpha_fb1 = 300;
[~,info_wl1] = freqwavelet({'cauchy', cauchyalpha_fb1},L_fb1);
fprintf('Approximated Q-factor of the Cauchy wavelet with alpha parameter %i: %.2f \n\n',...
    cauchyalpha_fb1, info_wl1.fc/info_wl1.bw);
% For filter bank 2, use a wavelet with a higher Q-factor, and approximate
% it for the length of filter bank 1.
cauchyalpha_fb2 = 700;
[~,info_wl2] = freqwavelet({'cauchy', cauchyalpha_fb2},L_fb1);
fprintf('Approximated Q-factor of the Cauchy wavelet with alpha parameter %i: %.2f \n\n',...
    cauchyalpha_fb2, info_wl2.fc/info_wl2.bw);

% Increasing the redundancy via lowering the downsampling rate of the filter 
% bank channels can, for non-uniformly sampled filter banks, improve the frame bounds.
% Here, the downsampling rate is decreased by factor 4. Note that the multiplier 
% is applied before(!) downsampling rates are adapted to the chosen
% decimation scheme. If the decimation scheme is not 'fractional',
% the final redundancy may be higher. 
redundancy_multiplier = 4;

% Construct the set of waveletfilters.
[g_fb2, a_fb2, fc_fb2, L_fb2, info_fb2] = ...
    waveletfilters(Ls,'bins', fs, fmin, fmax, bins, {'cauchy',cauchyalpha_fb2},...
   'redmul', redundancy_multiplier,'fractional');

% Compute the redundancy.
a_long_fb2 = a_fb2(:,1)./a_fb2(:,2);
red_fb2 = 1/a_long_fb2(1) + 2*sum(1./a_long_fb2(2:end));
fprintf('Redundancy (geom. f-spacing fb 2):                %.2f\n\n',red_fb2);

% Compute the filter bank coefficients, and plot them, along with the
% frequency responses of the single filters and the total filter bank
% response.
c_fb2=filterbank(f,g_fb2,a_fb2);

figure(2), clf
subplot(3,1,1)
plotfilterbank(c_fb2, a_fb2)
title('Filter bank coefficients (fb 2)')
subplot(3,1,2)
filterbankfreqz(g_fb2,a_fb2,L_fb2, 'plot', 'posfreq', 'dynrange', 60);
title('Frequency response single filters')
subplot(3,1,3)
filterbankresponse(g_fb2,a_fb2,L_fb2, 'plot', 'real', 'fs', fs);
title('Total filter bank response')

% Compute the frame bounds.
F = frame('filterbankreal',g_fb2, a_fb2, numel(g_fb2));
[A,B]=framebounds(F,L_fb2,'iter','tol',1e-4,'maxit',100,'pcgtol',1e-4,'pcgmaxit',150);
FB_ratio_fb2 = B/A;
fprintf('Frame bound ratio (geom. f-spacing fb 2):                %.2f\n\n',FB_ratio_fb2);

%% filter bank 3: invertible filter bank, specified via wavelet scales

% Besides the specification based on frequency range as above,
% `waveletfilters` allows for passing the scales directly. Scale values 
% greater than 1 make the wavelet wider in time, values lower than 1 make
% the wavelet narrower.
% Since in LTFAT, the normalized sampling frequency is fs = 2 per default, 
% and in |freqwavelet|, wavelets are specified such that their center frequency
% fc = 0.1/scale, the (default) Nyquist frequency of `waveletfilters` 
% is fn = fs_intern/2 with fs_intern = 2 * 1/0.1 = 20 Hz.
fs_intern = 20;

% Similarly, `waveletfilters` does not allow for wavelet center frequencies
% fc that exceed LTFAT's Nyquist frequency fn, hence the smallest possible 
% scale = 0.1/fn = 0.1.
scale_min = log2(0.1);

% The maximum scale can be directly derived from fmin.
fmin_intern = fmin*fs_intern/fs;
scale_max = log2(1/fmin_intern);

% The largest scale accepted by waveletfilters is subject to numerical 
% limitations. The bandwidth of wavelets at low frequencies may become too
% narrow to capture sufficient signal energy. Hence, `waveletfilters`
% per default adds a lowpass channel to the filter bank. See `help
% waveletfilters` for further options.
% To construct a filter bank with the same number of channels as before,
% that additional lowpass filter needs to be taken into account.
channels = numel(fc_fb1) - 1;
scales = 2.^linspace(scale_max,scale_min,channels);

% To ensure invertibility, for uniformly sampled filter banks, a sufficiently 
% high redundancy is advantageous. Note that, similar to the redundancy
% multiplier in fb 2, this is applied before(!) downsampling rates are adapted 
% to the chosen decimation scheme. If the decimation scheme is not 'fractional',
% the final redundancy may differ (usually by being higher). 
redundancy = 64;

% The 'energy' flag is used to scale each filter such that its energy is
% *1*. For further scaling options, see the help of the function |setnorm|.
[g_fb3, a_fb3,fc_fb3,L_fb3,info_fb3] = waveletfilters(Ls,scales,{'cauchy',cauchyalpha_fb1},...
    'redtar', redundancy, 'uniform');

% Compute the redundancy.
red = 1/a_fb3(1) + 2*sum(1./a_fb3(2:end));
fprintf('Redundancy (lin. f-spacing fb 3):                %.2f\n\n',red);

% Compute the filter bank coefficients, and plot them, along with the
% frequency responses of the single filters and the total filter bank
% response.
c_fb3=filterbank(f,g_fb3,a_fb3);

figure(3), clf
subplot(3,1,1)
plotfilterbank(c_fb3, a_fb3)
title('Filter bank coefficients (fb 3)')
subplot(3,1,2)
filterbankfreqz(g_fb3,a_fb3,Ls, 'plot', 'posfreq', 'dynrange', 60);
title('Frequency response single filters')
subplot(3,1,3)
filterbankresponse(g_fb3,a_fb3,Ls, 'plot', 'real', 'fs', fs);
title('Total filter bank response')

% Compute the frame bounds.
[A,B]=filterbankrealbounds(g_fb3,a_fb3,L_fb3);
FB_ratio_fb3 = B/A;
fprintf('Frame bound ratio (lin. f-spacing fb 3):                %.2f\n\n',FB_ratio_fb3);

% For inverting the filter bank, derive the dual system.
L = filterbanklength(Ls, a_fb3);
gd_fb3 = filterbankrealdual(g_fb3, a_fb3, L);
% The |ifilterbank| function targets complex filter bank coefficients. For
% energy preservation when inverting real coefficients, taking twice the 
% real part of ifilterbank's result is necessary.
frec = 2*real(ifilterbank(c_fb3, gd_fb3, a_fb3));

% Calculate the reconstruction error.
if length(frec) > length(f)
   err=norm(frec(1:length(f))-f);
else
   err=norm(frec-f(1:length(frec)));
end
fprintf('Reconstruction error (geom. f-spacing fb 3):       %e\n\n',err);

%% filter bank 4: invertible filter bank with linear frequency spacing

% An invertible wavelet filter bank with linear frequency spacing can be 
% constructed by adding a small delay to each wavelet sampling point. Here, 
% the delay is generated from a Kronecker low discrepancy sequence and passed 
% to waveletfilters as an anonymous function. For details, see
% ltfatnote057.
delay = lowdiscrepancy('kronecker');

% Specifying a new sampling frequency and redundancy
fs = 2;
redundancy = 4;

% Besides the default addition of a lowpass filter at the bottom of the
% filter bank mentioned above (see filter bank 3), it is also possible to
% specify the 'repeat' application of frequency shifted copies of the
% lowest-frequency wavelet. This can further improve the frame bounds at
% the loss of the constant-Q property in the lower frequency range of the
% filter bank. Here, 8 lowpass (compensation) channels are used.
M = 128;
MC = 8;
wlchannels = M-MC+1;
% Calculate the frequency range covered by the wavelets.
fmax = fs/2;
fmin = MC/M;

[g_fb4,a_fb4,fc_fb4, L_fb4, info_fb4]=waveletfilters(Ls,'linear',fs,fmin,fmax,...
    wlchannels,{'cauchy',cauchyalpha_fb1},'uniform','redtar', redundancy,...
    'repeat','delay',delay,'energy');

% Compute the redundancy.
red_fb4 = 1/a_fb4(1) + 2*sum(1./a_fb4(2:end));
fprintf('Redundancy (lin. f-spacing fb 4):                %.2f\n\n',red_fb4);

% Compute the filter bank coefficients, and plot them, along with the
% frequency responses of the single filters and the total filter bank
% response.
c_fb4=filterbank(f,g_fb4,a_fb4);

figure(4), clf
subplot(3,1,1)
plotfilterbank(c_fb4, a_fb4)
title('Filter bank coefficients (fb 4)')
subplot(3,1,2)
filterbankfreqz(g_fb4,a_fb4,L_fb4, 'plot', 'posfreq', 'dynrange', 60);
title('Frequency response single filters')
subplot(3,1,3)
filterbankresponse(g_fb4,a_fb4,L_fb4, 'plot', 'real', 'fs', fs);
title('Total filter bank response')

% Compute the frame bounds.
[A,B]=filterbankrealbounds(g_fb4,a_fb4,L_fb4);
FB_ratio_fb4 = B/A;
fprintf('Frame bound ratio (lin. f-spacing fb 4):                %.2f\n\n',FB_ratio_fb4);


% Inversion via the dual system, as for filter bank 3.
L = filterbanklength(Ls, a_fb4);
gd_fb4 = filterbankrealdual(g_fb4, a_fb4, L);
frec_fb4 = 2*real(ifilterbank(c_fb4, gd_fb4, a_fb4));

% Calculate the reconstruction error.
if length(frec_fb4) > length(f)
   err=norm(frec_fb4(1:length(f))-f);
else
   err=norm(frec_fb4-f(1:length(frec_fb4)));
end
fprintf('Reconstruction error (lin. f-spacing fb 4):       %e\n\n',err);


% Compare the frequency spacing of the filter banks.
figure(5), clf
plot(fc_fb3, 'x') %adjust for display purposes with the previous Nyquist f
hold on
plot(fc_fb4, 'o')
xlabel('Number of wavelet filters')
ylabel('Center frequencies [Hz]')
xlim([0 numel(fc_fb4)])
grid on
legend({'Conventional f-spacing', 'Linear f-spacing'}, 'location', 'northwest')