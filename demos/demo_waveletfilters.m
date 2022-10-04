%DEMO_WAVELETFILTERS  Introduction to grid-like wavelet sampling
%
%   This demo shows how to generate invertible wavelet filter 
%   banks.
%
%   The following wavelet filter banks are produced:
%
%     * a non-uniformly sampled invertible filter bank with conventional, 
%     geometric frequency spacing.
%
%     * a non-uniformly sampled invertible filter bank with conventional, 
%     geometric frequency spacing and adjusted redundancy and additional 
%     compensation filters for simpler invertibility.
%
%     * a uniformly sampled invertible filter bank with linear frequency
%     spacing. 
%
%     * a uniformly sampled invertible filter bank with linear frequency
%     spacing, where the scales are passed directly.
%
%   See also: waveletfilters, freqwavelet, lowdiscrepancy, filterbankrealdual

%the input signal
[f,fs]=greasy;
Ls = length(f);

%% filter bank 1: retrieve the filters and downsampling factors for a 
%  wavelet filterbank with the following parameters:
fmin = 40;
fmax = fs/2;
bins = 16;
[g_geo, a_geo, fc_geo, L_geo, info_geo] = waveletfilters(Ls,'bins', fs,...
    fmin, fmax, bins);

% compute redundancy
if size(a_geo,2) == 2
    a_temp = a_geo(:,1)./a_geo(:,2);
    red_geo = 1/a_temp(1) + 2*sum(1./a_temp(2:end));
else
    red_geo = 1/a_geo(1) + 2*sum(1./a_geo(2:end));
end

% compute coefficients
c_geo=filterbank(f,g_geo,a_geo);

% plot filter bank coefficients, frequency responses and total filter bank
% response
figure(1)
subplot(3,1,1)
plotfilterbank(c_geo, a_geo)
subplot(3,1,2)
filterbankfreqz(g_geo,a_geo,Ls, 'plot', 'posfreq', 'dynrange', 60);
subplot(3,1,3)
filterbankresponse(g_geo,a_geo,Ls, 'plot', 'real', 'fs', fs);

% compute frame bounds, the filter bank has nonuniform decimation and the
% wavelet is not bandlimited, such that the bounds are estimated by an
% iterative procedure. Iterative computation is not available directly in 
% the filterbank module, but in the frames module.

F = frame('filterbankreal',g_geo, a_geo, numel(g_geo));
[A,B]=framebounds(F,'iter');
FB_ratio_geo = B/A

%% filter bank 2: waveletfilters supports the direct setting of an 
%  (approximate) redundancy for the filter bank. the flag 'repeat' 
%  specifies the usage of more than one compensation filter in the 
%  low frequency range
redundancy = 32;
bins = 8;

[g_geored, a_geored, fc_geored, L_geored, info_geored] = waveletfilters(Ls,'bins',...
    fs, fmin, fmax, bins, 'redtar', redundancy, 'repeat');


if size(a_geo,2) == 2
    a_temp = a_geored(:,1)./a_geored(:,2);
    red_geored = 1/a_temp(1) + 2*sum(1./a_temp(2:end));
else
    red_geored = 1/a_geored(1) + 2*sum(1./a_geored(2:end));
end

fprintf('Default redundancy (geom. f-spacing fb 1):                %.2f\n',red_geo);
fprintf('Adjusted redundancy (geom. f-spacing fb 2):                %.2f\n',red_geored);


F = frame('filterbankreal',g_geored, a_geored, numel(g_geored));
[A,B]=framebounds(F,'iter');
FB_ratio_geored = B/A

c_geored=filterbank(f,g_geored,a_geored);
figure(2)
subplot(3,1,1)
plotfilterbank(c_geored, a_geored)
subplot(3,1,2)
filterbankfreqz(g_geored,a_geored,L_geored, 'plot', 'posfreq', 'dynrange', 60);
subplot(3,1,3)
filterbankresponse(g_geored,a_geored,L_geored, 'plot', 'real', 'fs', fs);

%% filter bank 3: construct a linearly spaced filterbank using twice the 
%  number of channels and a lower redundancy than before. its invertibility 
%  is achieved by the addition of a small delay to each wavelet sampling point. 
channels = numel(fc_geo)*2;
delay = lowdiscrepancy('kronecker');
fmin = 250;%the range below fmin is covered by compensation filters
cauchyalpha = 100;

[g_lin, a_lin, fc_lin, L_lin, info_lin] = waveletfilters(Ls,'linear',...
    fs, fmin, fmax, channels, {'cauchy', cauchyalpha}, 'uniform', 'delay',...
    delay, 'redtar', redundancy/6,  'repeat', 'energy');
% the 'energy' flag is used to scale each filter such that its energy is
% *1*. for further scaling options, see the help of the function |setnorm|.

F = frame('filterbankreal',g_lin, a_lin, numel(g_lin));
[A,B]=framebounds(F,'iter');
FB_ratio_lin = B/A

if size(a_lin,2) == 2
    a_temp = a_lin(:,1)./a_lin(:,2);
    red_lin = 1/a_temp(1) + 2*sum(1./a_temp(2:end));
else
    red_lin = 1/a_lin(1) + 2*sum(1./a_lin(2:end));
end
fprintf('Redundancy (lin. f-spacing fb 3):                %.2f\n',red_lin);

c_lin=filterbank(f,g_lin,a_lin);

% inversion via the dual system is possible
gd_lin = filterbankrealdual(g_lin, a_lin, L_lin);
% the |ifilterbank| function targets complex filter bank coefficients. for
% energy preservation when inverting real coefficients, taking 2 times the 
% real part of ifilterbank's result is necessary.
frec_lin = 2*real(ifilterbank(c_lin, gd_lin, a_lin));

% for the Cauchy wavelet, its first hyperparameter is proportional to its
% Q-factor
fprintf('Approximated Q-factor of the Cauchy wavelet with hyperparameter %i: %.2f \n',...
    cauchyalpha, info_lin.fc/info_lin.bw)

figure(3)
subplot(3,1,1)
plotfilterbank(c_lin, a_lin)
subplot(3,1,2)
filterbankfreqz(g_lin,a_lin,Ls, 'plot', 'posfreq', 'dynrange', 60);
subplot(3,1,3)
filterbankresponse(g_lin,a_lin,Ls, 'plot', 'real', 'fs', fs);

% compare the frequency spacing of the filter banks
figure(4)
plot(fc_geo, 'x')
hold on
plot(fc_lin, 'o')
xlabel('Number of wavelet filters')
ylabel('Center frequencies [Hz]')
xlim([0 numel(fc_lin)])
grid on
legend({'Conventional f-spacing', 'Linear f-spacing'}, 'location', 'northwest')

%% filter bank 4: to specify the wavelet scales directly, it is necessary 
%  to know waveletfilters' internal sampling frequency:
fs_intern = 20;

% the number of desired channels needs to be specified without the desired
% number of compensation channels; here, the same number of compensation
% channels as in filter bank 2 should be used
channels = 256-info_geored.lowpassstart;
fmax = fs/2; %maximum frequency to be covered in Hz
fmax_intern = fmax*fs_intern/fs;
freq_step = fmax_intern/channels;
fmin_intern = freq_step;
fmin = fs/fs_intern * fmin_intern; %minimum frequency to be covered in Hz

% adjust fmin with the start index from filter bank 2
start_index = info_geored.lowpassstart;
fmin = fmin*start_index;
fmin_intern = fmin_intern*start_index;
scales = 1./linspace(fmin_intern,fmax_intern,channels);

cauchyalpha = 200;

[g, a,fc,L,info] = waveletfilters(Ls,scales,{'cauchy',cauchyalpha},...
    'uniform','delay',delay,'energy', 'repeat', 'redtar', 4);

fprintf('Approximated Q-factor of the Cauchy wavelet with hyperparameter %i: %.2f \n',...
    cauchyalpha, info.fc/info.bw)
fprintf('Approximated Q-factor of the Cauchy wavelet with hyperparameter 300: %.2f \n',...
    info_geo.fc/info_geo.bw)

c=filterbank(f,g,a);

figure(5)
subplot(3,1,1)
plotfilterbank(c, a)
subplot(3,1,2)
filterbankfreqz(g,a,Ls, 'plot', 'posfreq', 'dynrange', 60);
subplot(3,1,3)
filterbankresponse(g,a,Ls, 'plot', 'real', 'fs', fs);

F = frame('filterbankreal',g, a, numel(g));
[A,B]=framebounds(F,'iter');
FB_ratio = B/A

if size(a,2) == 2
    a_temp = a(:,1)./a(:,2);
    red = 1/a_temp(1) + 2*sum(1./a_temp(2:end));
else
    red = 1/a(1) + 2*sum(1./a(2:end));
end

% inversion via the dual system is possible
gd = filterbankrealdual(g, a, L);
frec = 2*real(ifilterbank(c, gd, a));

% calculate the reconstruction error
if length(frec) > length(f)
   err=norm(frec(1:length(f))-f);
else
   err=norm(frec-f(1:length(frec)));
end
fprintf('Reconstruction error (lin. f-spacing fb 4):       %e\n',err);