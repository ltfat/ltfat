%%test for new user interface
 [f, ~] = gspi;
 Ls = length(f);
 scales = linspace(10,0.1,100);
 fs = 44100;
 fmin = 250;
 fmax = 20000;
 bins = 8;
 M=100;
 
%first call includes the option to set a starting frequency for the wavelet
%frequency range (can be applied to others too)
[g_scales,a_scales,fc_scales,L_scales, info_scales]=waveletfilters(Ls,scales, 'uniform');
[g_bins,a_bins,fc_bins,L_bins,info_bins] = waveletfilters(Ls,'bins', fs,fmin, fmax, bins, 'uniform', 'startfreq', 400);
[g_linear,a_linear,fc_linear,L_linear,info_linear] = waveletfilters(Ls,'linear', fs,fmin, fmax, M, 'uniform');
[g_log,a_log,fc_log,L_log,info_log] = waveletfilters(Ls,'logarithmic', fs,fmin, fmax, M, 'uniform');

Lscales = filterbanklength(L_scales, a_scales);
Lbins = filterbanklength(L_bins, a_bins);
Llinear = filterbanklength(L_linear, a_linear);
Llog = filterbanklength(L_log, a_log);

gd_scales=filterbankrealdual(g_scales,a_scales(:,1),Lscales, 'asfreqfilter');
gd_bins=filterbankrealdual(g_bins,a_bins(:,1),Lbins, 'asfreqfilter');
gd_linear=filterbankrealdual(g_linear,a_linear(:,1),Llinear, 'asfreqfilter');
gd_log=filterbankrealdual(g_log,a_log(:,1),Llog, 'asfreqfilter');
%figure; filterbankfreqz(g_scales,a_scales,Ls,fs,'plot','linabs','posfreq');
%figure; filterbankfreqz(g_bins,a_bins,Ls,fs,'plot','linabs','posfreq');
%figure; filterbankfreqz(g_linear,a_linear,Ls,fs,'plot','linabs','posfreq');
%figure; filterbankfreqz(g_log,a_log,Ls,fs,'plot','linabs','posfreq');


%% basic test waveletfilters
%for ensuring that the external behaviour of the refactored version equals the old one
%addpath('/run/media/clara/3143d7fe-0bef-4c9d-8983-732cfe02d2c9/ltfat/waveletfilters_ref');
pedantic = 1;%switch pedantic on to compare fb coefficients, but this may take a while
 [f, fs] = gspi;
 Ls = length(f);
 scales = linspace(10,0.1,100);

 alpha = 1-2/(1+sqrt(5));
 delays = @(n,a) a*(mod(n*alpha+.5,1)-.5);

[g,a,fc,L, info]=waveletfilters(Ls,scales, {'fbsp', 4, 3}, 'single', 'delay', delays);
[gtemp,atemp,fctemp,Ltemp, infotemp]=waveletfilterstemp(Ls,scales, {'fbsp', 4, 3}, 'single', 'delay', delays);

if pedantic
    %check coefficients
    c = filterbank(f, g, a);
    ctemp = filterbank(f, gtemp, atemp);
    assert(abs(mean(c{2,1}-ctemp{2,1})) == 0)
    assert(sum(c{1,1}-ctemp{1,1}) == 0)
end
%abs(sum(sum(c{1,1}/ctemp{1,1})))/numel(c{1,1})
%check basic params
assert(isequal(a(:,1), atemp(:,1)))
assert(isequal(fc, fctemp))
assert(isequal(L, Ltemp))
%check info struct
assert(isequal(info.fc, infotemp.fc))
assert(isequal(info.basefc, infotemp.basefc))
assert(isequal(info.foff, infotemp.foff))
assert(isequal(info.fsupp, infotemp.fsupp))
assert(isequal(info.scale, infotemp.scale))
assert(isequal(info.dilation, infotemp.dilation))
assert(isequal(info.bw, infotemp.bw))
assert(isequal(info.tfr, infotemp.tfr))
assert(isequal(info.aprecise, infotemp.aprecise))
assert(isequal(info.a_natural, infotemp.a_natural))
assert(isequal(info.cauchyAlpha, infotemp.cauchyAlpha))
assert(isequal(info.lowpassstart, infotemp.startindex))
assert(isequal(g{1,95}.delay, gtemp{1,95}.delay))
assert(isequal(g{1,1}.delay, gtemp{1,1}.delay))

disp('passed all tests')
