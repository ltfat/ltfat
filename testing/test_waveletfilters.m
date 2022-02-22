function test_failed=test_waveletfilters
%TEST_WAVELETFILTERS  Test the erbfilters filter generator

%%test for new user interface
[f, fs] = gspi;
f=f(1:1000);
Ls = length(f);
scales = linspace(10,0.1,100);
fmin = 250;
fmax = 20000;
bins = 8;
M=100;
alpha = 1-2/(1+sqrt(5)); % 1-1/(goldenratio) delay sequence
delays = @(n,a) a*(mod(n*alpha+.5,1)-.5);
 
%first call includes the option to set a starting frequency for the wavelet
%frequency range (can be applied to others too)
[g_scales,a_scales,fc_scales,L_scales, info_scales]=waveletfilters(Ls,scales, 'delay', delays, 'uniform', 'redtar', 4, 'complex');
[g_bins,a_bins,fc_bins,L_bins,info_bins] = waveletfilters(Ls,'bins', fs,fmin, fmax, bins, 'uniform', 'startfreq', 800);
[g_linear,a_linear,fc_linear,L_linear,info_linear] = waveletfilters(Ls,'linear', fs,fmin, fmax, M, 'uniform');


Lscales = filterbanklength(L_scales, a_scales);
Lbins = filterbanklength(L_bins, a_bins);
Llinear = filterbanklength(L_linear, a_linear);


gd_scales=filterbankdual(g_scales,a_scales,Lscales, 'asfreqfilter');
gd_bins_asf=filterbankrealdual(g_bins,a_bins,Lbins, 'asfreqfilter');
%gd_bins_e=filterbankrealdual(g_bins,a_bins,Lbins, 'econ');
gd_linear=filterbankrealdual(g_linear,a_linear,Llinear, 'asfreqfilter');

if 0
    % Inspect it: Dual windows, frame bounds and the response
    disp('Frame bounds scales:')
    [A,B]=filterbankbounds(g_scales,a_scales,Lscales);
    A
    B
    B/A
    filterbankresponse(g_scales,a_scales,Lscales,'real','plot');
    disp('Frame bounds bins:')
    [A,B]=filterbankrealbounds(g_bins,a_bins,Lbins);
    A
    B
    B/A
    filterbankresponse(g_bins,a_bins,Lbins,'real','plot');
    disp('Frame bounds linear:')
    [A,B]=filterbankrealbounds(g_linear,a_linear,Llinear);
    A
    B
    B/A
    filterbankresponse(g_linear,a_linear,Llinear,'real','plot');
    figure; filterbankfreqz(g_scales,a_scales,Ls,fs,'plot','linabs','posfreq');
    figure; filterbankfreqz(g_bins,a_bins,Ls,fs,'plot','linabs','posfreq');
    figure; filterbankfreqz(g_linear,a_linear,Ls,fs,'plot','linabs','posfreq');
end

c = filterbank(f, g_bins, a_bins);
fhat = 2*real(ifilterbank(c, gd_bins_asf, a_bins));

if length(fhat) > length(f)
    res=norm(fhat(1:length(f)) - f);
else
    res=norm(f(1:length(fhat)) - fhat);
end

[test_failed,fail]=ltfatdiditfail(res,0, 0.0001);
%s=sprintf(['WAVELETFILTER DUAL  %s %s %s %s L:%3i %0.5g %s'],realcomplex,warping,fracname,uniform,L,res,fail);    
%disp(s);