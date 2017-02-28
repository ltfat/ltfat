function tfr = getgausstfr_warp(freqtoscale,scaletofreq,fc,fs,L,varargin)

definput.keyvals.bwmul=1;
[flags,kv]=ltfatarghelper({},definput,varargin);

%Determine Gauss ERB width
% probelen = 10000;
% probebw = 0.01; 
% 
% % Determine where to truncate the window
% H = freqwin('gauss',probelen,probebw);
% winbw = norm(H).^2/(probebw*probelen/2);

winbw = 0.754;

basebw = 1.875657;

bw = (scaletofreq(freqtoscale(fc)+.5*kv.bwmul) ...
                               -scaletofreq(freqtoscale(fc)-.5*kv.bwmul))/winbw;
bw=2.*bw./fs;

tfr = (bw./basebw).^2.*L;
tfr = 1./tfr;
tfr(tfr == Inf) = 0;