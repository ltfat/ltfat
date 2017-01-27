function tfr = getgausstfr(fc,fs,L,varargin)

definput.flags.audscale={'erb','erb83','bark','mel','mel1000'};
definput.keyvals.bwmul=[];
[flags,kv]=ltfatarghelper({},definput,varargin);
switch flags.audscale
    case {'mel','mel1000'} % The mel scales are very fine, therefore default spacing is adjusted
        definput.keyvals.bwmul=100;
     otherwise
        definput.keyvals.bwmul=1;
end
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

bw=audfiltbw(fc,flags.audscale)/winbw*kv.bwmul;
bw=2.*bw./fs;

tfr = (bw./basebw).^2.*L;
tfr = 1./tfr;