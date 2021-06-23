function [gout, info] = waveletfilters(L, scales, varargin)
%WAVELETFILTERS Generates wavelet filters
%   Usage: H=freqwavelet(L,scales)
%          [H,info]=freqwavelet(...)
%
%   Input parameters:
%         L     : System length
%         scale : Vector of wavelet scales
%   Output parameters:
%         gout  : Cell arrary of wavelet filters
%         info  : Struct with additional outputs
%
%   `waveletfilters(L,scales)` is a wrapper for creating a system of wavelet
%   filters covering scales in the range *scales* for system length *L*. 
%   The output *gout* is ready to be used within the LTFAT filter bank 
%   framework. By default, a scale of 1 corresponds to a wavelet filter with
%   peak positioned at the frequency 0.1 relative to the Nyquist rate (fs/2).
%
%   See also: freqwavelet, filterbank

% AUTHORS: Zdenek Prusa, Nicki Holighaus, Guenther Koliander, Clara Hollomey


complainif_notenoughargs(nargin,2,upper(mfilename));
complainif_notposint(L,'L',upper(mfilename));

definput.flags.real        = {'real','complex'};
definput.keyvals.trunc_at  = 10^(-5);
definput.flags.wavelettype = getfield(arg_freqwavelet(),'flags','wavelettype');

[varargin,winCell] = arghelper_filterswinparser(definput.flags.wavelettype,varargin);
[flags,kv]=ltfatarghelper({},definput,varargin);
if isempty(winCell), winCell = {flags.wavelettype}; end

a = 1;

[gout,info] = freqwavelet(winCell,L,scales,'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1);
info.tfr = @(Ltrue) [1,info.tfr./L*Ltrue];
info.fc = [0,info.fc];

fsupp_LP = 0.2/info.scale(4); % 2 times center frequency of 8th scale
ratio = 1-info.scale(4)/info.scale(2); % Transition begins at 4th scale center frequency
Lw = @(L) min(ceil(fsupp_LP*L),L);

P0 = blfilter({'hann','taper',ratio},fsupp_LP*2,'fs',2,'inf');
temp_fbresp = @(L) filterbankresponse(gout,a,L,'real');
Hinv = @(L) sqrt(max(temp_fbresp(L))-temp_fbresp(L));

% Compute the final low-pass filter
glow.H = @(L) fftshift(long2fir(...
filterbankfreqz(P0,a,L).*Hinv(L),Lw(L)))*sqrt(a)/sqrt(2); % /sqrt(2) because real
glow.foff = @(L) -floor(Lw(L)/2);
glow.realonly = 0;
glow.delay = 0;

gout = [{glow},gout];


