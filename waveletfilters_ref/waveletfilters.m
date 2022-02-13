%function [gout,a,fc,L,info] = waveletfilters(fs, fmin, fmax, M0, Ls, varargin)
function [gout,a,fc,L,info] = waveletfilters(Ls, scales, varargin)
%WAVELETFILTERS Generates wavelet filters
%   Usage: H=freqwavelet(Ls,scales)
%          [H,info]=freqwavelet(...)
%
%   Input parameters:
%         Ls     : System length
%         scales : Vector of wavelet scales
%   Output parameters:
%         gout  : Cell arrary of wavelet filters
%         a     : Downsampling rate for each channel.
%         fc    : Center frequency of each channel.
%         L     : Next admissible length suitable for the generated filters.
%         info  : Struct with additional outputs
%
%   `waveletfilters(Ls,scales)` constructs a system of wavelet filters covering 
%   scales in the range *scales* for system length *Ls*. A scale of 1 corresponds 
%   to a wavelet filter with peak positioned at the frequency 0.1 relative 
%   to the Nyquist rate.
%
%   Wavelet types
%   --------------------
%
%   The following wavelets can be passed as a flag:
%
%   'cauchy'     Cauchy wavelet (default parameters [alpha beta gamma] = [300 0 3])
%
%   'morse'      Generalized morse wavelet (default parameters [alpha beta gamma] = [300 0 3])
%
%   'morlet'     Morlet wavelet (default parameters sigma = [4])
%
%   'fbsp'       Frequency B-spline wavelet (default parameters [order fb] = [4 2])
%
%   'analyticsp' Analytic spline wavelet (default parameters [order fb] = [4 2])
%
%   'cplxsp'     Complex spline wavelet (default parameters [order fb] = [4 2])
%
%   For more details on the construction of the wavelets and the available
%   wavelet types, please see |freqwavelet|. 
%
%   By default, wavelet filters are peak normalized before being adjusted
%   to the proposed downsampling factor. The peak normalization can be 
%   overridden by forwarding any norm flag accepted by |normalize|.
%
%   Downsampling factors
%   --------------------
%
%   The integer downsampling rates of the channels must all divide the
%   signal length, |filterbank| will only work for input signal lengths
%   being multiples of the least common multiple of the downsampling rates.
%   See the help of |filterbanklength|. 
%   The fractional downsampling rates restrict the filterbank to a single
%   length *L=Ls*.
%
%   `[gout,a]=waveletfilters(...,'regsampling')` constructs a non-uniform
%   filterbank with integer subsampling factors. This is the default.
%
%   `[gout,a]=waveletfilters(...,'uniform')` constructs a uniform filterbank
%   where the integer downsampling rate is the same for all the channels. This
%   results in the most redundant representation which produces nice plots.
%
%   `[gout,a]=waveletfilters(...,'fractional')` constructs a filterbank with
%   fractional downsampling rates *a*. 
%   This results in the least redundant system.
%
%   `[gout,a]=waveletfilters(...,'fractionaluniform')` constructs a filterbank 
%   with fractional downsampling rates *a*, which are uniform for all filters
%   except the "filling" low-pass filter which can have different
%   fractional downsampling rates. This is useful when uniform subsampling
%   and low redundancy at the same time are desirable.
%
%   Lowpass filters
%   --------------------
%
%   `[gout,a]=waveletfilters(...,'single')` uses a single lowpass filter
%   for covering the range from zero frequency to the center frequency of the
%   largest scale specified. This is the default.
%
%   `[gout,a]=waveletfilters(...,'repeat')` constructs frequency-shifted
%   copies of the largest scale wavelet to cover the range from zero frequency 
%   to the center frequency of the largest scale specified.
%
%   `[gout,a]=waveletfilters(...,'none')` foregoes the construction of a
%   lowpass filter. This option cannot be expected to yield an invertible 
%   filterbank.
%
%   Additional parameters
%   ---------------------
%
%   `waveletfilter` accepts the following optional parameters:
%
%     'redmul',redmul    Redundancy multiplier. Increasing the value of this
%                        will make the system more redundant by lowering the
%                        channel downsampling rates. It is only used if the
%                        filterbank is a non-uniform filterbank. Default
%                        value is *1*. If the value is less than one, the
%                        system may no longer be painless.
% 
%     'redtar',redtar    Target redundancy. The downsampling factors will be
%                        adjusted to achieve a redundancy as close as possible
%                        to 'redtar'.
%
%     'trunc_at',thr     Applies hard thresholding of the wavelet filters 
%                        at the specified threshold value to reduce their 
%                        support size. 
%                        The default value is *trunc_at=10e-5*. When no 
%                        truncation is desired, *trunc_at=0* should be chosen.
%
%     'delay',delay      A scalar, numeric vector of function handle that 
%                        specifies delays for the wavelet filters. A
%                        numeric vector must have at least as many entries
%                        as there are filters in the filterbank. A function
%                        handle must accept two inputs *(k-1,a(k))*, where 
%                        *k* is the channel index and *a* are the
%                        downsampling rates. If a function handle is given
%                        and 'redtar' is specified, delays are computed
%                        based on the final value of *a*.
%
%     'real'             Allows positive scales with center frequencies up 
%                        to Nyquist. This is the default.
%
%     'complex'          Allows positive scales with center frequencies up 
%                        to Nyquist, which are also mirrored to cover
%                        negative scales.
%
%     'analytic'         Allows positive scales with center frequencies up 
%                        to twice the Nyquist frequency. This setting is
%                        suitable for the analysis of analytic signals.
%
%   Examples:
%   ---------
%
%   In the first example, we analyze a glockenspiel signal with a
%   regularly sampled wavelet filterbank using a frequency B-spline
%   wavelet of order 4 and with parameter fb=3 and visualize the result:::
%
%     [f,fs]=gspi;  % Get the test signal
%     Ls = length(f);
%     scales = linspace(10,0.1,100);
%     [g,a,fc,L]=waveletfilters(Ls,scales, {'fbsp', 4, 3}, 'repeat');
%     c=filterbank(f,g,a);
%     plotfilterbank(c,a,fc,fs,90);
%
%   In the second example, we construct a wavelet filterbank with a
%   lowpass channels based on a Cauchy wavelet and verify it.
%   The plot shows the frequency responses of
%   filters used for analysis (top) and synthesis (bottom). :::
%
%     [f,fs]=greasy;  % Get the test signal
%     Ls = length(f);
%     M0 = 511; %Desired number of channels (without 0 Hz-lowpass channel)
%     max_freqDiv10 = 10;  % 10 corresponds to the nyquist frequency
%     freq_step = max_freqDiv10/M0;
%     rate = 44100;
%     start_index = 1;
%     min_freqHz = rate/10*freq_step
%     min_scale_freq = min_freqHz*start_index
%     min_freqDiv10 = freq_step*start_index; %1/25; % By default, the reference scale for freqwavelet has center frequency 0.1
%     scales = 1./linspace(min_freqDiv10,max_freqDiv10,M0);
%     alpha = 1-2/(1+sqrt(5)); % 1-1/(goldenratio) delay sequence
%     delays = @(n,a) a*(mod(n*alpha+.5,1)-.5);
%     CauchyAlpha = 600;
%     [g, a,fc,L,info] = waveletfilters(Ls,scales,{'cauchy',CauchyAlpha},'uniform','single','energy', 'delay',delays, 'redtar', 8);
%
%     c=filterbank(f,{'realdual',g},a);
%     r=2*real(ifilterbank(c,g,a));
%     if length(r) > length(f)
%         norm(r(1:length(f))-f)
%     else
%         norm(r-f(1:length(r)))
%      end
%     % Plot frequency responses of individual filters
%     gd=filterbankrealdual(g,a,L);
%     figure(1);
%     subplot(2,1,1);
%     filterbankfreqz(gd,a,L,fs,'plot','linabs','posfreq');
%
%     subplot(2,1,2);
%     filterbankfreqz(g,a,L,fs,'plot','linabs','posfreq');
% 
%   See also: freqwavelet, filterbank, normalize

% AUTHORS: Nicki Holighaus, Zdenek Prusa, Guenther Koliander, Clara Hollomey

complainif_notenoughargs(nargin,2,upper(mfilename));
complainif_notposint(Ls,'Ls',upper(mfilename));

definput.import={'normalize'};
definput.importdefaults={'null'};
definput.flags.real = {'real','complex','analytic'};
definput.flags.lowpass  = {'single','repeat','none'};
definput.flags.painless  = {'regular','painless'};
definput.flags.sampling = {'regsampling','uniform',...
                           'fractional','fractionaluniform'};
definput.flags.wavelettype = getfield(arg_freqwavelet(),'flags','wavelettype');
definput.keyvals.redmul=1;
definput.keyvals.redtar=[];
definput.keyvals.delay = 0;
definput.keyvals.trunc_at  = 10^(-5);
definput.keyvals.fs = 2;
definput.keyvals.M0 = 512;
nyquist = 10;%because freqwavelet has internally 0.1 as a reference and fs = 2;
fstep = nyquist/definput.keyvals.M0;
%minf = fs/nyquist*fstep;
%maxf = fs/nyq

%definput.keyvals.scales = 1./linspace(minf,maxf,definput.keyvals.M0);%linearly spaced f?

[varargin,winCell] = arghelper_filterswinparser(definput.flags.wavelettype,varargin);
[flags,kv]=ltfatarghelper({},definput,varargin);
%fs = kv.fs;
%scales = kv.scales;

if isempty(winCell), winCell = {flags.wavelettype}; end

if ~isa(kv.delay,'function_handle') && ~isnumeric(kv.delay)
    error('%s: delay must be a function handle or numeric.',upper(mfilename));
end

if ~isscalar(kv.redmul) || kv.redmul <= 0
    error('%s: redmul must be a positive scalar.',upper(mfilename));
end

if ~isempty(kv.redtar)
    if ~isscalar(kv.redtar) || kv.redtar <= 0
        error('%s: redtar must be a positive scalar.',upper(mfilename));
    end
end

if ~isnumeric(scales) || any(scales <= 0)
    error('%s: scales must be positive and numeric.',upper(mfilename));
end
    
if size(scales,2)>1
    if size(scales,1)==1
        % scales was a row vector.
        scales=scales(:);
    else
        error('%s: scales must be a vector.',upper(mfilename));
    end
end


%% Generate mother wavelet to determine parameters from
[~,info] = freqwavelet(winCell,Ls,1,'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1);
basea = info.aprecise;


%% Determine total number of filters and natural subsampling factor for lowpass
[aprecise, M, lowpass_number, lowpass_at_zero] = c_det_lowpass(Ls, scales, basea, flags, kv);

%% Compute the downsampling rate
[a, L] = c_comp_downsampling(Ls, M, scales, aprecise, lowpass_at_zero, flags, kv);


%% Compute the scaling of the filters and the numeric delay vector

if isa(kv.delay,'function_handle')
    delayvec = zeros(M,1);
    for kk = 1:M
        delayvec(kk) = kv.delay(kk-1,a(kk,1)./a(kk,2));
    end
elseif numel(kv.delay) == 1
    delayvec = repmat(kv.delay,M,1);
elseif ~isempty(kv.delay) && size(kv.delay,2) > 1
    delayvec = kv.delay(:);
else
    error('%s: delay must be scaler or have enough elements to cover all channels.',upper(mfilename));
end
scal=sqrt(a(:,1)./a(:,2));

if flags.do_complex
    
    if lowpass_at_zero
       a=[a;flipud(a(2:end,:))];
       scal=[scal;flipud(scal(2:end))];
       delayvec=[delayvec;flipud(delayvec(2:end))];
    else
        a=[a;flipud(a)];
        scal=[scal;flipud(scal)];
        delayvec=[delayvec;flipud(delayvec)];
    end
    
    [gout_positive,info_positive] = freqwavelet(winCell,L,scales,...
        'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,...
        'scal',scal(lowpass_number+1:M),'delay', delayvec(lowpass_number+1:M),flags.norm);
    [gout_negative,info_negative] = freqwavelet(winCell,L,-flipud(scales),...
        'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,...
        'negative','scal',scal(M+1:M+numel(scales)),'delay', delayvec(M+1:M+numel(scales)), flags.norm);
    gout = [gout_positive,gout_negative];
    fields = fieldnames(info_positive);
    info = struct();
    for kk = 1:length(fields)
            info.(fields{kk}) = [info_positive.(fields{kk}),info_negative.(fields{kk})];
    end
elseif flags.do_analytic
    [gout,info] = freqwavelet(winCell,L,scales,...
        'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,...
        'analytic','scal',scal(lowpass_number+1:M),'delay', delayvec(lowpass_number+1:M),flags.norm);
else
    if lowpass_at_zero
        % Scale the lowpass filters
        scal(1)=scal(1)/sqrt(2);
    end
    
    [gout,info] = freqwavelet(winCell,L,scales,'asfreqfilter','efsuppthr',...
        kv.trunc_at,'basefc',0.1,'scal',scal(lowpass_number+1:M),'delay', delayvec(lowpass_number+1:M),flags.norm);
end
    
%% Generate lowpass filters if desired
[gout, info] = c_make_filters(winCell, gout, a, L, info, scales, scal, lowpass_number, lowpass_at_zero, kv, flags);

% Apply delays (these are now for the lowpasses only)
%for kk = 1:numel(gout)
%    gout{kk}.delay = delayvec(kk);
%end
info.lowpassstart = lowpass_number + 1;%startindex of actual wavelets (tentative)
% Assign fc and adjust for sampling rate 
fc = (kv.fs/2).*info.fc;
end

