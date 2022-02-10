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
%   BAUSTELLE---Something about the default wavelet and parameters, also refer to 
%   |freqwavelet| for more information---
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
%   results in most redundant representation which produces nice plots.
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
%   BAUSTELLE: Add some examples, be inspired by audfilters or cqtfilters
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

[varargin,winCell] = arghelper_filterswinparser(definput.flags.wavelettype,varargin);
[flags,kv]=ltfatarghelper({},definput,varargin);
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

% Sort scales for later use
scales_sorted = sort(scales,'descend');

%% Determine total number of filters and natural subsampling factor for lowpass
[aprecise, M, lowpass_number, lowpass_at_zero] = c_det_lowpass(Ls, scales, basea, flags, kv);

%% Compute the downsampling rate
[a, L] = c_comp_downsampling(Ls, M, scales, aprecise, lowpass_at_zero, flags, kv);
% Get an expanded "a"
%afull=comp_filterbank_a(a,M2,struct());

%% Check or compute numeric delay vector
delayvec = c_comp_delay(kv, M, a);

%% Compute the scaling of the filters
[scal, a, delayvec] = c_comp_scaling(a, delayvec, lowpass_at_zero, flags);

if flags.do_complex
    [gout_positive,info_positive] = freqwavelet(winCell,L,scales,'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,'scal',scal(lowpass_number+1:M2),flags.norm);
    [gout_negative,info_negative] = freqwavelet(winCell,L,-flipud(scales),'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,'negative','scal',scal(M2+1:M2+M),flags.norm);
    gout = [gout_positive,gout_negative];
    fields = fieldnames(info_positive);
    info = struct();
    for kk = 1:length(fields)
            info.(fields{kk}) = [info_positive.(fields{kk}),info_negative.(fields{kk})];
    end
elseif flags.do_analytic
    [gout,info] = freqwavelet(winCell,L,scales,'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,'analytic','scal',scal(lowpass_number+1:M2),flags.norm);
else
    [gout,info] = freqwavelet(winCell,L,scales,'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,'scal',scal(lowpass_number+1:M),flags.norm);
end
    
%% Generate lowpass filters if desired
if flags.do_single % Compute single lowpass from frequency response
    lowpass_bandwidth = 0.2/scales_sorted(4); % Twice the center frequency of the fourth scale
    taper_ratio = 1-scales_sorted(4)/scales_sorted(2); % Plateau width is twice the center frequency of the second scale
    [glow,infolow] = wavelet_lowpass(gout,a,L,lowpass_bandwidth,taper_ratio,scal(1),flags);
    gout = [{glow},gout];
    fields = fieldnames(info);
    for kk = 1:length(fields) % Concatenate info and infolow
            info.(fields{kk}) = [infolow.(fields{kk}),info.(fields{kk})];
    end
elseif flags.do_repeat % Lowpass filters are created by repeating smallest scale wavelet with shifted center frequency
    if numel(scales_sorted) < 2 %This function fails for fewer than two elements in scales
        error('%s: Lowpass generation requires at least two scales.',upper(mfilename));
    end
    [glow,infolow] = wavelet_lowpass_repeat(winCell,scales_sorted(1:2),L,lowpass_number,lowpass_at_zero,scal(1:lowpass_number),kv,flags);
    gout = [glow,gout];
    fields = fieldnames(info);
    if flags.do_complex  
        [ghigh,infohigh] = wavelet_lowpass_repeat(winCell,-scales_sorted(1:2),L,lowpass_number,lowpass_at_zero,scal(end:-1:end-lowpass_number+1),kv,flags);
        gout = [gout,ghigh];
        for kk = 1:length(fields) % Concatenate info with infolow and infohigh
                info.(fields{kk}) = [infolow.(fields{kk}),info.(fields{kk}),infohigh.(fields{kk})];
        end
    else
        for kk = 1:length(fields) % Concatenate info and infolow
                info.(fields{kk}) = [infolow.(fields{kk}),info.(fields{kk})];
        end
    end
elseif flags.do_none % No lowpass, do nothing
    %Do nothing
else
    error('%s: This should not happen.',upper(mfilename));
end


% Apply delays (these are now for the lowpasses only)
for kk = 1:numel(gout)
    gout{kk}.delay = delayvec(kk);
end
info.lowpassstart = lowpass_number + 1;%startindex of actual wavelets (tentative)
% Assign fc and adjust for sampling rate 
fc = (kv.fs/2).*info.fc;
end

%% Lowpass filters
function [glow,infolow] = wavelet_lowpass(g,a,L,lowpass_bandwidth,taper_ratio,scal,flags)

Lw = @(L) min(ceil(lowpass_bandwidth*L),L);

% Compute tapered window
P0 = blfilter({'hann','taper',taper_ratio},lowpass_bandwidth,'fs',2,'inf');
% Compute the necessary compensation of the filterbank response
if flags.do_real
    temp_fbresp = @(L) filterbankresponse(g,a(2:end,:),L,'real');
else
    temp_fbresp = @(L) filterbankresponse(g,a(2:end,:),L);
end
Hinv = @(L) sqrt(max(temp_fbresp(L))-temp_fbresp(L));

% Compute the final low-pass filter
glow.H = @(L) fftshift(long2fir(...
filterbankfreqz(P0,1,L).*Hinv(L),Lw(L)))*scal; 
glow.foff = @(L) -floor(Lw(L)/2);
glow.realonly = 0;
glow.delay = 0;

% Initialize and populate infolow
infolow = struct();
infolow.fc = 0;
infolow.foff = glow.foff(L);
infolow.fsupp = Lw(L);
infolow.basefc = 0; % This value has no meaning and is only assigned to prevent errors.
infolow.scale = 0; % This value has no meaning and is only assigned to prevent errors.
infolow.dilation = 0; % This value has no meaning and is only assigned to prevent errors.
infolow.bw = lowpass_bandwidth;
infolow.tfr = 1; % This value has no meaning and is only assigned to prevent errors.
infolow.aprecise = lowpass_bandwidth*L;
infolow.a_natural = [L ceil(infolow.aprecise)];
infolow.a_natural = infolow.a_natural';
infolow.cauchyAlpha = 1; % This value has no meaning and is only assigned to prevent errors.

end

function [glow,infolow] = wavelet_lowpass_repeat(winCell,scales_sorted,L,lowpass_number,lowpass_at_zero,scal,kv,flags)

% Compute frequency range to be covered and step between lowpass filters
LP_range = 0.1/scales_sorted(1);
LP_step = abs(0.1/scales_sorted(2)-LP_range);

% Distinguish between positive and negative scales
if scales_sorted(1)>0
    [glow,infolow] = freqwavelet(winCell,L,repmat(scales_sorted(1),1,lowpass_number),...
        'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,'scal',scal,flags.norm);
    if ~iscell(glow)
        glow = {glow};
    end
    for kk = 1:lowpass_number % Correct foff
        %Flooring sometimes introduces rounding problems when the filter is
        %too slim. 
        %glow{lowpass_number-kk+1}.foff = @(L) glow{lowpass_number-kk+1}.foff(L) - floor(L*kk*LP_step/2);
        glow{lowpass_number-kk+1}.foff = @(L) glow{lowpass_number-kk+1}.foff(L) - round(L*kk*LP_step/2);
        infolow.foff(lowpass_number-kk+1) = glow{lowpass_number-kk+1}.foff(L);
    end
    infolow.fc = LP_range-(lowpass_number:-1:1)*LP_step; % Correct fc
    %if infolow.fc(1) < LP_step/2 %Experimental
    %    glow{1}.H = @(L) glow{1}.H(L)/sqrt(2);
    %end
elseif scales_sorted(1)<0
    if lowpass_at_zero
        lowpass_number = lowpass_number-1;
    end
    [glow,infolow] = freqwavelet(winCell,L,repmat(scales_sorted(1),1,lowpass_number),...
        'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,'scal',scal,'negative',flags.norm);
    if ~iscell(glow)
        glow = {glow};
    end
    for kk = 1:lowpass_number % Correct foff
        glow{kk}.foff = @(L) glow{kk}.foff(L) + floor(L*kk*LP_step/2);
        infolow.foff(kk) = glow{kk}.foff(L);
    end
    infolow.fc = LP_range+(1:lowpass_number)*LP_step; % Correct fc   
else
    error('%s: This should not happen.',upper(mfilename));
end

if numel(infolow.fc) == 0 % This value has no meaning and is only assigned to prevent errors.
    infolow.cauchyAlpha = [];
end

% Correction of TFR
%infolow.tfr = @(L) ones(LP_num,1); %This should be unnecessary, as infolow.tfr is expected to be correct.
  
end

