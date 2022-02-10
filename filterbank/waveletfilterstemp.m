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
%     [g,a,fc,L]=waveletfilters(Ls,scales, {'fbsp', 4, 3});
%     c=filterbank(f,g,a);
%     plotfilterbank(c,a,fc,fs,90);
%
%   In the second example, we construct a wavelet filterbank with several
%   lowpass channels based on a Cauchy wavelet and verify it.
%   The plot shows frequency responses of
%   filters used for analysis (top) and synthesis (bottom). :::
%
%     [f,fs]=greasy;  % Get the test signal
%     Ls = length(f);
%     M0 = 511; %Desired number of channels (without 0 Hz-lowpass channel)
%     max_freqDiv10 = 10;  % 10 corresponds to the nyquist frequency
%     freq_step = max_freqDiv10/M0;
%     rate = 44100;
%     min_freqHz = rate/20*freq_step
%     start_index = 10;
%     min_scale_freq = min_freqHz*start_index
%     min_freqDiv10 = freq_step*start_index; %1/25; % By default, the reference scale for freqwavelet has center frequency 0.1
%     scales = 1./linspace(min_freqDiv10,max_freqDiv10,M0-start_index+1);
%     alpha = 1-2/(1+sqrt(5)); % 1-1/(goldenratio) delay sequence
%     delays = @(n,a) a*(mod(n*alpha+.5,1)-.5);
%     CauchyAlpha = 600;
%     [g, a,fc,L,info] = waveletfilters(Ls,scales,{'cauchy',CauchyAlpha},'uniform','repeat','energy', 'delay',delays, 'redtar', 8);
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

% We probably do not care about this. 
% if kv.redtar <= 1
%     warning('%s: redtar is very low; the resulting system might be unstable.',upper(mfilename));
% end

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
% Sort scales for later use
scales_sorted = sort(scales,'descend');

% Get number of scales
M = numel(scales); 

%% Generate mother wavelet to determine parameters from
[~,info] = freqwavelet(winCell,Ls,1,'asfreqfilter', 'efsuppthr',kv.trunc_at,'basefc',0.1);
basea = info.aprecise;

lowpass_at_zero = 0; %Set to 1 if there is a lowpass filter at zero frequency
%% Determine total number of filters and natural subsampling factor for lowpass
if flags.do_repeat 
    lowpass_number = scales_sorted(2)/(scales_sorted(1)-scales_sorted(2)); % Maybe adjust this to not guarantee some distance between first filter and zero frequency.
    if abs(lowpass_number - round(lowpass_number)) < eps*10^3
        lowpass_number = round(lowpass_number);
        lowpass_at_zero = 1;
    end
    lowpass_number = floor(lowpass_number);
    if lowpass_number == 0
        lowpass_number = 1;
    end
    M2 = M + lowpass_number;
    aprecise = (basea.*scales_sorted(1))*ones(lowpass_number,1);
elseif flags.do_single
    lowpass_number = 1;
    lowpass_at_zero = 1;
    M2 = M+1;
    aprecise = (0.2./scales_sorted(4))*Ls; % This depends on how single lowpass is called (l.195ff). Maybe automate. Adapt if necessary!!!
else
    lowpass_number = 0;%this should never happen
    M2 = M;
    aprecise = [];
end

%% Get subsampling factors
aprecise = [aprecise;basea.*scales];

if any(aprecise<1)
    error(['%s: Bandwidth of one of the filters is bigger than fs. ',...
           'This should not happen'],upper(mfilename));
end

aprecise=aprecise/kv.redmul;
if any(aprecise<1)
    error('%s: The maximum redundancy mult. for this setting is %5.2f',...
         upper(mfilename), min(basea./scales));
end

%% Compute the downsampling rate
if flags.do_regsampling % This should only be used for lowpass = single!!!
    a = ones(M2,1);
    
    [lower_scale,~] = max(scales);
    [upper_scale,~] = min(scales);
    lower_scale = floor(log2(1/lower_scale));
    upper_scale = floor(log2(1/upper_scale));
    
    % Find minimum a in each octave and floor23 it.
    for kk = lower_scale:upper_scale
        tempidx = find( floor(log2(1./scales)) == kk );
        [~,tempminidx] = min(1/scales(tempidx));
        idx = tempidx(tempminidx);
        
        % Deal the integer subsampling factors
        a(tempidx) = floor23(aprecise(idx));
    end   
    
    % Determine the minimal transform length lcm(a)
    L = filterbanklength(Ls,a);
    
    % Heuristic trying to reduce lcm(a)
    while L>2*Ls && ~(all(a==a(1)))
        maxa = max(a);
        a(a==maxa) = 0;
        a(a==0) = max(a);
        L = filterbanklength(Ls,a);
    end
    
elseif flags.do_fractional
    L = Ls;
    N=ceil(Ls./aprecise);
    a=[repmat(Ls,M2,1),N];
elseif flags.do_fractionaluniform
    L = Ls;
    if lowpass_at_zero
        aprecise(2:end)= min(aprecise(2:end));
    else 
        aprecise= repmat(min(aprecise),numel(aprecise),1);
    end
    N=ceil(Ls./aprecise);
    a=[repmat(Ls,M2,1),N];
elseif flags.do_uniform
    a=floor(min(aprecise));
    L=filterbanklength(Ls,a);
    a = repmat(a,M2,1);
end

% Get an expanded "a"
afull=comp_filterbank_a(a,M2,struct());

%% Check or compute numeric delay vector
if isa(kv.delay,'function_handle')
    delayvec = zeros(M2,1);
    for kk = 1:M2
        delayvec(kk) = kv.delay(kk-1,afull(kk,1)./afull(kk,2));
    end
elseif numel(kv.delay) == 1
    delayvec = repmat(kv.delay,M2,1);
elseif any(size(kv.delay,2)) == 1 && numel(kv.delay) >= M2
    delayvec = kv.delay(:);
else
    error('%s: delay must either be a scalar or have enough elements to cover all channels.',upper(mfilename));
end

%% Compute the scaling of the filters
% Individual filter peaks are made square root of the subsampling factor
scal=sqrt(afull(:,1)./afull(:,2));

if flags.do_real && lowpass_at_zero
    % Scale the lowpass filters
    scal(1)=scal(1)/sqrt(2);
elseif flags.do_complex
    % Replicate the scaling, sampling rates and delays, except at zero
    % frequency
    if lowpass_at_zero
       a=[a;flipud(a(2:end,:))];
       scal=[scal;flipud(scal(2:end))];
       delayvec=[delayvec;flipud(delayvec(2:end))];
    else
        a=[a;flipud(a)];
        scal=[scal;flipud(scal)];
        delayvec=[delayvec;flipud(delayvec)];
    end
end

%% Adjust the downsampling rates in order to achieve 'redtar' [It may be possible to do this up front and not run filterbankscale]
if ~isempty(kv.redtar)
    
    if size(a,2) == 2
        a_old = a(:,1)./a(:,2);
    else
        a_old = a;
    end
    
    if ~flags.do_real
        org_red = sum(1./a_old);
    elseif lowpass_at_zero
        org_red = 1./a_old(1) + sum(2./a_old(2:end));
    else
        org_red = sum(2./a_old);
    end
    
    a_new = floor(a*org_red/kv.redtar);
    scal_new = org_red/kv.redtar*ones(1,M2);

    % Adjust function handle generated delays
    if isa(kv.delay,'function_handle')
        delayvec = zeros(M2,1);
        for kk = 1:M2
            delayvec(kk) = kv.delay(kk-1,a_new(kk));
        end
        if flags.do_complex
            % Replicate the delays, except at zero
            % frequency
            if lowpass_at_zero
                delayvec=[delayvec;flipud(delayvec(2:end))];
            else
                delayvec=[delayvec;flipud(delayvec)];
            end
        end
    end
    
    if ~flags.do_uniform
        N_old = ceil(L./a_old);
        N_new=ceil(L./a_new);
        a_new=[repmat(L,numel(N_new),1),N_old];
    else 
        L = filterbanklength(L,a_new);
    end
    
end

%% retrieve the wavelets

if flags.do_complex
    [gout_positive,info_positive] = freqwavelet(winCell,L,scales,'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,'scal',scal(lowpass_number+1:M2), 'delay', delayvec(lowpass_number+1:M2+M), flags.norm);
    [gout_negative,info_negative] = freqwavelet(winCell,L,-flipud(scales),'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,'negative','scal',scal(M2+1:M2+M),'delay', delayvec(M2+1:M2+M), flags.norm);
    gout = [gout_positive,gout_negative];
    fields = fieldnames(info_positive);
    info = struct();
    for kk = 1:length(fields)
        if strcmp('a_natural',fields{kk})
            info.(fields{kk}) = [info_positive.(fields{kk});info_negative.(fields{kk})];
        else
            info.(fields{kk}) = [info_positive.(fields{kk}),info_negative.(fields{kk})];
        end
    end
elseif flags.do_analytic
    [gout,info] = freqwavelet(winCell,L,scales,'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,'analytic','scal',scal(lowpass_number+1:M2),'delay', delayvec(lowpass_number+1:M2),flags.norm);
else
    [gout,info] = freqwavelet(winCell,L,scales,'asfreqfilter','efsuppthr',kv.trunc_at,'basefc',0.1,'scal',scal(lowpass_number+1:M2),'delay', delayvec(lowpass_number+1:M2),flags.norm);
end
    
%% Generate lowpass filters if desired
if flags.do_single % Compute single lowpass from frequency response
    if numel(scales_sorted) < 4 % Single lowpass generation currently fails for fewer than four elements in scales
        error('%s: Lowpass generation requires at least 4 scales.',upper(mfilename));
    end
    lowpass_bandwidth = 0.2/scales_sorted(4); % Twice the center frequency of the fourth scale
    taper_ratio = 1-scales_sorted(4)/scales_sorted(2); % Plateau width is twice the center frequency of the second scale
    [glow,infolow] = wavelet_lowpass(gout,a,L,lowpass_bandwidth,taper_ratio,scal(1),flags);
    %[glow,infolow] = wavelet_lowpass(gout,a,L,lowpass_bandwidth,taper_ratio,scal(1)/3,flags);
    gout = [{glow},gout];
    fields = fieldnames(info);
    for kk = 1:length(fields) % Concatenate info and infolow
        %if strcmp('a_natural',fields{kk})
        %    info.(fields{kk}) = [infolow.(fields{kk});info.(fields{kk})];
        %else
            info.(fields{kk}) = [infolow.(fields{kk}),info.(fields{kk})];
        %end
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
           % if strcmp('a_natural',fields{kk})
           %     info.(fields{kk}) = [infolow.(fields{kk});info.(fields{kk});infohigh.(fields{kk})];
           % else
                info.(fields{kk}) = [infolow.(fields{kk}),info.(fields{kk}),infohigh.(fields{kk})];
           % end
        end
    else
        for kk = 1:length(fields) % Concatenate info and infolow
          %  if strcmp('a_natural',fields{kk})
          %      info.(fields{kk}) = [infolow.(fields{kk});info.(fields{kk})];
          %  else
                info.(fields{kk}) = [infolow.(fields{kk}),info.(fields{kk})];
          %  end
        end
    end
elseif flags.do_none % No lowpass, do nothing
    %Do nothing
else
    error('%s: This should not happen.',upper(mfilename));
end

% %% Adjust the downsampling rates in order to achieve 'redtar' [It may be possible to do this up front and not run filterbankscale]
 if ~isempty(kv.redtar)
    
    if size(a,2) == 2
        a_old = a(:,1)./a(:,2);
    else
        a_old = a;
    end
    
    if ~flags.do_real
        org_red = sum(1./a_old);
    elseif lowpass_at_zero
        org_red = 1./a_old(1) + sum(2./a_old(2:end));
    else
        org_red = sum(2./a_old);
    end
    
    a_new = floor(a*org_red/kv.redtar);
    scal_new = org_red/kv.redtar*ones(1,numel(gout));
    
    % Adjust function handle generated delays
    if isa(kv.delay,'function_handle')
        delayvec = zeros(M2,1);
        for kk = 1:M2
            delayvec(kk) = kv.delay(kk-1,a_new(kk));
        end
        if flags.do_complex
            % Replicate the delays, except at zero
            % frequency
            if lowpass_at_zero
                delayvec=[delayvec;flipud(delayvec(2:end))];
            else
                delayvec=[delayvec;flipud(delayvec)];
            end
        end
    end
    
    if ~flags.do_uniform
        N_old = ceil(L./a_old);
        N_new=ceil(L./a_new);
        a_new=[repmat(L,numel(N_new),1),N_old];
    else 
        L = filterbanklength(L,a_new);
    end
%     
     g_new = filterbankscale(gout,sqrt(scal_new));
     a = a_new;
     gout = g_new;
%     
%     if 0 % This is just for testing.
%         if size(a_new,2) == 2
%             a_test = a_new(:,1)./a_new(:,2);
%         else
%             a_test = a_new;
%         end
%         if ~flags.do_real
%             new_red = sum(1./a_test);
%         elseif lowpass_at_zero
%             new_red = 1./a_test(1) + sum(2./a_test(2:end));
%         else
%             new_red = sum(2./a_test);
%         end
%         % Compute and display redundancy for verification
%         fprintf('Original redundancy: %g \n', org_red);
%         fprintf('Target redundancy: %g \n', kv.redtar);
%         fprintf('Actual redundancy: %g \n', new_red);
%     end
 end

% Apply delays (might want to move this to freqwavelet instead)
for kk = 1:lowpass_number
   gout{kk}.delay = delayvec(kk);
end
info.startindex = lowpass_number + 1;%startindex of actual wavelets (tentative)
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

