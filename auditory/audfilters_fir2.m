 function [g,a,fc,L,info]=audfilters_fir2(fs,Ls,varargin)
%AUDFILTERS Generates filters equidistantly spaced on auditory frequency scales
%   Usage:  [g,a,fc,L]=audfilters(fs,Ls);
%           [g,a,fc,L]=audfilters(fs,Ls,...);
%
%   Input parameters:
%      fs    : Sampling rate (in Hz).
%      Ls    : Signal length.
%   Output parameters:
%      g     : Cell array of filters.
%      a     : Downsampling rate for each channel.
%      fc    : Center frequency of each channel.
%      L     : Next admissible length suitable for the generated filters.
%
%   `[g,a,fc,L]=audfilters(fs,Ls)` constructs a set of filters *g* that are
%   equidistantly spaced on a perceptual frequency scale (see |freqtoaud|) between
%   0 and the Nyquist frequency. The filter bandwidths are proportional to the 
%   critical bandwidth of the auditory filters |audfiltbw|. The filters are intended 
%   to work with signals with a sampling rate of *fs*. The signal length *Ls* is 
%   mandatory, since we need to avoid too narrow frequency windows.
%
%   By default the ERB scale is chosen but other frequency scales are
%   possible. Currently supported scales are 'erb', 'erb83', 'bark', 'mel'
%   and 'mel1000', and can be changed by passing the associated string as 
%   an optional parameter. See |freqtoaud| for more information on the
%   supported frequency scales.
%
%   By default, a Hann window shape is chosen as prototype frequency 
%   response for all filters. The prototype frequency response can be 
%   changed by passing any of the window types from |firwin| or |freqwin| 
%   as an optional parameter.
%
%   `[g,a,fc,L]=audfilters(fs,Ls,fmin,fmax)` constructs a set of filters 
%   between *fmin* and *fmax*. The filters are equidistantly spaced on the 
%   selected frequency scale. One additional filter will be positioned at 
%   the 0 and Nyquist frequencies each, so as to cover the full range of 
%   positive frequencies. 
%   The values of *fmin* and *fmax* can be instead specified using a 
%   key/value pair as::
%
%       [g,a,fc,L]=audfilters(fs,Ls,...,'fmin',fmin,'fmax',fmax)
%
%   Default values are *fmin=0* and *fmax=fs/2*. 
%
%   For more details on the construction of the filters, please see the
%   given references.
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
%   `[g,a]=audfilters(...,'regsampling')` constructs a non-uniform
%   filterbank with integer subsampling factors.
%
%   `[g,a]=audfilters(...,'uniform')` constructs a uniform filterbank
%   where the integer downsampling rate is the same for all the channels. This
%   results in most redundant representation which produces nice plots.
%
%   `[g,a]=audfilters(...,'fractional')` constructs a filterbank with
%   fractional downsampling rates *a*. 
%   This results in the least redundant system.
%
%   `[g,a]=audfilters(...,'fractionaluniform')` constructs a filterbank with
%   fractional downsampling rates *a*, which are uniform for all filters
%   except the "filling" low-pass and high-pass filters which can have different
%   fractional downsampling rates. This is useful when uniform subsampling
%   and low redundancy at the same time are desirable.
%
%   Additional parameters
%   ---------------------
%
%   `audfilters` accepts the following optional parameters:
%
%     'spacing',b        Specify the spacing between the filters, measured in
%                        scale units. Default value is *b=1* for the scales
%                        'erb', 'erb83' and 'bark'; the default is *b=100* for
%                        'mel' and 'mel1000'.
%
%     'bwmul',bwmul      Bandwidth of the filters relative to the bandwidth
%                        returned by |audfiltbw|. Default value is *bwmul=1* for 
%                        the scales 'erb', 'erb83' and 'bark'; the default is 
%                        *b=100* for 'mel' and 'mel1000'.
%
%     'redmul',redmul    Redundancy multiplier. Increasing the value of this
%                        will make the system more redundant by lowering the
%                        channel downsampling rates. Default
%                        value is *1*. If the value is less than one, the
%                        system may no longer be painless.
% 
%     'redtar',redtar    Target redundancy. The downsampling factors will be
%                        adjusted to achieve a redundancy as close as possible
%                        to 'redtar'.
%
%     'M',M              Specify the total number of filters between *fmin* and 
%                        *fmax*. If this parameter is specified, it overwrites the
%                        `'spacing'` parameter.
%
%     'symmetric'        Create filters that are symmetric around their centre
%                        frequency. This is the default.
%
%     'warped'           Create asymmetric filters that are symmetric on the
%                        auditory scale. 
%
%     'complex'          Construct a filterbank that covers the entire
%                        frequency range instead of just the positive 
%                        frequencies this allows the analysis of complex
%                        valued signals.
%
%     'nosubprec'        Disable subsample window positions.
%
%     'trunc_at'         When using a prototype defined in |freqwin|, a hard 
%                        thresholding of the filters at the specified threshold 
%                        value is performed to reduce their support size. 
%                        The default value is *trunc_at=10e-5*. When no 
%                        truncation is desired, *trunc_at=0* should be chosen.
%                        This value is ignored when a prototype shape from
%                        |firwin| was chosen.
%
%     'min_win',min_win  Minimum admissible window length (in samples).
%                        Default is *4*. This restrict the windows not
%                        to become too narrow when *L* is low.
%
%   Examples:
%   ---------
%
%   In the first example, we construct a highly redudant uniform
%   filterbank on the ERB scale and visualize the result:::
%
%     [f,fs]=greasy;  % Get the test signal
%     [g,a,fc,L]=audfilters(fs,length(f),'uniform','M',100);
%     c=filterbank(f,g,a);
%     plotfilterbank(c,a,fc,fs,90,'audtick');
%
%   In the second example, we construct a non-uniform filterbank with
%   fractional sampling that works for this particular signal length, and
%   test the reconstruction. The plot displays the response of the
%   filterbank to verify that the filters are well-behaved both on a
%   normal and an ERB-scale. The second plot shows frequency responses of
%   filters used for analysis (top) and synthesis (bottom). :::
%
%     [f,fs]=greasy;  % Get the test signal
%     L=length(f);
%     [g,a,fc]=audfilters(fs,L,'fractional');
%     c=filterbank(f,{'realdual',g},a);
%     r=2*real(ifilterbank(c,g,a));
%     norm(f-r)
%
%     % Plot the response
%     figure(1);
%     subplot(2,1,1);
%     R=filterbankresponse(g,a,L,fs,'real','plot');
%
%     subplot(2,1,2);
%     semiaudplot(linspace(0,fs/2,L/2+1),R(1:L/2+1));
%     ylabel('Magnitude');
%
%     % Plot frequency responses of individual filters
%     gd=filterbankrealdual(g,a,L);
%     figure(2);
%     subplot(2,1,1);
%     filterbankfreqz(gd,a,L,fs,'plot','linabs','posfreq');
%
%     subplot(2,1,2);
%     filterbankfreqz(g,a,L,fs,'plot','linabs','posfreq');
%
%
%   See also: filterbank, ufilterbank, ifilterbank, ceil23
%
%   References: ltfatnote027 nehobaprpide18

% Authors: Peter L. Søndergaard (original 'erbfilters' function)
% Modified by: Thibaud Necciari, Nicki Holighaus
% Comments updated by: Nicki Holighaus (09.05.22)

% Date: 16.12.16

complainif_notenoughargs(nargin,2,upper(mfilename));
complainif_notposscalar(fs,'fs',upper(mfilename));
complainif_notposint(Ls,'Ls',upper(mfilename));

firwinflags=getfield(arg_firwin,'flags','wintype');
freqwinflags=getfield(arg_freqwin,'flags','wintype');

definput.flags.wintype = [ firwinflags, freqwinflags];
definput.keyvals.M=[];
definput.keyvals.redmul=1;
definput.keyvals.min_win = 800;
definput.keyvals.bwmul=[];
definput.keyvals.spacing=[];
definput.keyvals.trunc_at=10^(-5);
definput.keyvals.fmin=0;
definput.keyvals.fmax=fs/2;
definput.keyvals.redtar=[];
definput.flags.subprec={'subprec','nosubprec'};
definput.flags.audscale={'erb','erb83','bark','mel','mel1000'};
definput.flags.warp     = {'symmetric','warped'};
definput.flags.real     = {'real','complex'};
definput.flags.sampling = {'regsampling','uniform','fractional',...
                           'fractionaluniform'};

[varargin,winCell] = arghelper_filterswinparser(definput.flags.wintype,varargin);

[flags,kv]=ltfatarghelper({'fmin','fmax'},definput,varargin);
if isempty(winCell), winCell = {flags.wintype}; end

switch flags.audscale
    case {'mel','mel1000'} % The mel scales are very fine, therefore default spacing is adjusted
        definput.keyvals.bwmul=100;
        definput.keyvals.spacing=100;
    otherwise
        definput.keyvals.bwmul=1;
        definput.keyvals.spacing=1;
end
[flags,kv]=ltfatarghelper({'fmin','fmax','redtar'},definput,varargin);

if flags.do_bark && (fs > 44100)
    error(['%s: Bark scale is not suitable for sampling rates higher than 44.1 kHz. ',...
    'Please choose another scale.'],upper(mfilename));
end 

if ~isscalar(kv.bwmul) || kv.bwmul <= 0
    error('%s: bwmul must be a positive scalar.',upper(mfilename));
end

if ~isscalar(kv.redmul) || kv.redmul <= 0
    error('%s: redmul must be a positive scalar.',upper(mfilename));
end

if ~isempty(kv.redtar)
    if ~isscalar(kv.redtar) || kv.redtar <= 0
        error('%s: redtar must be a positive scalar.',upper(mfilename));
    end
end

if kv.redtar <= 1
    warning('%s: redtar is very low; the resulting system might be unstable.',upper(mfilename));
end

if kv.fmax <= kv.fmin || kv.fmin < 0 || kv.fmax > fs/2
    error('%s: fmax must be bigger than fmin and in the range [0,fs/2].',upper(mfilename));
end

if kv.trunc_at > 1 || kv.trunc_at < 0
    error('%s: trunc_at must be in range [0,1].',upper(mfilename));
end

if ~isscalar(kv.min_win) || rem(kv.min_win,1) ~= 0 || kv.min_win < 1
    error('%s: min_win must be an integer bigger or equal to 1.',upper(mfilename));
end

if ~isempty(kv.M)
    complainif_notposint(kv.M,'M',upper(mfilename));
    kv.spacing = (freqtoaud(kv.fmax,flags.audscale) - freqtoaud(kv.fmin,flags.audscale))/(kv.M-1);
end

% Construct function handle for filter prototype and determine its ERB-type
% bandwidth
[filterfunc,winbw] = helper_filtergeneratorfunc_fir(...
                          flags.wintype,winCell,fs,kv.bwmul,kv.min_win,kv.trunc_at,...
                          flags.audscale,flags.do_subprec,flags.do_symmetric,flags.do_warped);

% Construct the AUD filterbank
fmin = max(kv.fmin,audtofreq(kv.spacing,flags.audscale));
fmax = min(kv.fmax,fs/2);

innerChanNum = floor((freqtoaud(fmax,flags.audscale)-freqtoaud(fmin,flags.audscale))/kv.spacing)+1;

fmax = audtofreq(freqtoaud(fmin,flags.audscale)+(innerChanNum-1)*kv.spacing,flags.audscale);

% Make sure that fmax < fs/2, and F_ERB(fmax) = F_ERB(fmin)+k/spacing, for
% some k.
count = 0;
while fmax >= fs/2
    count = count+1;
    fmax = audtofreq(freqtoaud(fmin,flags.audscale)+(innerChanNum-count-1)*kv.spacing,flags.audscale);    
end
innerChanNum = innerChanNum-count;

if fmax <= fmin || fmin > fs/4 || fmax < fs/4
    error(['%s: Bad combination of fs, fmax and fmin.'],upper(mfilename));
end

%%%%%%%%%%%%%%%%%% CHANGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Center frequencies are given as equidistantly spaced points on auditory 
% scale 
fc = audspace(fmin,fmax,innerChanNum,flags.audscale).';
fc = [fc;fs/2];
M2 = innerChanNum+1;

ind = (1:M2-1)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute the frequency support
% fsupp is measured in Hz 

fsupp=zeros(M2,1);
aprecise=zeros(M2,1);

fsupp(ind)=audfiltbw(fc(ind),flags.audscale)/winbw*kv.bwmul;

% % Generate lowpass filter parameters     
% fsupp(1) = 0; % Placeholder value 
% % Determine border of passband
% fps0 = audtofreq(freqtoaud(fc(2),flags.audscale)+3*kv.spacing,flags.audscale);% f_{p,s}^{-}
% % Determine lowpass width
% fsupp_temp1 = audfiltbw(fps0,flags.audscale)/winbw*kv.bwmul;
% % Determine bandwidth-adapted decimation factor
% aprecise(1) = max(fs./(2*max(fps0,0)+fsupp_temp1*kv.redmul),1);  

% Generate highpass filter parameters     
fsupp(end) = 0; % Placeholder value
% Determine border of passband
fps0 = audtofreq(freqtoaud(fc(end-1),flags.audscale)-3*kv.spacing,flags.audscale);% f_{p,s}^{+}
% Determine highpass width
fsupp_temp1 = audfiltbw(fps0,flags.audscale)/winbw*kv.bwmul;
% Determine bandwidth-adapted decimation factor
aprecise(end) = max(fs./(2*(fc(end)-min(fps0,fs/2))+fsupp_temp1*kv.redmul),1);

% Find suitable channel subsampling rates
aprecise(ind)=fs./fsupp(ind)/kv.redmul;
aprecise=aprecise(:);

if any(aprecise<1)
    error('%s: Invalid subsampling rates. Decrease redmul.',upper(mfilename))
end

%% Compute the downsampling rate
if flags.do_regsampling
    % Shrink "a" to the next composite number
    a=floor23(aprecise);

    % Determine the minimal transform length
    L=filterbanklength(Ls,a);

    % Heuristic trying to reduce lcm(a)
    while L>2*Ls && ~(all(a)==a(1))
        maxa = max(a);
        a(a==maxa) = 0;
        a(a==0) = max(a);
        L = filterbanklength(Ls,a);
    end

% Determine true decimation factors from "aprecise" according to chosen
% subsampling scheme (see help above)
elseif flags.do_fractional
    L = Ls;
    N=ceil(Ls./aprecise);
    a=[repmat(Ls,M2,1),N];
elseif flags.do_fractionaluniform
    L = Ls;
    N=ceil(Ls./min(aprecise));
    a= repmat([Ls,N],M2,1);
elseif flags.do_uniform
    a=floor(min(aprecise));
    L=filterbanklength(Ls,a);
    a = repmat(a,M2,1);
end

% Get an expanded "a" / Convert "a" to LTFAT 2-column fractional format
afull=comp_filterbank_a(a,M2,struct());

%% Compute the scaling of the filters
% Filters are scaled such that the energy of the subband coefficients
% remains approximately constant independent of the decimation factor
scal=sqrt(afull(:,1)./afull(:,2));

%% Construct the real or complex filterbank

if flags.do_real
    % Scale the first and last channels
    scal(1)=scal(1)/sqrt(2);
    scal(M2)=scal(M2)/sqrt(2);
else
    % Replicate the centre frequencies and sampling rates, except the first and
    % last
    a=[a;flipud(a(2:M2-1,:))];
    scal=[scal;flipud(scal(2:M2-1))];
    fc  =[fc; -flipud(fc(2:M2-1))];
    fsupp=[fsupp;flipud(fsupp(2:M2-1))];
    ind = [ind;numel(fc)+2-(M2-1:-1:2)'];
end


%%%%%%%%%%%%%%%%%%%%% NEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not allow longer support length than keyvals.min_win

% time supports (w/o high-pass)
tsupp = round(winbw./fsupp(1:end-1)*fs*10.64); 

% all channels that have to be shortened
ind_crit = find(tsupp(1:end-1) >= kv.min_win); 

 if ~isempty(ind_crit)
    % the center frequency for the last valid channel
    fc_crit = fc(ind_crit(end)+1);


    % the frequency step from the previous to the last valid fc 
    LP_step = (fc(ind_crit(end)+2) - fc_crit);

    % number of bands needed to get to 0
    LP_num = floor(fc_crit/LP_step); 

    % center frequencies of the lp and hp filters
    fc_high = fc(ind_crit(end)+1:end);
    fc_low = flip((fc_high(1)-(1:LP_num)*LP_step)'); 
    fc = [fc_low;fc_high];

    tsuppmax = tsupp(ind_crit(end)+1);
    %tsuppmax = kv.min_win;
    tsupp_low = ones(LP_num,1)*tsuppmax;
    tsupp_high = tsupp(ind_crit(end)+1:end);
    tsupp = [tsupp_low;tsupp_high];
end

M2 = length(fc);
ind = (1:M2-1)';

%LP_range = 0.1/scales_sorted(1); % freqmin
%LP_step = abs(0.1/scales_sorted(2)-LP_range); % b = (fmin+1 - fmin) wenn der letzte dann zu nah ist (<b/2) oder firekt bei 0 dann mit 1/sqrt(2)


%% Compute the filters
% This is actually much faster than the vectorized call.
g = cell(1,numel(fc));
for m=ind.'
    g{m}=filterfunc(tsupp(m),fc(m));
    g{m}.h = sqrt(a(m))*g{m}.h;
end

% if fc(1) < fc(2)/2
%     g{1}.h = g{1}.h/sqrt(1.5);
% end

a = a(1:M2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate lowpass filter
%g{1} = audlowpassfilter(g(1:M2),a(1:M2,:),fc(1:M2),fs,scal(1),kv,flags);

% Generate highpass filter
scal_hp = scal(end)*sqrt(2);
g{M2} = audhighpassfilter(g(1:M2),a(1:M2,:),fc(1:M2),fs,scal_hp,kv,flags);

% if g{M2}.H(0) == 0
%     g{M2-1}.h = g{M2-1}.h*1.075;
% end

winbwrat = winbw/0.754;
basebw = 1.875657;

info.fc  = 2*fc/fs;
info.tfr = @(L)(1/L)*1./((2*fsupp*winbwrat/fs)./basebw).^2;


    
function ghigh = audhighpassfilter(g,a,fc,fs,scal,kv,flags)

% Make a probe, compute the restricted filter bank response to check if hi-pass filter is needed
Lprobe = 10000;
FBresp0 = filterbankresponse(g(2:end-1),a(2:end-1,:),Lprobe,'real');
eps_thr = 1e-3;
ind_f1 = floor(fc(2)*Lprobe/fs);
ind_fK = floor(fc(end-1)*Lprobe/fs);
if ind_f1 == 0 ||...
    min(FBresp0(ind_fK:floor(Lprobe/2))) >= (1-eps_thr)*min(FBresp0(ind_f1:ind_fK))
%       Not required
    ghigh.H = @(L) 0;
    ghigh.foff = @(L) 0;
    ghigh.realonly = 0;
    ghigh.delay = 0;
    ghigh.fs = g{2}.fs;
else

    %     Compute the transition frequencies f_{p,s}^{+} and f_{p,e}^{+}
    % Determines the width of the plateau
    fps = audtofreq(freqtoaud(fc(end-1),flags.audscale)-3*kv.spacing,flags.audscale);
    % Determines the cosine transition frequency
    fpe = audtofreq(freqtoaud(fc(end-1),flags.audscale)-4*kv.spacing,flags.audscale);

    %     plateauWidth = 2*(fs/2-fps);
    fsupp_HP = 2*(fs/2-fpe);
    ratio = 2*(fps-fpe)/fsupp_HP;
    Lw = @(L) min(ceil(fsupp_HP*L/fs),L);

    PK = blfilter({'hann','taper',ratio},fsupp_HP,'fc',fs/2,'fs',fs,'inf','min_win',kv.min_win);
    temp_fbresp = @(L) filterbankresponse(g(2:end-1),a(2:end-1,:),L,'real');
    Hinv = @(L) sqrt(max(temp_fbresp(L))-temp_fbresp(L));

    %     Compute the final high-pass filter
    ghigh.H = @(L) fftshift(long2fir(fftshift(...
        filterbankfreqz(PK,a(1,:),L).*Hinv(L)),Lw(L)))*scal;

    ghigh.foff = @(L) ceil(L/2)-floor(Lw(L)/2)-1;
    ghigh.realonly = 0;
    ghigh.delay = 0;
    ghigh.fs = g{2}.fs;
end


