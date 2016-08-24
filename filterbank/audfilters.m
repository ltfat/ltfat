function [g,a,fc,L]=audfilters(fs,Ls,varargin)
% AUDFILTERS generates AUD-spaced filters
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
%   0 and the Nyquist frequency and with bandwidths that are proportional to the critical 
%   bandwidth of the auditory filters |audfiltbw|. The filters are intended to work
%   with signals with a sampling rate of *fs*. The signal length *Ls* is mandatory,
%   since we need to avoid too narrow frequency windows.
% 
%   By default the ERB scale is chosen but other frequency scales are
%   possible. See 'freqtoaud' for all available options. The most scales
%   are 'erb', 'bark', and 'mel'.
%
%   By default, a Hann window on the frequency side is chosen, but the
%   window can be changed by passing any of the window types from
%   |firwin| as an optional parameter.
%   Run `getfield(getfield(arg_firwin,'flags'),'wintype')` to get a cell
%   array of window types available.
%
%   The integer downsampling rates of the channels must all divide the
%   signal length, |filterbank| will only work for input signal lengths
%   being multiples of the least common multiple of the downsampling rates.
%   See the help of |filterbanklength|. 
%   The fractional downsampling rates restrict the filterbank to a single
%   length *L=Ls*.
% 
%   `[g,a,fc,L]=audfilters(fs,Ls,flow,figh)` constructs a set of filters that are
%   equidistantly spaced between flow and fhigh. In that case two
%   additional filters will be positioned at the 0 and Nyquist frequencies
%   so as to cover the full spectrum. The values of *flow* and *fhigh* can
%   be instead specified using a key/value pair as::
%
%       [g,a,fc,L]=audfilters(fs,Ls,...,'flow',flow,'fhigh',figh)
%
%   `[g,a,fc,L]=audfilters(...,'regsampling')` constructs a non-uniform
%   filterbank with integer subsampling factors.
%
%   `[g,a,fc,L]=audfilters(...,'uniform')` constructs a uniform filterbank
%   where the integer downsampling rate is the same for all the channels. This
%   results in most redundant representation which produces nice plots.
%
%   `[g,a,fc,L]=audfilters(...,'fractional')` constructs a filterbank with
%   fractional downsampling rates *a*. 
%   This results in the least redundant system.
%
%   `[g,a,fc,L]=audfilters(...,'fractionaluniform')` constructs a filterbank with
%   fractional downsampling rates *a*, which are uniform for all filters
%   except the "filling" low-pass and high-pass filters can have different
%   fractional downsampling rates. This is usefull when uniform subsampling
%   and low redundancy at the same time are desirable.
%
%   `audfilters` accepts the following optional parameters:
%
%     'spacing',b     Specify the spacing in Ecritical bandwidth (ERB or Bark
%                     depending on the scale) between the filters. Default value
%                     is *b=1*.
%
%     'M',M           Specify the total number of filters between 'flow' and 'fhigh', *M*.
%                     If this parameter is specified, it overwrites the
%                     `'spacing'` parameter.
%
%     'redmul',redmul  Redundancy multiplier. Increasing the value of this
%                      will make the system more redundant by lowering the
%                      channel downsampling rates. It is only used if the
%                      filterbank is a non-uniform filterbank. Default
%                      value is *1*. If the value is less than one, the
%                      system may no longer be painless.
%
%     'symmetric'     Create filters that are symmetric around their centre
%                     frequency. This is the default.'sqrtsquare','sqrtrect'
%
%     'warped'        Create asymmetric filters that are asymmetric on the
%                     ERB scale. The warping does not work with other
%                     scales yet.
%
%     'complex'       Construct a filterbank that covers the entire
%                     frequency range.
%
%
%     'bwmul',bwmul   Bandwidth of the filters relative to the bandwidth
%                     returned by |audfiltbw|. Default is $bwmul=1$.
%
%     'min_win',min_win     Minimum admissible window length (in samples).
%                           Default is *4*. This restrict the windows not
%                           to become too narrow when *L* is low.
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
%     [g,a,fc,L]=audfilters(fs,L,'fractional');
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
%   See also: erbfilters, filterbank, ufilterbank, ifilterbank, ceil23
%
%   References: ltfatnote027

% Authors: Peter L. Søndergaard (original 'erbfilters' function)
% Modified by: Thibaud Necciari
% Date: 30.09.15

%% ------ Checking of input parameters ---------

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

complainif_notposint(fs,'fs');
complainif_notposint(Ls,'Ls');

definput.import={'freqtoaud','firwin'};
definput.keyvals.M=[];
definput.keyvals.bwmul=1;
definput.keyvals.redmul=1;
definput.keyvals.min_win = 4;
definput.keyvals.spacing=1;
definput.keyvals.flow=0;
definput.keyvals.fhigh=fs/2;
definput.flags.warp     = {'symmetric','warped'};
definput.flags.real     = {'real','complex'};
definput.flags.sampling = {'regsampling','uniform','fractional',...
                           'fractionaluniform'};

[flags,kv,flow,fhigh]=ltfatarghelper({'flow','fhigh'},definput,varargin);


% Default parameters.

if ~isnumeric(flow) || ~isscalar(flow) 
  error('%s: flow must be a scalar.',upper(mfilename));
end;

if ~isnumeric(fhigh) || ~isscalar(fhigh) 
  error('%s: fhigh must be a scalar.',upper(mfilename));
end;

if flow>fhigh
  error('%s: flow must be less than or equal to fhigh.',upper(mfilename));
end;

if fhigh>fs/2
  error('%s: fhigh must be smaller or equal to the Nyquist frequency.',upper(mfilename));
end;


% Get the bandwidth of the chosen window by doing a probe
winbw=norm(firwin(flags.wintype,1000)).^2/1000;% This is the ERB at 1000 Hz

% Construct the AUD filterbank
if flags.do_real
    if isempty(kv.M)
        M2=ceil((freqtoaud(fhigh,flags.audscale)-freqtoaud(flow,flags.audscale))/kv.spacing)+1;
%         M=M2;
    else
        M=kv.M;
        M2=M;
    end;
else
    if isempty(kv.M)
        M2=ceil(freqtoaud(fs/2,flags.audscale)/kv.spacing)+1;
        M=2*(M2-1);
    else
        M=kv.M;
        if rem(M,2)>0
            error(['%s: M must be even for full frequency range ' ...
                   'filterbanks.',upper(mfilename)]);
        end;
        M2=M/2+1;
    end;
end;
% The spacing will be used below to compute "virtual" filters...
spac = (freqtoaud(fhigh,flags.audscale)-freqtoaud(flow,flags.audscale))/M2;

% Compute center frequencies on the perceptual scale
fc=audspace(flow,fhigh,M2,flags.audscale).';

% Compute bandwidths on the corresponding auditory scale
if flags.do_erb || flags.do_erb83
    cb = audfiltbw(fc,'erb');
elseif flags.do_bark
    cb = audfiltbw(fc,'bark');
% else
%     No auditory bandwidth concept applies to other scales, the frequency support
%     will then be computed below so as to achieve approx. 50% overlap between channels
end

%% Compute the frequency support
if flags.do_symmetric
    % fsupp is measured in Hz
    if flags.do_erb || flags.do_erb83 || flags.do_bark
        fsupp=round(cb/winbw*kv.bwmul);
%         Add additional filters if required
        if flow > 0 
            fc = [0;fc];
            fsupp = [2*flow;fsupp];
            M2 = M2 + 1;
        end
        if fhigh < fs/2
            fc = [fc;fs/2];
            fsupp = [fsupp;2*(fs/2-fhigh)];
            M2 = M2 + 1;
        end
    else
%         First check for frequency range and add additional filters if required
        if flow > 0 
            fc = [0;fc];
            M2 = M2 + 1;
        end
        if fhigh < fs/2
            fc = [fc;fs/2];
            M2 = M2 + 1;
        end
        fsupp = zeros(M2,1);
        fsupp(1) = 2*fc(find(fc,1));
        for k = 2:M2-1
            fsupp(k) = (fc(k+1)-fc(k-1));
        end
        if fhigh < fs/2
            % Correct the bandwidth of the last filter positioned at fhigh
            % Frequency of the (virtual) filter positioned at fhigh+spac
            f1 = min(fs/2,audtofreq(freqtoaud(fhigh,flags.audscale)+spac,flags.audscale));
            fsupp(M2-1)= f1-fc(end-1);
        end
        fsupp(M2) = 2*(fs/2-fc(M2-1));
    end
    
else
    if flags.do_erb || flags.do_erb83
        % fsupp_erb is measured in Erbs
        % The scaling is incorrect, it does not account for the warping
        fsupp_erb=1/winbw*kv.bwmul;

        % Convert fsupp into the correct widths in Hz, necessary to compute
        % "a" in the next if-statement
        fsupp=audtofreq(freqtoaud(fc,flags.audscale)+fsupp_erb/2,flags.audscale)...
            -audtofreq(freqtoaud(fc,flags.audscale)-fsupp_erb/2,flags.audscale);
    else
%         [FIXME] WARPING ON OTHER SCALES?
        error('%s: Warped asymmetric filters can only be achieved on the ERB scale.',upper(mfilename));
    end
end;

% Do not allow lower bandwidth than keyvals.min_win
fsuppmin = kv.min_win/Ls*fs;
for ii = 1:numel(fsupp)
    if fsupp(ii) < fsuppmin;
        fsupp(ii) = fsuppmin;
    end
end

% Find suitable channel subsampling rates
% If flow > 0 and fhigh < fs/2, then do not apply redmul to channels 1 and M2
% as it produces unecessarily badly conditioned frames!
aprecise=fs./fsupp/kv.redmul;
if flow > 0
    aprecise(1)=aprecise(1)*kv.redmul;
end
if fhigh < fs/2
    aprecise(end)=aprecise(end)*kv.redmul;
end
aprecise=aprecise(:);

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
end;

% Get an expanded "a"
afull=comp_filterbank_a(a,M2,struct());

%% Compute the scaling of the filters
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
    if flags.do_symmetric
        fsupp=[fsupp;flipud(fsupp(2:M2-1))];
    end;

end;


%% Compute the filters
if flags.do_symmetric
    % This is actually much faster than the vectorized call.
    g = cell(1,numel(fc));
    if flow > 0
        % Frequency of the filter that would be positioned left below flow on the perceptual scale
        f0 = max(0,audtofreq(freqtoaud(flow,flags.audscale)-spac,flags.audscale));
        if flags.do_erb || flags.do_erb83 || flags.do_bark
            cb0 = audfiltbw(f0,flags.audscale);
            fsupp0 = round(cb0/winbw*kv.bwmul);
        else
            f00 = audtofreq(freqtoaud(flow,flags.audscale)-2*spac,flags.audscale);
            fsupp0= fc(2)-f00;
        end
        % Compute the filter amd middle-pad it to cover the full frequency
        % range not covered by the filter bank
        g{1} = blfilter({flags.wintype,'taper',fsupp0/(2*fsupp(1))},fsupp(1),fc(1),'fs',fs,'scal',scal(1),...
                   'inf','min_win',kv.min_win);   
    else
        g{1}=blfilter(flags.wintype,fsupp(1),fc(1),'fs',fs,'scal',scal(1),...
                   'inf','min_win',kv.min_win);
    end
    for m=2:numel(g)-1
        g{m}=blfilter(flags.wintype,fsupp(m),fc(m),'fs',fs,'scal',scal(m),...
                   'inf','min_win',kv.min_win);
    end
    if fhigh < fs/2
        % Do the same on the high frequency side, frequency of the filter right next
%         above fhigh on the perceptual scale
        f0 = min(fs/2,audtofreq(freqtoaud(fhigh,flags.audscale)+spac,flags.audscale));
        if flags.do_erb || flags.do_erb83 || flags.do_bark
            cb0 = audfiltbw(f0,flags.audscale);
            fsupp0 = round(cb0/winbw*kv.bwmul);
        else
            f1 = audtofreq(freqtoaud(fhigh,flags.audscale)+2*spac,flags.audscale);
            fsupp0= f1-fc(end-1);
        end
        % Compute the filter and middle-pad it
        g{end} = blfilter({flags.wintype,'taper',fsupp0/(2*fsupp(end))},fsupp(end),fc(end),'fs',fs,'scal',scal(end),...
                   'inf','min_win',kv.min_win);
    else
        g{end}=blfilter(flags.wintype,fsupp(end),fc(end),'fs',fs,'scal',scal(end),...
                   'inf','min_win',kv.min_win);
    end
else
    g = cell(1,numel(fc));
    for m=1:numel(g)
        g{m}=warpedblfilter(flags.wintype,fsupp_erb,fc(m),fs,@freqtoaud,@audtofreq, ...
                     'scal',scal(m),'inf');
    end
end;

end

