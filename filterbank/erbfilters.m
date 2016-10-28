function [g,a,fc,L]=erbfilters(fs,Ls,varargin)
%ERBFILTERS   ERB-spaced filters
%   Usage:  [g,a,fc]=erbfilters(fs,Ls);
%           [g,a,fc]=erbfilters(fs,Ls,...);
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
%   `[g,a,fc]=erbfilters(fs,Ls)` constructs a set of filters *g* that are
%   equidistantly spaced on the ERB-scale (see |freqtoerb|) with bandwidths
%   that are proportional to the width of the auditory filters
%   |audfiltbw|. The filters are intended to work with signals with a
%   sampling rate of *fs*. The signal length *Ls* is mandatory, since we
%   need to avoid too narrow frequency windows.
%
%   By default, a Hann window on the frequency side is choosen, but the
%   window can be changed by passing any of the window types from
%   |firwin| or |freqwin| as an optional parameter.
%   Run `getfield(getfield(arg_firwin,'flags'),'wintype')` or 
%   `getfield(getfield(arg_freqwin,'flags'),'wintype')` to get a cell
%   array of window types available.
%
%   The integer downsampling rates of the channels must all divide the
%   signal length, |filterbank| will only work for input signal lengths
%   being multiples of the least common multiple of the downsampling rates.
%   See the help of |filterbanklength|. 
%   The fractional downsampling rates restrict the filterbank to a single
%   length *L=Ls*.
%
%   `[g,a]=erbfilters(...,'regsampling')` constructs a non-uniform
%   filterbank with integer subsampling factors.
%
%   `[g,a]=erbfilters(...,'uniform')` constructs a uniform filterbank
%   where the integer downsampling rate is the same for all the channels. This
%   results in most redundant representation which produces nice plots.
%
%   `[g,a]=erbfilters(...,'fractional')` constructs a filterbank with
%   fractional downsampling rates *a*. 
%   This results in the least redundant system.
%
%   `[g,a]=erbfilters(...,'fractionaluniform')` constructs a filterbank with
%   fractional downsampling rates *a*, which are uniform for all filters
%   except the "filling" low-pass and high-pass filters can have different
%   fractional downsampling rates. This is usefull when uniform subsampling
%   and low redundancy at the same time are desirable.
%
%   `erbfilters` accepts the following optional parameters:
%
%     'spacing',b     Specify the spacing in ERBS between the
%                     filters. Default value is *b=1*.
%
%     'M',M           Specify the number of filters, *M*. If this
%                     parameter is specified, it overwrites the
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
%                     frequency. This is the default.
%
%     'warped'        Create asymmetric filters that are symmetric on the
%                     Erb-scale.
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
%   filterbank and visualize the result:::
%
%     [f,fs]=greasy;  % Get the test signal
%     [g,a,fc]=erbfilters(fs,length(f),'uniform','M',100);
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
%     [g,a,fc]=erbfilters(fs,L,'fractional');
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
%   References: ltfatnote027

% Authors: Peter L. Søndergaard, Zdenek Prusa, Nicki Holighaus

complainif_notenoughargs(nargin,2,upper(mfilename));
complainif_notposint(Ls,'Ls',upper(mfilename));

if ~isscalar(fs) || fs <= 0
    error('%s: fs must be positive scalar.',upper(mfilename))
end

firwinflags=getfield(getfield(arg_firwin,'flags'),'wintype');
freqwinflags=getfield(getfield(arg_freqwin,'flags'),'wintype');

definput.flags.wintype = [ firwinflags, freqwinflags];
definput.keyvals.M=[];
definput.keyvals.bwmul=1;
definput.keyvals.redmul=1;
definput.keyvals.min_win = 4;
definput.keyvals.spacing=1;
definput.keyvals.trunc_at=10^(-5);
definput.keyvals.flow=0;
definput.keyvals.fhigh=fs/2;
definput.flags.warp     = {'symmetric','warped'};
definput.flags.real     = {'real','complex'};
definput.flags.sampling = {'regsampling','uniform','fractional',...
                           'fractionaluniform'};

% Search for window given as cell array
candCellId = cellfun(@(vEl) iscell(vEl) && any(strcmpi(vEl{1},definput.flags.wintype)),varargin);

winCell = {};
% If there is such window, replace cell with function name so that 
% ltfatarghelper does not complain
if ~isempty(candCellId) && any(candCellId)
    winCell = varargin{candCellId(end)};
    varargin{candCellId} = [];
    varargin{end+1} = winCell{1};
end
                       
[flags,kv]=ltfatarghelper({},definput,varargin);
if isempty(winCell), winCell = {flags.wintype}; end

if ~isscalar(kv.bwmul) || kv.bwmul <= 0
    error('%s: bwmul must be a positive scalar.',upper(mfilename));
end

if ~isscalar(kv.redmul) || kv.redmul <= 0
    error('%s: bwmul must be a positive scalar.',upper(mfilename));
end

if kv.fhigh < kv.flow || kv.flow < 0 || kv.fhigh > fs/2
    error('%s: fhigh must be bigger than flow and in the range [0,fs/2].',upper(mfilename));
end

if kv.trunc_at > 1 || kv.trunc_at < 0
    error('%s: trunc_at must be in range [0,1].',upper(mfilename));
end

if ~isscalar(kv.min_win) || rem(kv.min_win,1) ~= 0 || kv.min_win < 1
    error('%s: min_win must be an integer bigger or equal to 1.',upper(mfilename));
end

if ~isempty(kv.M)
    complainif_notposint(kv.M,'M',upper(mfilename));
    kv.spacing = (freqtoerb(kv.fhigh) - freqtoerb(kv.flow))/kv.M;
end

probelen = 10000;

switch flags.wintype
    case firwinflags
        winbw=norm(firwin(flags.wintype,probelen)).^2/probelen;
        
        if flags.do_symmetric
            filterfunc = @(fsupp,fc,scal)... 
                         blfilter(winCell,fsupp,fc,'fs',fs,'scal',scal,...
                                  'inf','min_win',kv.min_win);
        else
            fsupp_erb=1/winbw*kv.bwmul;
            filterfunc = @(fsupp,fc,scal)...
                         warpedblfilter(winCell,fsupp_erb,fc,fs,...
                                        @freqtoerb,@erbtofreq,'scal',scal,'inf');
        end
        bwtruncmul = 1;
    case freqwinflags
        if flags.do_warped
            error('%s: TODO: Warping is not supported for windows from freqwin.',...
                upper(mfilename));
        end
        
        probebw = 0.01;
 
        % Determine where to truncate the window
        H = freqwin(winCell,probelen,probebw);
        winbw = norm(H).^2/(probebw*probelen/2);
        bwrelheight = 10^(-3/10);
        
        if kv.trunc_at <= eps
            bwtruncmul = inf;
        else
            try
                bwtruncmul = winwidthatheight(abs(H),kv.trunc_at)/winwidthatheight(abs(H),bwrelheight);
            catch
                bwtruncmul = inf;
            end
        end
        
        filterfunc = @(fsupp,fc,scal)...
                     freqfilter(winCell, fsupp, fc,'fs',fs,'scal',scal,...
                                'inf','min_win',kv.min_win,...
                                'bwtruncmul',bwtruncmul);        
end

% Construct the Erb filterbank
flow = max(0,kv.flow);
fhigh = min(kv.fhigh,fs/2);

innerChanNum = ceil((freqtoerb(fhigh)-freqtoerb(flow))/kv.spacing)+1;

fhigh = erbtofreq(freqtoerb(flow)+(innerChanNum-1)*kv.spacing);

% Make sure that fhigh <= fs/2, and F_ERB(fhigh) = F_ERB(flow)+k/spacing, for
% some k.
count = 1;
while fhigh > fs/2
    count = count+1;
    fhigh = erbtofreq(freqtoerb(flow)+(innerChanNum-count)*kv.spacing);
end

if flags.do_real
    if isempty(kv.M)
        M2=innerChanNum;
        M=M2;
    else
        M=kv.M;
        M2=M;
    end;
else
    if isempty(kv.M)
        M2=innerChanNum;
        M=2*(M2-2);
    else
        M=kv.M;
        if rem(M,2)>0
            error(['%s: M must be even for full frequency range ' ...
                   'filterbanks.',upper(mfilename)]);
        end;
        M2=M/2+1;
    end;

end;

fc=erbspace(flow,fhigh,M2).';
if flow > 0
    fc = [0;fc];
    M2 = M2+1;
    M=2*(M2-2);
end

if fhigh < fs/2
    fc = [fc;fs/2];
    M2 = M2+1;
    M=2*(M2-2);
end

ind = (1:M2)';
%% Compute the frequency support
% fsupp is measured in Hz 

fsupp=zeros(M2,1);
aprecise=zeros(M2,1);
if flow > 0
    ind = ind(2:end);
end
if fhigh < fs/2
    ind = ind(1:end-1);
end

if flags.do_symmetric
    fsupp(ind)=audfiltbw(fc(ind))/winbw*kv.bwmul;
    
    if flow > 0 
        fsupp(1) = 0;%(2*fc(2)+4*audfiltbw(fc(2))/winbw);
        fc_temp = max(erbtofreq(freqtoerb(fc(2))-kv.spacing),0);
        fsupp_temp = audfiltbw(fc_temp)/winbw*kv.bwmul;
        aprecise(1) = max(fs./(2*fc_temp+fsupp_temp*kv.redmul),1);%/kv.redmul,1);        
    end
    if fhigh < fs/2        
        fsupp(end) = 0;%(2*(fc(end)-fc(end-1))+4*audfiltbw(fc(end-1))/winbw);
        fc_temp = min(erbtofreq(freqtoerb(fc(end-1))+kv.spacing),fs/2);
        fsupp_temp = audfiltbw(fc_temp)/winbw*kv.bwmul;
        aprecise(end) = max(fs./(2*(fc(end)-fc_temp)+fsupp_temp*kv.redmul),1);%/kv.redmul,1);       
    end    
else    
    % fsupp_erb is measured in Erbs
    % The scaling is incorrect, it does not account for the warping
    fsupp_erb=1/winbw*kv.bwmul;

    % Convert fsupp into the correct widths in Hz, necessary to compute
    % "a" in the next if-statement
    fsupp(ind)=erbtofreq(freqtoerb(fc(ind))+fsupp_erb/2)-erbtofreq(freqtoerb(fc(ind))-fsupp_erb/2);
    
    if flow > 0 
        fsupp(1) = 0;%(2*fc(2)+4*audfiltbw(fc(2))/winbw);
        fc_temp = max(erbtofreq(freqtoerb(fc(2))-kv.spacing),0);
        fsupp_temp=2*(erbtofreq(freqtoerb(fc_temp)+fsupp_erb/2)-fc_temp);        
        aprecise(1) = max(fs./(2*fc_temp+fsupp_temp*kv.redmul),1);%/kv.redmul,1);        
    end
    if fhigh < fs/2        
        fsupp(end) = 0;%(2*(fc(end)-fc(end-1))+4*audfiltbw(fc(end-1))/winbw);
        fc_temp = min(erbtofreq(freqtoerb(fc(end-1))+kv.spacing),fs/2);
        fsupp_temp=2*(fc_temp-erbtofreq(freqtoerb(fc_temp)-fsupp_erb/2));        
        aprecise(end) = max(fs./(2*(fc(end)-fc_temp)+fsupp_temp*kv.redmul),1);%/kv.redmul,1);       
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
aprecise(ind)=fs./fsupp(ind)/kv.redmul;
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
    %if flags.do_symmetric
        fsupp=[fsupp;flipud(fsupp(2:M2-1))];
    %end;
    ind = [ind;numel(fc)+2-(M2-1:-1:2)'];
end;


%% Compute the filters
% This is actually much faster than the vectorized call.
g = cell(1,numel(fc));
for m=ind.'
    g{m}=filterfunc(fsupp(m),fc(m),scal(m));
end

if flow > 0
    g{1} = erblowpassfilter(filterfunc,g,a,flow,fs,winbw,scal(1),M2,bwtruncmul,kv,flags);
end

if fhigh < fs/2
    g{M2} = erbhighpassfilter(filterfunc,g,a,fhigh,fs,winbw,scal(M2),M2,bwtruncmul,kv,flags);
end

function glow = erblowpassfilter(filterfunc,g,a,flow,fs,winbw,scal,M2,bwtruncmul,kv,flags)
    next_fc = max(erbtofreq(freqtoerb(flow)-kv.spacing),0);
    temp_fc = [];    
% Temporary center frequencies and bandwidths for lowpass -----------------
% (only makes sense for symmetric filters)
%     ERB = @(f) 9.265.*log(1+f./228.8455);
%     iERB = @(E) 228.8455.*(exp(E./9.265)-1);
%     bwERB = @(f) 24.7+f./9.265;
%     while next_fc >= iERB(-4)
%         temp_fc(end+1) = next_fc;
%         next_fc = iERB(ERB(next_fc)-kv.spacing);
%     end
%     
%     temp_bw=bwERB(temp_fc)/winbw*kv.bwmul;
% [end] -------------------------------------------------------------------

% Temporary center frequencies and bandwidths for lowpass - Simplified ----
    while next_fc >= erbtofreq(-4)
        temp_fc(end+1) = next_fc;
        next_fc = erbtofreq(freqtoerb(next_fc)-kv.spacing);
    end
    
    if flags.do_symmetric
        temp_bw=audfiltbw(abs(temp_fc))/winbw*kv.bwmul;
        plateauWidth = max(2*temp_fc(1),0);
        Lw = @(L) min(ceil((plateauWidth+temp_bw(1)*bwtruncmul)*L/fs),L);
    else
        fsupp_erb=1/winbw*kv.bwmul;
        temp_bw = erbtofreq(freqtoerb(temp_fc)+fsupp_erb/2)-erbtofreq(freqtoerb(temp_fc)-fsupp_erb/2);
        Lw = @(L) min(ceil(2*erbtofreq(freqtoerb(temp_fc(1))+fsupp_erb/2)*L/fs),L);
    end
% Simplified [end]---------------------------------------------------------
    
    temp_g = cell(1,numel(temp_fc));
    for m=1:numel(temp_g)
        temp_g{m}=filterfunc(temp_bw(m),temp_fc(m),1);%scal.^2);
    end    
    

    %Lw = @(L) L;
    
    temp_fun = @(L) filterbankresponse(temp_g,1,L);
    if 1%flags.do_symmetric
        temp_fbresp = @(L) filterbankresponse(g(2:M2-1),a(2:M2-1,:),L);
        glow.H = @(L) fftshift(long2fir(...
                    sqrt(postpad(ones(ceil(L/2),1),L).*temp_fun(L) + ...
                    flipud(postpad(ones(L-ceil(L/2),1),L).*circshift(temp_fun(L),-1)) - ...
                    (flipud(postpad(ones(round(L/4)-1,1),L)).*temp_fbresp(L) + ...
                    postpad(ones(round(L/4),1),L).*flipud(temp_fbresp(L)))),...
                    Lw(L)))*scal;
    else
        glow.H = @(L) fftshift(long2fir(...
                    sqrt(postpad(ones(ceil(L/2),1),L).*temp_fun(L) + ...
                    flipud(postpad(ones(L-ceil(L/2),1),L).*circshift(temp_fun(L),-1))),...
                    Lw(L)))*scal;
    end
                    
    glow.foff = @(L) -floor(Lw(L)/2);
    glow.realonly = 0;
    glow.delay = 0;
    glow.fs = g{2}.fs;
    %glow.H = @(L) glow.H(L)*scal;%*sqrt(2);
    
function ghigh = erbhighpassfilter(filterfunc,g,a,fhigh,fs,winbw,scal,M2,bwtruncmul,kv,flags)
    next_fc = min(erbtofreq(freqtoerb(fhigh)+kv.spacing),fs/2);
    temp_fc = [];    
    while next_fc <= erbtofreq(freqtoerb(fs/2)+4)
        temp_fc(end+1) = next_fc;
        next_fc = erbtofreq(freqtoerb(next_fc)+kv.spacing);
    end
    
    if flags.do_symmetric
        temp_bw=audfiltbw(abs(temp_fc))/winbw*kv.bwmul;
        plateauWidth = max(2*(fs/2-temp_fc(1)),0);
        Lw = @(L) min(ceil((plateauWidth+temp_bw(1)*bwtruncmul)*L/fs),L);
    else
        fsupp_erb=1/winbw*kv.bwmul;
        temp_bw = erbtofreq(freqtoerb(temp_fc)+fsupp_erb/2)-erbtofreq(freqtoerb(temp_fc)-fsupp_erb/2);
        Lw = @(L) min(ceil(2*(fs/2-erbtofreq(freqtoerb(temp_fc(1))-fsupp_erb/2))*L/fs),L);
    end
    
    temp_g = cell(1,numel(temp_fc));
    for m=1:numel(temp_g)
        temp_g{m}=filterfunc(temp_bw(m),temp_fc(m),1);%scal.^2);
    end
    
    %plateauWidth = max(2*(fs/2-temp_fc(1)),0);
    %Lw = @(L) min(ceil((plateauWidth+temp_bw(1)*bwtruncmul)*L/fs),L);
    
    temp_fun = @(L) filterbankresponse(temp_g,1,L);
    if 1%flags.do_symmetric
        temp_fbresp = @(L) filterbankresponse(g(2:M2-1),a(2:M2-1,:),L);
        ghigh.H = @(L) fftshift(long2fir(...
                         fftshift(sqrt(postpad(ones(ceil(L/2),1),L).*temp_fun(L) + ...
                         flipud(postpad(ones(L-ceil(L/2),1),L).*circshift(temp_fun(L),-1)) - ...
                         (postpad([zeros(ceil(L/2),1);ones(round(L/4),1)],L).*temp_fbresp(L) + ...
                         flipud(postpad([zeros(L-ceil(L/2),1);ones(round(L/4),1)],L).*temp_fbresp(L))))),...
                         Lw(L)))*scal;
    else
        ghigh.H = @(L) fftshift(long2fir(...
                         fftshift(sqrt(postpad(ones(ceil(L/2),1),L).*temp_fun(L) + ...    
                         flipud(postpad(ones(L-ceil(L/2),1),L).*circshift(temp_fun(L),-1)))),...
                         Lw(L)))*scal;
    end
    
    ghigh.foff = @(L) ceil(L/2)-floor(Lw(L)/2)-1;
    ghigh.realonly = 0;
    ghigh.delay = 0;
    ghigh.fs = g{2}.fs;
%    ghigh.H = @(L) ghigh.H(L)*scal;%*sqrt(2);


function width = winwidthatheight(gnum,atheight)

width = zeros(size(atheight));
for ii=1:numel(atheight)
    gl = numel(gnum);
    gmax = max(gnum);
    frac=  1/atheight(ii);
    fracofmax = gmax/frac;
    
    
    ind =find(gnum(1:floor(gl/2)+1)==fracofmax,1,'first');
    if isempty(ind)
        %There is no sample exactly half of the height
        ind1 = find(gnum(1:floor(gl/2)+1)>fracofmax,1,'last');
        ind2 = find(gnum(1:floor(gl/2)+1)<fracofmax,1,'first');
%         if isempty(ind2)
%            width(ii) = gl;
%         else 
           rest = 1-(fracofmax-gnum(ind2))/(gnum(ind1)-gnum(ind2));
           width(ii) = 2*(ind1+rest-1);
%        end        
    else
        width(ii) = 2*(ind-1);
    end
end

