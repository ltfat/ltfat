function [g,a,fc,L]=erbfilters(fs,varargin)
%ERBFILTERS   ERB-spaced filters
%   Usage:  [g,a,fc]=erbfilters(fs);
%           [g,a,fc]=erbfilters(fs,...);
%
%   Input parameters:
%      fs    : Sampling rate (in Hz).
%   Output parameters:
%      g     : Cell array of filters.
%      a     : Downsampling rate for each channel.
%      fc    : Center frequency of each channel.
% 
%   `[g,a,fc]=erbfilters(fs)` constructs a set of filters *g* that are
%   equidistantly spaced on the ERB-scale (see |freqtoerb|) with bandwidths
%   that are proportional to the width of the auditory filters
%   |audfiltbw|. The filters are intended to work with signals with a
%   sampling rate of *fs*.
%
%   By default, a Hann window on the frequency side is choosen, but the
%   window can be changed by passing any of the window types from
%   |firwin| as an optional parameter.
%
%   Because the downsampling rates of the channels must all divide the
%   signal length, |filterbank| will only work for multiples of the
%   least common multiple of the downsampling rates. See the help of
%   |filterbanklength|.
%
%   `[g,a,fc]=erbfilters(fs,L,'fractional')` constructs a filterbank with
%   fractional downsampling rates *a*. The rates are constructed such
%   that the filterbank can handle signal length that are multiples of
%   *L*, so the benefit of the fractional downsampling is that you get to
%   choose the value returned by |filterbanklength|.
%
%   `[g,a,fc]=erbfilters(fs,'uniform')` constructs a uniform filterbank
%   where the downsampling rate is the same for all channels.
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
%     'nonuniform'    Construct a non-uniform filterbank. This is the
%                     default.
%    
%     'uniform'       Construct a uniform filterbank.
%
%     'symmetric'     Create filters that are symmetric around their centre
%                     frequency. This is the default.
%
%     'warped'        Create asymmetric filters that are symmetric on the
%                     Erb-scale.
%
%     'real'          Construct a filterbank that works for real-valued
%                     signals only (the filters cover only the positive
%                     frequencies). This is the default.
%
%     'complex'       Construct a filterbank that covers the entire
%                     frequency range.
%
%     'regsampling'   Choose the downsampling rates to be products of 2
%                     and 3 (see |floor23| and |ceil23|). This is the
%                     default.
%
%     'fractional'    Use fractional downsampling. If this flag is
%                     specified, you must also specify the `'L'` parameter.
%
%     'L',L           Specify a transform length for which the fractional
%                     sampling rates must match up.
%
%     'bwmul',bwmul   Bandwidth of the filters relative to the bandwidth
%                     returned by |audfiltbw|. Default is $bwmul=1$.
%
%   Examples:
%   ---------
%
%   In the first example, we construct a highly redudant uniform
%   filterbank and visualize the result:::
%
%     [f,fs]=greasy;  % Get the test signal
%     [g,a,fc]=erbfilters(fs,'uniform','M',100);
%     c=filterbank(f,g,a);
%     plotfilterbank(c,a,fc,fs,90,'audtick');
%
%   In the second example, we construct a non-uniform filterbank with
%   fractional sampling that works for this particular signal length, and
%   test the reconstruction. The plot displays the response of the
%   filterbank to verify that the filters are well-behaved both on a
%   normal and an ERB-scale:::
%
%     [f,fs]=greasy;  % Get the test signal
%     L=length(f);
%     [g,a,fc]=erbfilters(fs,'fractional','L',L);
%     c=filterbank(f,{'realdual',g},a);
%     r=2*real(ifilterbank(c,g,a));
%     norm(f-r)
%
%     % Plot the response
%     subplot(2,1,1);     
%     R=filterbankresponse(g,a,L,fs,'real','plot');
%
%     subplot(2,1,2);
%     semiaudplot(linspace(0,fs/2,L/2+1),R(1:L/2+1));
%     ylabel('Magnitude');
%
%   See also: filterbank, ufilterbank, ifilterbank, ceil23
%
%   References: ltfatnote027

% Authors: Peter L. Søndergaard

definput.import = {'firwin'};

definput.keyvals.L=[];
definput.keyvals.M=[];
definput.keyvals.bwmul=1;
definput.keyvals.redmul=1;
definput.keyvals.spacing=1;

definput.flags.warp     = {'symmetric','warped'};
definput.flags.uniform  = {'nonuniform','uniform'};
definput.flags.real     = {'real','complex'};
definput.flags.sampling = {'regsampling','fractional'};

[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

% Get the bandwidth of the choosen window by doing a probe
winbw=norm(firwin(flags.wintype,1000)).^2/1000;

% Construct the Erb filterbank


if flags.do_real
    if isempty(kv.M)
        M2=ceil(freqtoerb(fs/2)/kv.spacing)+1;
        M=M2;
    else
        M=kv.M;
        M2=M;
    end;
else
    if isempty(kv.M)
        M2=ceil(freqtoerb(fs/2)/kv.spacing)+1;
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

fc=erbspace(0,fs/2,M2).';

    
%% Compute the frequency support
if flags.do_symmetric
    % fsupp is measured in Hz
    fsupp=round(audfiltbw(fc)/winbw*kv.bwmul);
else
    % fsupp_erb is measured in Erbs
    % The scaling is incorrect, it does not account for the warping
    fsupp_erb=1/winbw*kv.bwmul;
    
    % Convert fsupp into the correct widths in Hz, necessary to compute
    % "a" in the next if-statement
    fsupp=erbtofreq(freqtoerb(fc)+fsupp_erb/2)-erbtofreq(freqtoerb(fc)-fsupp_erb/2);
    
end;

%% Compute the downsampling rate
if flags.do_nonuniform
    % Find suitable channel subsampling rates
    aprecise=fs./fsupp/kv.redmul; 
    aprecise=aprecise(:);
    
    if flags.do_fractional
       if isempty(L)
          error('%s: ''fractional'' flag refuires L to be specified.',upper(mfilename));
       end
        N=ceil(L./aprecise);
        a=[repmat(L,M2,1),N];                
        
    else
        a=floor23(aprecise); % Shrink "a" to the next composite number
        
        % Determine the minimal transform length
        L=filterbanklength(1,a);
    end;
    
else
    % Totally heuristic, no justification for this choice
    a=4;    
end;

% Get an expanded "a"
afull=comp_filterbank_a(a,M2,struct());

%% Compute the scaling of the filters
% Improve the scaling of the first and last channel
scal=ones(M2,1);
for ii=1:M2
    scal(ii)=sqrt(afull(ii,1)/afull(ii,2));
end;

%% Construct the real or complex filterbank

if flags.do_real
    % Scale the first and last channels
    scal(1)=scal(1)/sqrt(2);
    scal(M2)=scal(M2)/sqrt(2);
else
    % Replicate the centre frequencies and sampling rates, except the first and
    % last
    if flags.do_nonuniform
        a=[a;flipud(a(2:M2-1,:))];
    end;
    scal=[scal;flipud(scal(2:M2-1))];
    fc  =[fc; -flipud(fc(2:M2-1))];
    if flags.do_symmetric
        fsupp=[fsupp;flipud(fsupp(2:M2-1))];
    end;
    
end;


%% Compute the filters
if flags.do_symmetric
    g=blfilter(flags.wintype,fsupp,fc,'fs',fs,'scal',scal,'inf');    
else
    g=warpedblfilter(flags.wintype,fsupp_erb,fc,fs,@freqtoerb,@erbtofreq, ...
                     'scal',scal,'inf'); 
end;

