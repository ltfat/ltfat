function [g,a,fc]=erbfilters(fs,varargin)
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
%     'N',N           Specify the number of filters, *N*. If this
%                     parameter is specified, it overwrites the
%                     `'spacing'`
%                     parameter.
%
%     'nonuniform'    Construct a non-uniform filterbank. This is the
%                     default.
%    
%     'uniform'       Construct a uniform filterbank.
%
%     'real'          Construct a filterbank that works for real-valued
%                     signals only (the filters cover only the positive
%                     frequencies). This is the default.
%
%     'complex'       Construct a filterbank that covers the entire
%                     frequency range.
%
%     'refsampling'   Choose the downsampling rates to be products of 2
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
%     [g,a,fc]=erbfilters(fs,'uniform','N',100);
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


definput.keyvals.L=[];
definput.keyvals.N=[];
definput.keyvals.bwmul=1;
definput.flags.uniform  = {'nonuniform','uniform'};
definput.flags.real     = {'real','complex'};
definput.flags.sampling = {'regsampling','fractional'};

[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

% Construct the Erb filterbank

N=kv.N;
if isempty(N)
    N=ceil(freqtoerb(fs/2))+1;
end;

fc=erbspace(0,fs/2,N);
% "*3" is just a heuristic, no justification
fsupp=round(audfiltbw(fc)*4*kv.bwmul);

% Improve the scaling of the first and last channel
scal=ones(1,N);
scal(1)=scal(1)/sqrt(2);
scal(N)=scal(N)/sqrt(2);


if flags.do_nonuniform
    % Do the non-uniform case
    % Energy scaling works best
    g=blfilter('hanning',fsupp,fc,'fs',fs,'scal',scal,'2');
    
    % Find suitable channel subsampling rates
    aprecise=round(fs./fsupp/2); % "/2" is just a heuristic, no justification
    aprecise=aprecise(:);
    
    if flags.do_fractional
        Nfilts=round(L./aprecise);
        a=[repmat(L,N,1),Nfilts];                
        
    else
        a=ceil23(aprecise); % Grow "a" to the next composite number
    
        % Determine the minimal transform length
        L=filterbanklength(1,a);
    end;

else
    % Do the uniform case
    % Peak-frequency scaling works best
    g=blfilter('hanning',fsupp,fc,'fs',fs,'scal',scal,'inf');
    a=3;

end;


