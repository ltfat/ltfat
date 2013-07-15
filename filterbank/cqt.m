function [c,Ls,g,shift,M] = cqt(f,fmin,fmax,bins,fs,varargin)
%CQT  Constant-Q nonstationary Gabor filterbank
%   Usage: [c,Ls,g,shift,M] = cqt(f,fmin,fmax,bins,fs,M)
%          [c,Ls,g,shift,M] = cqt(f,fmin,fmax,bins,fs)
%          [c,Ls,g,shift] = cqt(...)
%          [c,Ls] = cqt(...)
%          c = cqt(...)
%
%   Input parameters: 
%         f         : The signal to be analyzed (For multichannel
%                     signals, input should be a matrix which each
%                     column storing a channel of the signal).
%         fmin      : Minimum frequency (in Hz)
%         fmax      : Maximum frequency (in Hz)
%         bins      : Vector consisting of the number of bins per octave
%         fs        : Sampling rate (in Hz)
%         M         : Number of time channels (optional)
%                     If M is constant, the output is converted to a
%                     matrix
%   Output parameters:
%         c         : Transform coefficients (matrix or cell array)
%         Ls        : Original signal length (in samples)
%         g         : Cell array of Fourier transforms of the analysis 
%                     windows
%         shift     : Vector of frequency shifts
%         M         : Number of time channels
%
%   This function computes a constant-Q transform via nonstationary Gabor
%   filterbanks. Given the signal *f*, the constant-Q parameters *fmin*,
%   *fmax* and *bins*, as well as the sampling rate *fs* of *f*, the
%   corresponding constant-Q coefficients *c* are given as output. For
%   reconstruction, the length of *f* and the filterbank parameters can
%   be returned also.
% 
%   The transform produces phase-locked coefficients in the
%   sense that each filter is considered to be centered at
%   0 and the signal itself is modulated accordingly.
%
%   Optional input arguments arguments can be supplied like this::
%       
%       cqt(f,fmin,fmax,bins,fs,'min_win',min_win)
%
%   The arguments must be character strings followed by an
%   argument:
%
%     'min_win',min_win        Minimum admissible window length 
%                              (in samples) 
%
%     'Qvar',Qvar              Bandwidth variation factor
%
%     'M_fac',M_fac            Number of time channels are rounded to 
%                              multiples of this
%
%     'winfun',winfun          Filter prototype (see |firwin| for available 
%                              filters)
%     'fractional'             Allow fractional shifts and bandwidths
%
%
%   Example:
%   --------
%
%   The following example shows analysis and synthesis with |cqt| and |icqt|:::
%
%     f = gspi;
%     [c,Ls,g,shift,M] = cqt(f,50,20000,12,44100);
%     fr = icqt(c,g,shift,Ls);
%     rel_err = norm(f-fr)/norm(f)
%     plotfilterbank(c,shift,'dynrange',60);
%
%   See also:  icqt, firwin
% 
%   References:  dogrhove11 dogrhove12

% Authors: Nicki Holighaus, Gino Velasco
% Date: 10.04.13

%% Check input arguments
if nargin < 5
    error('Not enough input arguments');
end

[f,Ls,W]=comp_sigreshape_pre(f,upper(mfilename),0);

% Set defaults

definput.keyvals.usrM = [];
definput.keyvals.Qvar = 1;
definput.keyvals.M_fac = 1;
definput.keyvals.min_win = 4;
definput.keyvals.winfun = 'hann';
definput.flags.fractype = {'nofractional','fractional'};

% Check input arguments

[flags,keyvals,usrM]=ltfatarghelper({'usrM'},definput,varargin);

%% Create the CQ-NSGT dictionary

nf = fs/2;

if fmax > nf
    fmax = nf;
end

b = ceil(log2(fmax/fmin))+1;

if length(bins) == 1;
    bins = bins*ones(b,1);
elseif length(bins) < b
    if size(bins,1) == 1
        bins=bins.';
    end
    bins( bins<=0 ) = 1;
    bins = [bins ; bins(end)*ones(b-length(bins),1)];
end

fbas = zeros(sum(bins),1);

ll = 0;
for kk = 1:length(bins);
    fbas(ll+(1:bins(kk))) = ...
        fmin*2.^(((kk-1)*bins(kk):(kk*bins(kk)-1)).'/bins(kk));
    ll = ll+bins(kk);
end

temp = find(fbas>=fmax,1);
if fbas(temp) >= nf
    fbas = fbas(1:temp-1);
else
    fbas = fbas(1:temp);
end

Lfbas = length(fbas);

fbas = [0;fbas];
fbas(Lfbas+2) = nf;
fbas(Lfbas+3:2*(Lfbas+1)) = fs-fbas(Lfbas+1:-1:2);

fbas = fbas*(Ls/fs);

% Set bandwidths

bw = zeros(2*Lfbas+2,1);

bw(1) = 2*fmin*(Ls/fs);
bw(2) = (fbas(2))*(2^(1/bins(1))-2^(-1/bins(1)));

for k = [3:Lfbas , Lfbas+2]
    bw(k) = (fbas(k+1)-fbas(k-1));
end

bw(Lfbas+1) = (fbas(Lfbas+1))*(2^(1/bins(end))-2^(-1/bins(end)));
bw(Lfbas+3:2*Lfbas+2) = bw(Lfbas+1:-1:2);

posit = zeros(size(fbas));
posit(1:Lfbas+2) = floor(fbas(1:Lfbas+2));
posit(Lfbas+3:end) = ceil(fbas(Lfbas+3:end));

bw = keyvals.Qvar*bw;

if flags.do_fractional
    warning(['Fractional sampling might lead to a warning when ', ...
        'computing the dual system']);
    fprintf('');
    corr_shift = fbas-posit;
    M = ceil(bw+1);
else
    bw = round(bw);
    M = bw;
end

for ii = 1:2*(Lfbas+1)
    if bw(ii) < keyvals.min_win;
        bw(ii) = keyvals.min_win;
        M(ii) = bw(ii);
    end
end

if flags.do_fractional
    g = arrayfun(@(x,y,z) ...
        firwin(keyvals.winfun,([0:ceil(z/2),-floor(z/2):-1]'-x)/y)/sqrt(y),corr_shift,...
        bw,M,'UniformOutput',0);
else
    g = arrayfun(@(x) firwin(keyvals.winfun,x)/sqrt(x),...
        bw,'UniformOutput',0);
end

M = keyvals.M_fac*ceil(M/keyvals.M_fac);

% Setup Tukey window for 0- and Nyquist-frequency
for kk = [1,Lfbas+2]
    if M(kk) > M(kk+1);
        g{kk} = ones(M(kk),1);
        g{kk}((floor(M(kk)/2)-floor(M(kk+1)/2)+1):(floor(M(kk)/2)+...
            ceil(M(kk+1)/2))) = firwin('hann',M(kk+1));
        g{kk} = g{kk}/sqrt(M(kk));
    end
end

N = length(posit);  % The number of frequency channels

if ~isempty(usrM)
    if numel(usrM) == 1
        M = usrM*ones(N,1);
    else
        M = usrM;
    end    
end

%% The CQ-NSG transform

% some preparation
f = fft(f);

c=cell(N,1); % Initialisation of the result

% Obtain input type
ftype = assert_classname(f);
% The actual transform

for ii = 1:N
    Lg = length(g{ii});
    
    idx = [ceil(Lg/2)+1:Lg,1:ceil(Lg/2)];
    win_range = mod(posit(ii)+(-floor(Lg/2):ceil(Lg/2)-1),Ls)+1;
    
    if M(ii) < Lg % if the number of frequency channels is too small,
        % aliasing is introduced
        col = ceil(Lg/M(ii));
        temp = zeros(col*M(ii),W,ftype);
        
        temp([end-floor(Lg/2)+1:end,1:ceil(Lg/2)],:) = ...
            bsxfun(@times,f(win_range,:),g{ii}(idx));
        temp = reshape(temp,M(ii),col,W);
        
        c{ii} = squeeze(ifft(sum(temp,2)));
        
        % Using c = cellfun(@(x) squeeze(ifft(x)),c,'UniformOutput',0);
        % outside the loop instead does not provide speedup; instead it is
        % slower in most cases.
    else
        temp = zeros(M(ii),W,ftype);
        temp([end-floor(Lg/2)+1:end,1:ceil(Lg/2)],:) = ...
            bsxfun(@times,f(win_range,:),g{ii}(idx));
        
        c{ii} = ifft(temp);
    end
end

if max(M) == min(M)
    c = cell2mat(c);
    c = reshape(c,M(1),N,W);
end

if nargout > 3
    shift = [mod(-posit(end),Ls); diff(posit)];
end
