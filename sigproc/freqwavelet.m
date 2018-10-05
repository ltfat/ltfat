function [H,info] = freqwavelet(name,L,varargin)
%FREQWAVELET  Wavelet in the freq. domain
%   Usage: H=freqwavelet(name,L)
%          H=freqwavelet(name,L,scale)
%          [H,info]=freqwavelet(...)
%
%   Input parameters:
%         name  : Name of the wavelet
%         L     : System length
%         scale : Wavelet scale
%   Output parameters:
%         H     : Frequency domain wavelet
%         info  : Struct with additional outputs
%
%   `freqwavelet(name,L)` returns peak-normalized "mother" frequency-domain
%   wavelet `name` for system length *L*. The basic scale is selected such
%   that the peak is positioned at the frequency 0.1 relative to the Nyquist
%   rate (fs/2).
%
%   The supported wavelets that can be used in place of `name` (the actual
%   equations can be found at the end of this help):
%
%     'cauchy'    Cauchy wavelet with alpha=100. Custom order=(alpha-1)/2
%                 (with alpha>1) can be set by `{'cauchy',alpha}`. The
%                 higher the order, the narrower and more symmetric the
%                 wavelet. A numericaly stable formula is used in order
%                 to allow support even for very high `alpha` e.g.
%                 milions.
%
%     'morse'     Morse wavelet with alpha=100 and gamma=3. Both parameters
%                 can be set directly by `{'morse',alpha,gamma}`. `alpha`
%                 has the same role as for 'cauchy', `gamma` is the
%                 'skewness' parameter. 'morse' becomes 'cauchy' for
%                 gamma=1.
%
%   `freqwavelet(name,L,scale)` returns a "dilated" version of the wavelet.
%   The positive scalar `scale` controls both the bandwidth and the center
%   frequency. Values greater than 1 make the wavelet wider (and narrower 
%   in the frequency domain), values lower than 1 make the wavelet narrower.
%   The center frequency is moved to`0.1/scale`. The center frequency is 
%   limitted to the range of "positive" frequencies ]0,1] (up to Nyquist rate).
%   If `scale` is a vector, the output is a `L x numel(scale)` matrix with 
%   the individual wavelets as columns.
%
%   The following additional flags and key-value parameters are available:
%
%     'basefc',fc      Normalized center frequency of the mother wavelet
%                      (scale=1). The default is 0.1.
%
%     'bwthr',bwthr    The height at which the bandwidth is computed.
%                      The default value is 10^(-3/10) (~0.5).
%
%     'efsuppthr',thr  The threshold determinig the effective support of
%                      the wavelet. The default value is 10^(-5).
%
%   The format of the output is controlled by the following flags:
%   'full' (default),'econ','asfreqfilter':
%
%     'full'           The output is a `L x numel(scale)` matrix.
%
%     'econ'           The output is a `numel(scale)` cell array with
%                      individual freq. domain wavelets truncated to the
%                      length of the effective support given by parameter 'efsuppthr'.
%
%     'asfreqfilter'   As `'econ'`, but the elements of the cell-array are
%                      filter structs with fields .H and .foff as in 
%                      |blfilter| to be used in |filterbank| and related.
%
%   `[H,info]=freqwavelet(...)` additionally returns a struct with the
%   following fields:
%
%     .fc             Normalized center frequency.
%
%     .foff           Index of the first sample above the effective
%                     support threshold (minus one). 
%                     It can be directly used to convert the 'econ'
%                     output to 'full' by `circshift(postpad(H,L),foff)`.
%
%     .fsupplen       Length of the effective support (with values above efsuppthr.).
%
%     .scale          The scale used.
%
%     .dilation       The actual dilation used in the formula.
%
%     .bw             Relative bandwidth at -3 dB (half of the height).
%
%     .tfr            Time-frequency ratio of a Gaussian with the same
%                     bandwidth as the wavelet.
%
%     .a_natural      Fractional natural subsampling factors in the
%                     format acceptable by |filterbank| and related.
%
%   Additionally, the function accepts flags to normalize the output.
%   Please see the help of |normalize|. By default, no normaliazation is
%   applied.
%
%   Wavelet definitions
%   -------------------
%
%   `C` is a normalization constant.
%
%   Cauchy wavelet
%
%       .. H = C \xi^{\frac{\alpha-1}{2}} exp( -2\pi\xi )
%
%       .. math:: H = C \xi^{\frac{\alpha-1}{2}} exp( -2\pi\xi )
%
%   Morse wavelet
%
%       .. H = C \xi^{\frac{\alpha-1}{2\gamma}} exp( -2\pi\xi^{\gamma} )
%
%       .. math:: H = C \xi^{\frac{\alpha-1}{2\gamma}} exp( -2\pi\xi^{\gamma} )
%
%   See also: normalize, filterbank, blfilter

% AUTHORS: Zdenek Prusa, Nicki Holighaus, Gunther Koliander

complainif_notenoughargs(nargin,2,upper(mfilename));

if ~isscalar(L)
    error('%s: L must be a scalar',upper(mfilename));
end

if ~iscell(name), name = {name}; end

freqwavelettypes = getfield(arg_freqwavelet(),'flags','wavelettype');

if ~ischar(name{1}) || ~any(strcmpi(name{1},freqwavelettypes))
  error('%s: First input argument must the name of a supported window.',...
        upper(mfilename));
end

winArgs = name(2:end);
winName = lower(name{1});

definput.import={'normalize'};
definput.importdefaults={'null'};
definput.keyvals.scale = 1;
definput.keyvals.basefc = 0.1;
definput.keyvals.bwthr = 10^(-3/10);
definput.keyvals.efsuppthr = 10^(-5);
definput.flags.outformat = {'full','econ','asfreqfilter'};
definput.keyvals.scal = 1;
[flags,kv,scale]=ltfatarghelper({'scale'},definput,varargin,'freqwavelet');

if ~isnumeric(L), error('%s: scale must be numeric',upper(mfilename)); end
if any(scale <= 0), error('%s: scale must be positive.', upper(mfilename)); end
if kv.efsuppthr < 0, error('%s: efsuppthr must be nonegative',upper(mfilename)); end
if kv.bwthr < 0, error('%s: bwthr must be nonegative',upper(mfilename)); end
if kv.bwthr < kv.efsuppthr, error('%s: bwthr must be lower than efsuppthr.',upper(mfilename)); end

switch winName
    case 'cauchy'
       definputcauchy.keyvals.alpha=100;
       [~,~,alpha]=ltfatarghelper({'alpha'},definputcauchy,winArgs);
       winName = 'morse';
       winArgs = {'alpha',alpha,'gamma',1};
end

fs = 2;
step = fs/L;
M = numel(scale);
info.fc    = zeros(1,M);
info.foff  = zeros(1,M);
info.fsupp = L*ones(1,M);
info.scale = zeros(1,M);
info.bw    = zeros(1,M);
info.tfr   = zeros(1,M);
info.a_natural = L*ones(M,2);

if flags.do_full
    H = zeros(L,M);
else
    H = cell(1,M);
end

for m = 1:M
    switch winName
       case 'morse'
            definputgenmorse.keyvals.alpha=100;
            definputgenmorse.keyvals.gamma=3;
            [~,~,alpha,gamma]=ltfatarghelper({'alpha','gamma'},definputgenmorse,winArgs);

            if alpha <= 1
                error('%s: Alpha must be larger than 1 (passed alpha=%.2f).',...
                      upper(mfilename),alpha);
            end

            if gamma <= 0
                error('%s: Gamma must be larger than 0 (passed gamma=%.2f).',...
                      upper(mfilename),gamma);
            end

            order = (alpha-1)/(2*gamma);
            peakpos = ( order/(2*pi*gamma) )^(1/(gamma));
            basedil = peakpos/(kv.basefc);

            info.fc(m) = peakpos/(basedil*scale(m));
            info.scale(m) = scale(m);
            info.dilation(m) = basedil*scale(m);

            if info.fc(m) <= 0 || info.fc(m) > 1
                error(['%s: fc out of range [0,1[. Decrease alpha ',...
                       '(passed %.2f) and/or scale (passed %.2f)'],...
                       upper(mfilename),alpha,scale(m));
            end

            freqatheightasc = @(thr) real( (-order/(2*pi*gamma)...
                                     *octave_lambertw( 0, ...
                                     -thr^(gamma/order)/exp(1)))^(1/(gamma)) )...
                                      /basedil/scale(m);
            freqatheightdesc= @(thr) real( (-order/(2*pi*gamma)...
                                     *octave_lambertw(-1, ...
                                     -thr^(gamma/order)/exp(1)))^(1/gamma) )...
                                     /basedil/scale(m);
            fsupp = [0 0 info.fc(m) fs fs];

            if kv.efsuppthr > 0
                fsupp(1) = max( 0,freqatheightasc(kv.efsuppthr));
                fsupp(5) = min(fs,freqatheightdesc(kv.efsuppthr));
            end
            fsupp(2) = max( 0,freqatheightasc(kv.bwthr));
            fsupp(4) = min(fs,freqatheightdesc(kv.bwthr));

            fsuppL = fsuppL_inner(fsupp,fs,L,1:5);

            morsefun = @(y) (y > 0).*exp(-2*pi*y.^gamma + order*log(y) ...
                             + ( order/gamma - order/gamma*log(order/(2*pi*gamma)) ));

            info.foff(m) = fsuppL(1);
            info.fsupp(m) = fsuppL(end) - fsuppL(1) + 1;
            if info.fsupp(m) <= 0, info.fsupp(m) = 0; end

            info.bw(m)  = (fsupp(4) - fsupp(2));
            bwinsamples = info.bw(m)/step;
            info.a_natural(m,2) = ceil(bwinsamples);
            %info.tfr(m) = 1/(pi/(4*log(1/kv.bwthr))*bwinsamples^2/L);
            info.tfr(m) = (alpha - 1)/(pi*info.fc(m)^2*L);
                         
            if flags.do_full
                y = ((0:L-1)').*basedil*step*scale(m);
                H(:,m) = kv.scal*normalize(morsefun(y), flags.norm);
            elseif flags.do_econ
                y = ((fsuppL(1):fsuppL(end)-1)').*basedil*step*scale(m);
                H{m} = kv.scal*normalize(morsefun(y), flags.norm);
            elseif flags.do_asfreqfilter
                y = @(L) ((fsuppL_inner(fsupp,fs,L,1):fsuppL_inner(fsupp,fs,L,5)-1)').*basedil*scale(m)*fs/L;
                H{m} = struct('H',@(L) kv.scal*normalize(morsefun(y(L)),flags.norm),'foff',@(L)fsuppL_inner(fsupp,fs,L,1),'realonly',0);
            end


        otherwise
            error('%s: SENTINEL. Unknown window.',upper(mfilename));
    end
end

if M==1 && iscell(H)
    H = H{1};
end

function fsuppL = fsuppL_inner(fsupp,fs,L,idx)
fsuppL_all = [ ceil(fsupp(1:2)/fs*L), round(fsupp(3)/fs*L), floor(fsupp(4:5)/fs*L) ];
fsuppL = fsuppL_all(idx);


function w = octave_lambertw(b,z)
% Copyright (C) 1998 by Nicol N. Schraudolph <schraudo@inf.ethz.ch>
%
% @deftypefn {Function File} {@var{x} = } lambertw (@var{z})
% @deftypefnx {Function File} {@var{x} = } lambertw (@var{n}, @var{z})
% Compute the Lambert W function of @var{z}.
%
% This function satisfies W(z).*exp(W(z)) = z, and can thus be used to express
% solutions of transcendental equations involving exponentials or logarithms.
%
% @var{n} must be integer, and specifies the branch of W to be computed;
% W(z) is a shorthand for W(0,z), the principal branch.  Branches
% 0 and -1 are the only ones that can take on non-complex values.
%
% If either @var{n} or @var{z} are non-scalar, the function is mapped to each
% element; both may be non-scalar provided their dimensions agree.
%
% This implementation should return values within 2.5*eps of its
% counterpart in Maple V, release 3 or later.  Please report any
% discrepancies to the author, Nici Schraudolph <schraudo@@inf.ethz.ch>.

if (nargin == 1)
    z = b;
    b = 0;
else
    %% some error checking
    if (nargin ~= 2)
        print_usage;
    else
        if (any(round(real(b)) ~= b))
            usage('branch number for lambertw must be integer')
        end
    end
end

%% series expansion about -1/e
%
% p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
% w = (11/72)*p;
% w = (w - 1/3).*p;
% w = (w + 1).*p - 1
%
% first-order version suffices:
%
w = (1 - 2*abs(b)).*sqrt(2*exp(1)*z + 2) - 1;

%% asymptotic expansion at 0 and Inf
%
v = log(z + double(~(z | b))) + 2*pi*1i*b;
v = v - log(v + double(v==0));

%% choose strategy for initial guess
%
c = abs(z + 1/exp(1));
c = (c > 1.45 - 1.1*abs(b));
c = c | (b.*imag(z) > 0) | (~imag(z) & (b == 1));
w = (1 - c).*w + c.*v;

%% Halley iteration
%
for n = 1:10
    p = exp(w);
    t = w.*p - z;
    f = (w ~= -1);
    t = f.*t./(p.*(w + f) - 0.5*(w + 2.0).*t./(w + f));
    w = w - t;
    if (abs(real(t)) < (2.48*eps)*(1.0 + abs(real(w))) ...
        && abs(imag(t)) < (2.48*eps)*(1.0 + abs(imag(w))))
        return
    end
end

%error('PRECISION:iteration limit reached, result of lambertw may be inaccurate');



