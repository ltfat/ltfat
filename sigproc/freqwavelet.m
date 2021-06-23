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
%                 10^6.
%
%     'morse'     Morse wavelet with alpha=100 and gamma=3. Both parameters
%                 can be set directly by `{'morse',alpha,gamma}`. `alpha`
%                 has the same role as for 'cauchy', `gamma` is the
%                 'skewness' parameter. 'morse' becomes 'cauchy' for
%                 gamma=1.
%
%     'morlet'    Morlet wavelet with sigma=4. The parameter `sigma` 
%                 is the center frequency to standard deviation ratio of
%                 the main generating Gaussian distribution. Note that the 
%                 true peak and standard deviation of the Morlet wavelet 
%                 differ, in particular for low values of `sigma` (<5). 
%                 This is a consequence of the correction factor. The 
%                 parameter can be set directly by `{'morlet',sigma}`.
% 
%    'fbsp'       Frequency B-spline wavelet of order m=3, center frequency
%                 to bandwidth ratio fb = 2. The parameters can be set
%                 by `{'fbsp',m,fb}`. Note that `m` must be integer and
%                 greater than or equal to 1, and `fb` must be greater than 
%                 or equal to 2. 
%
%    'analyticsp' Positive frequency part of cosine-modulated B-spline 
%                 wavelet of  order order=3, with center frequency to main 
%                 lobe width ratio fb = 1. The parameters can be set by 
%                 `{'analyticsp',order,fb}`. Note that `order` and `fb` 
%                 must be integer and greater than or equal to 1.
%
%    'cplxsp'     Complex-modulated B-spline of order order=3, with center
%                 frequency to main lobe width ratio fb = 1. The parameters 
%                 can be set by `{'cplxsp',order,fb}`. Note that `order` 
%                 and `fb` must be integer and greater than or equal to 1.
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
%     'waveletParams'   a vector containing the respective wavelet parameters
%                       [alpha, beta, gamma] for cauchy and morse wavelets
%                       [sigma] for morlet wavelets
%                       [order, fb] for splines 
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
%     'asfreqfilter'   As 'econ', but the elements of the cell-array are
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
%     .fsupp       Length of the effective support (with values above efsuppthr.).basefc
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
%       H = C \xi^{\frac{\alpha-1}{2}} exp( -2\pi\xi )
%
%       .. math:: H = C \xi^{\frac{\alpha-1}{2}} exp( -2\pi\xi )
%
%   Morse wavelet
%
%       H = C \xi^{\frac{\alpha-1}{2\gamma}} exp( -2\pi\xi^{\gamma} )
%
%       .. math:: H = C \xi^{\frac{\alpha-1}{2\gamma}} exp( -2\pi\xi^{\gamma} )
%
%   See also: normalize, filterbank, blfilter



complainif_notenoughargs(nargin,2,upper(mfilename));

%set default input parameters
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

definput.import={'normalize', 'freqwavelet'};
definput.importdefaults={'null'};
definput.keyvals.scale = 1;
%definput.keyvals.scal = 1;
definput.keyvals.basefc = 0.1;
definput.keyvals.bwthr = 10^(-3/10);
definput.keyvals.efsuppthr = 10^(-5);
definput.flags.outformat = {'full','econ','asfreqfilter'};
definput.keyvals.fs = 2;
definput.keyvals.alphaStep = definput.keyvals.fs/L;

switch lower(winName)
  case 'cauchy'
    definput.keyvals.waveletParams = [100,0,3];
  case 'morse'
    definput.keyvals.waveletParams = [100,0,3];  
  case 'morlet'
    definput.keyvals.waveletParams = [4];  
  case 'fbsp'
    definput.keyvals.waveletParams = [4, 2];  
  case 'analyticsp'
    definput.keyvals.waveletParams = [4, 2];
  case 'cplxsp'
    definput.keyvals.waveletParams = [4, 2];  
end 

[flags,kv,scale]=ltfatarghelper({'scale'},definput,varargin,'freqwavelet');

%check L
if ~isscalar(L)
    error('%s: L must be a scalar',upper(mfilename));
end

if ~isnumeric(L), error('%s: scale must be numeric',upper(mfilename)); end


if any(scale <= 0), error('%s: scale must be positive.', upper(mfilename)); end
if kv.efsuppthr < 0, error('%s: efsuppthr must be nonegative',upper(mfilename)); end
if kv.bwthr < 0, error('%s: bwthr must be nonegative',upper(mfilename)); end
if kv.bwthr < kv.efsuppthr, error('%s: bwthr must be lower than efsuppthr.',upper(mfilename)); end

%initialize the output info

M = numel(scale);
info.fc    = zeros(1,M);
info.foff  = zeros(1,M);
info.fsupp = L*ones(1,M);
info.scale = zeros(1,M);
info.dilation = 0;
info.bw    = zeros(1,M);
info.tfr   = zeros(1,M);
info.a_natural = L*ones(M,2);

if flags.do_full
    H = zeros(L,M);
else
    H = cell(1,M);
end

for m = 1:M
% Generalized Morse wavelets
        if strcmp('morse', lower(winName)) || strcmp('cauchy', lower(winName))
            
            definputgenmorse.keyvals.alpha = kv.waveletParams(1);
            definputgenmorse.keyvals.beta = kv.waveletParams(2);
            definputgenmorse.keyvals.gamma = kv.waveletParams(3);
            
            [~,~,alpha,beta,gamma]=ltfatarghelper({'alpha','beta','gamma'},definputgenmorse,winArgs);

            if alpha <= 1
                error('%s: Alpha must be larger than 1 (passed alpha=%.2f).',...
                      upper(mfilename),alpha);
            end

            if gamma <= 0
                error('%s: Gamma must be larger than 0 (passed gamma=%.2f).',...
                      upper(mfilename),gamma);
            end

            %derive the wavelet function
            order = (alpha-1)/(2*gamma);
            peakpos = ( order/(2*pi*gamma) )^(1/(gamma));
            basedil = peakpos/(kv.basefc);

            for m = 1:M
              freqatheightasc = @(thr) real( (-order/(2*pi*gamma)...
                                     *octave_lambertw( 0, ...
                                     -thr^(gamma/order)/exp(1)))^(1/(gamma)) )...
                                      /basedil/scale(m);
              freqatheightdesc= @(thr) real( (-order/(2*pi*gamma)...
                                     *octave_lambertw(-1, ...
                                     -thr^(gamma/order)/exp(1)))^(1/gamma) )...
                                     /basedil/scale(m);


              fun = @(y) (y > 0).*exp(-2*pi*y.^gamma + (order - 1i*beta)*log(y) ...
                             + ( order/gamma - order/gamma*log(order/(2*pi*gamma)) ));
              
              info.fc(m) = peakpos/(basedil*scale(m));
              info.scale(m) = scale(m);
              info.dilation(m) = basedil*scale(m);
              
              fsupp = [0 0 info.fc(m) kv.fs kv.fs];

              if kv.efsuppthr > 0
                fsupp(1) = max( 0,freqatheightasc(kv.efsuppthr));
                fsupp(5) = min(kv.fs,freqatheightdesc(kv.efsuppthr));
              end
              fsupp(2) = max( 0,freqatheightasc(kv.bwthr));
              fsupp(4) = min(kv.fs,freqatheightdesc(kv.bwthr));

              fsuppL = fsuppL_inner(fsupp,kv.fs,L,1:5);
              
              info.foff(m) = fsuppL(1);
              info.fsupp(m) = fsuppL(end) - fsuppL(1) + 1;
              if info.fsupp(m) <= 0, info.fsupp(m) = 0; end

              info.bw(m)  = (fsupp(4) - fsupp(2));
              bwinsamples = info.bw(m)/kv.alphaStep;
              info.a_natural(m,2) = ceil(bwinsamples);
 
%              CauchyAlpha = wpghi_findalpha({'morse',order},0.2);
              CauchyAlpha = alpha;
              info.tfr(m) = (CauchyAlpha - 1)/(pi*info.fc(m)^2*L);
              info.CauchyAlpha(m) = CauchyAlpha;
              
              if flags.do_full
                y = ((0:L-1)').*basedil*kv.alphaStep*scale(m);
                H(:,m) = scale(m)*normalize(fun(y), flags.norm);
              elseif flags.do_econ
                y = ((fsuppL(1):fsuppL(end)-1)').*basedil*kv.alphaStep*scale(m);
                H{m} = scale(m)*normalize(fun(y), flags.norm);
              elseif flags.do_asfreqfilter
                y = @(L) ((fsuppL_inner(fsupp,kv.fs,L,1):fsuppL_inner(fsupp,kv.fs,L,5)-1)').*basedil*scale(m)*kv.fs/L;
                H{m} = struct('H',@(L) scale(m)*normalize(fun(y(L)),flags.norm),'foff',@(L)fsuppL_inner(fsupp,kv.fs,L,1),'realonly',0);
              end
            end 
        
% Morlet wavelets
      elseif strcmp('morlet', lower(winName))
                  
            definputmorlet.keyvals.sigma=kv.waveletParams(1);
            [~,~,sigma]=ltfatarghelper({'sigma'},definputmorlet,winArgs);

            if sigma <= 1
                error('%s: Sigma must be larger than 1 (passed sigma=%.2f).',...
                      upper(mfilename),sigma);
            end

            %derive the wavelet parameters
            % Fixed point iteration to find the maximum of the Morlet
            % wavelet
            peakpos = sigma;
            peakpos_tmp = 0;
            while abs(peakpos-peakpos_tmp) > 1e-6
                peakpos_tmp = peakpos;
                peakpos = sigma./(1-exp(-sigma*peakpos));
            end
            basedil = peakpos/(kv.basefc);
            
            fun = @(y) ( exp(-0.5*(sigma-y).^2) - exp(-0.5*( sigma.^2+y.^2 )) )...
                               ./ ( exp(-0.5*(sigma-peakpos).^2) - exp(-0.5*( sigma.^2+peakpos.^2 )) );            

            for m = 1:M
              freqatheightdesc = @(thr) determine_freqatheight(fun,peakpos,thr,1)/basedil/scale(m);
              freqatheightasc= @(thr) determine_freqatheight(fun,peakpos,thr,0)/basedil/scale(m);
              
              info.fc(m) = peakpos/(basedil*scale(m));
              info.scale(m) = scale(m);
              info.dilation(m) = basedil*scale(m);
              
              fsupp = [0 0 info.fc(m) kv.fs kv.fs];

              if kv.efsuppthr > 0
                fsupp(1) = max( 0,freqatheightasc(kv.efsuppthr));
                fsupp(5) = min(kv.fs,freqatheightdesc(kv.efsuppthr));
              end
              fsupp(2) = max( 0,freqatheightasc(kv.bwthr));
              fsupp(4) = min(kv.fs,freqatheightdesc(kv.bwthr));

              fsuppL = fsuppL_inner(fsupp,kv.fs,L,1:5);
              
              info.foff(m) = fsuppL(1);
              info.fsupp(m) = fsuppL(end) - fsuppL(1) + 1;
              if info.fsupp(m) <= 0, info.fsupp(m) = 0; end

              info.bw(m)  = (fsupp(4) - fsupp(2));
              bwinsamples = info.bw(m)/kv.alphaStep;
              info.a_natural(m,2) = ceil(bwinsamples);
 
              CauchyAlpha = wpghi_findalpha({'morlet',sigma},0.2);
              info.tfr(m) = (CauchyAlpha - 1)/(pi*info.fc(m)^2*L);
              info.CauchyAlpha(m) = CauchyAlpha;
              
              if flags.do_full
                y = ((0:L-1)').*basedil*kv.alphaStep*scale(m);
                H(:,m) = scale(m)*normalize(fun(y), flags.norm);
              elseif flags.do_econ
                y = ((fsuppL(1):fsuppL(end)-1)').*basedil*kv.alphaStep*scale(m);
                H{m} = scale(m)*normalize(fun(y), flags.norm);
              elseif flags.do_asfreqfilter
                y = @(L) ((fsuppL_inner(fsupp,kv.fs,L,1):fsuppL_inner(fsupp,kv.fs,L,5)-1)').*basedil*scale(m)*kv.fs/L;
                H{m} = struct('H',@(L) scale(m)*normalize(fun(y(L)),flags.norm),'foff',@(L)fsuppL_inner(fsupp,kv.fs,L,1),'realonly',0);
              end
            end 

       elseif strcmp('sp', lower(winName(end-1:end)))
%if it is a spline...
           definputfbsp.keyvals.order=kv.waveletParams(1);
           definputfbsp.keyvals.fb=kv.waveletParams(2);
           [~,~,order,fb]=ltfatarghelper({'order','fb'},definputfbsp,winArgs);

            if order < 1 || order > 5 || round(order) ~= order
                error('%s: order must be integer and between 1 and 5 (passed order=%.2f).',...
                      upper(mfilename),order);
            end
            
            if fb < 1 || round(fb) ~= fb
              if strcmp('fbsp', winName) && fb < 2
                      error('%s: fb must be an integer and at least 2 (passed fb=%.2f).',...
                      upper(mfilename),fb);
              else
                error('%s: fb must be an integer and at least 1 (passed fb=%.2f).',...
                      upper(mfilename),fb);
              end
            end
            
            peakpos = 1;
            basedil = peakpos/(kv.basefc);
            
        switch lower(winName)
          
          case 'fbsp'
                   
            switch order
                case 1
                  prefun = @(x) ( x >= 0 ).*( x < 1 ) .* 1; 
                  
                case 2
                  prefun = @(x) ( x >= 0 ).*( x < 1 ) .* ...
                                   (x) ...
                                + ( x >= 1 ).*( x < 2 ) .* ...
                                   (2-x);
                                   
                case 3
                  prefun = @(x)  ( x >= 0 ).*( x < 1 ) .* ...
                                   (.5*x.^2) ...
                                + ( x >= 1 ).*( x < 2 ) .* ...
                                   (-x.^2 + 3.*x -1.5) ...
                                + ( x >= 2 ).*( x < 3 ) .* ...
                                   (.5*x.^2 - 3.*x + 4.5);
                case 4
                  prefun = @(x)  ( x >= 0 ).*( x < 1 ) .* ...
                                   (x.^3./6) ...
                                + ( x >= 1 ).*( x < 2 ) .* ...
                                   (-x.^3./2 + 2.*x.^2 - 2.*x + 2/3) ...
                                + ( x >= 2 ).*( x < 3 ) .* ...
                                   (x.^3./2 - 4.*x.^2 + 10.*x - 22/3) ...
                                + ( x >= 3 ).*( x < 4 ) .* ...
                                   (-x.^3./6 + 2.*x.^2 - 8.*x + 32/3); 
                case 5
                  prefun = @(x) ( x >= 0 ).*( x < 1 ) .* ...
                                   (x.^4./24) ...
                                + ( x >= 1 ).*( x < 2 ) .* ...
                                   (-x.^4./6 + 5.*x.^3./6 - 5.*x.^2./4 + 5.*x./6 - 5/24) ...
                                + ( x >= 2 ).*( x < 3 ) .* ...
                                   (x.^4./4 - 5.*x.^3./2 + 35.*x.^2./4 - 25.*x./2 + 155/24) ... 
                                + ( x >= 3 ).*( x < 4 ) .* ...
                                   (-x.^4./6 + 5.*x.^3./2 - 55.*x.^2./4 + 65.*x./2 - 655/24) ... 
                                + ( x >= 4 ).*( x < 5 ) .* ...
                                   (x.^4./24 -5.*x.^3./6 + 25.*x.^2./4 - 125.*x./6 + 625/24);                                   
                          
               end
            
            fun = @(y) prefun((y-1).*fb.*order./2+order/2)./prefun(order/2);
            
            
            for m = 1:M
              freqatheightdesc = @(thr) determine_freqatheight(fun,peakpos,thr,1)/basedil/scale(m);
              freqatheightasc= @(thr) determine_freqatheight(fun,peakpos,thr,0)/basedil/scale(m);
    
              info.fc(m) = peakpos/(basedil*scale(m));
              info.scale(m) = scale(m);
              info.dilation(m) = basedil*scale(m);
    
              fsupp = [0 0 info.fc(m) kv.fs kv.fs];

              if kv.efsuppthr > 0
                fsupp(1) = max( 0,freqatheightasc(kv.efsuppthr));
                fsupp(5) = min(kv.fs,freqatheightdesc(kv.efsuppthr));
              end
              fsupp(2) = max( 0,freqatheightasc(kv.bwthr));
              fsupp(4) = min(kv.fs,freqatheightdesc(kv.bwthr));

              fsuppL = fsuppL_inner(fsupp,kv.fs,L,1:5);
              
              info.foff(m) = fsuppL(1);
              info.fsupp(m) = fsuppL(end) - fsuppL(1) + 1;
              if info.fsupp(m) <= 0, info.fsupp(m) = 0; end

              info.bw(m)  = (fsupp(4) - fsupp(2));
              bwinsamples = info.bw(m)/kv.alphaStep;
              info.a_natural(m,2) = ceil(bwinsamples);
 
              CauchyAlpha = wpghi_findalpha({'fbsp',order,fb},0.2);
              info.tfr(m) = (CauchyAlpha - 1)/(pi*info.fc(m)^2*L);
              info.CauchyAlpha(m) = CauchyAlpha;
              
              if flags.do_full
                y = ((0:L-1)').*basedil*kv.alphaStep*scale(m);
                H(:,m) = scale(m)*normalize(fun(y), flags.norm);
              elseif flags.do_econ
                y = ((fsuppL(1):fsuppL(end)-1)').*basedil*kv.alphaStep*scale(m);
                H{m} = scale(m)*normalize(fun(y), flags.norm);
              elseif flags.do_asfreqfilter
                y = @(L) ((fsuppL_inner(fsupp,kv.fs,L,1):fsuppL_inner(fsupp,kv.fs,L,5)-1)').*basedil*scale(m)*kv.fs/L;
                H{m} = struct('H',@(L) scale(m)*normalize(fun(y(L)),flags.norm),'foff',@(L)fsuppL_inner(fsupp,kv.fs,L,1),'realonly',0);
              end
            end
           
         case 'analyticsp' % Positive frequency part of cosine-modulated B-spline
            
            fun = @(y) (y>0).* (sinc( fb.*(y - 1) ).^order + sinc( fb.*( y + 1) ).^order);
            
            for m = 1:M    
              heightfun = @(y) min(1,(y>20).* ( 1./abs(fb.*(pi.*y - pi)+eps).^order + 1./abs(fb.*( pi.*y + pi )).^order ));
              freqatheightdesc = @(thr) determine_freqatheight(heightfun,peakpos,thr,1)/basedil/scale(m);
              freqatheightasc= @(thr) determine_freqatheight(heightfun,peakpos,thr,0)/basedil/scale(m);
              
              info.fc(m) = peakpos/(basedil*scale(m));
              info.scale(m) = scale(m);
              info.dilation(m) = basedil*scale(m);
              
              fsupp = [0 0 info.fc(m) kv.fs kv.fs];

              if kv.efsuppthr > 0
                fsupp(1) = max( 0,freqatheightasc(kv.efsuppthr));
                fsupp(5) = min(kv.fs,freqatheightdesc(kv.efsuppthr));
              end
              fsupp(2) = max( 0,freqatheightasc(kv.bwthr));
              fsupp(4) = min(kv.fs,freqatheightdesc(kv.bwthr));

              fsuppL = fsuppL_inner(fsupp,kv.fs,L,1:5);
     
              info.foff(m) = fsuppL(1);
              info.fsupp(m) = fsuppL(end) - fsuppL(1) + 1;
              if info.fsupp(m) <= 0, info.fsupp(m) = 0; end

              info.bw(m)  = (fsupp(4) - fsupp(2));
              bwinsamples = info.bw(m)/kv.alphaStep;
              info.a_natural(m,2) = ceil(bwinsamples);
 
              CauchyAlpha = wpghi_findalpha({'analyticsp',order,fb},0.2);
              info.tfr(m) = (CauchyAlpha - 1)/(pi*info.fc(m)^2*L); 
              info.CauchyAlpha(m) = CauchyAlpha;    
              
              if flags.do_full
                y = ((0:L-1)').*basedil*kv.alphaStep*scale(m);
                H(:,m) = scale(m)*normalize(fun(y), flags.norm);
              elseif flags.do_econ
                y = ((fsuppL(1):fsuppL(end)-1)').*basedil*kv.alphaStep*scale(m);
                H{m} = scale(m)*normalize(fun(y), flags.norm);
              elseif flags.do_asfreqfilter
                y = @(L) ((fsuppL_inner(fsupp,kv.fs,L,1):fsuppL_inner(fsupp,kv.fs,L,5)-1)').*basedil*scale(m)*kv.fs/L;
                H{m} = struct('H',@(L) scale(m)*normalize(fun(y(L)),flags.norm),'foff',@(L)fsuppL_inner(fsupp,kv.fs,L,1),'realonly',0);
              end
            end
            
         case 'cplxsp' % Complex-modulated B-Spline
            
            fun = @(y) sinc( fb.*(y - 1) ).^order;
            
            for m = 1:M  
              heightfun = @(y) min(1,1./abs(fb.*(pi*y - pi)+eps).^order);
              freqatheightdesc = @(thr) determine_freqatheight(heightfun,peakpos,thr,1)/(basedil+eps)/scale(m);
              freqatheightasc= @(thr) determine_freqatheight(heightfun,peakpos,thr,0)/(basedil+eps)/scale(m);
              
              info.fc(m) = peakpos/(basedil*scale(m));
              info.scale(m) = scale(m);
              info.dilation(m) = basedil*scale(m);
              
              fsupp = [0 0 info.fc(m) kv.fs kv.fs];

              if kv.efsuppthr > 0
                fsupp(1) = max( 0,freqatheightasc(kv.efsuppthr));
                fsupp(5) = min(kv.fs,freqatheightdesc(kv.efsuppthr));
              end
              fsupp(2) = max( 0,freqatheightasc(kv.bwthr));
              fsupp(4) = min(kv.fs,freqatheightdesc(kv.bwthr));

              fsuppL = fsuppL_inner(fsupp,kv.fs,L,1:5);
              
              info.foff(m) = fsuppL(1);
              info.fsupp(m) = fsuppL(end) - fsuppL(1) + 1;
              if info.fsupp(m) <= 0, info.fsupp(m) = 0; end

              info.bw(m)  = (fsupp(4) - fsupp(2));
              bwinsamples = info.bw(m)/kv.alphaStep;
              info.a_natural(m,2) = ceil(bwinsamples);
 
              CauchyAlpha = wpghi_findalpha({'cplxsp',order,fb},0.2);
              info.tfr(m) = (CauchyAlpha - 1)/(pi*info.fc(m)^2*L);
              info.CauchyAlpha(m) = CauchyAlpha;
              
              if flags.do_full
                y = ((0:L-1)').*basedil*kv.alphaStep*scale(m);
                H(:,m) = scale(m)*normalize(fun(y), flags.norm);
              elseif flags.do_econ
                y = ((fsuppL(1):fsuppL(end)-1)').*basedil*kv.alphaStep*scale(m);
                H{m} = scale(m)*normalize(fun(y), flags.norm);
              elseif flags.do_asfreqfilter
                y = @(L) ((fsuppL_inner(fsupp,kv.fs,L,1):fsuppL_inner(fsupp,kv.fs,L,5)-1)').*basedil*scale(m)*kv.fs/L;
                H{m} = struct('H',@(L) scale(m)*normalize(fun(y(L)),flags.norm),'foff',@(L)fsuppL_inner(fsupp,kv.fs,L,1),'realonly',0);
              end
            end 
                         
        end 

        end
        

end 


if M==1 && iscell(H)
    H = H{1};
end
end 


function fsuppL = fsuppL_inner(fsupp,fs,L,idx)
fsuppL_all = [ ceil(fsupp(1:2)/fs*L), round(fsupp(3)/fs*L), floor(fsupp(4:5)/fs*L) ];
fsuppL = fsuppL_all(idx);
end 


% function alpha = determine_alpha_from_bandwidth(bwatthr,bwthr,basefc,steps)
% % This function computes alpha from a bandwidth `bwatthr` at a reference height `bwthr`, together with
% % a given base center frequency `basefc`.
%    
% cauchybwatthr = @(alph) basefc * ...
%                           ( octave_lambertw(0, -bwthr^(2/(alph-1))/exp(1))...
%                            -octave_lambertw(-1,-bwthr^(2/(alph-1))/exp(1)) );
% 
% alpha_current = 10;
% cauchybw_current = cauchybwatthr(alpha_current);
% 
% % Find initial guess
% if cauchybw_current > bwatthr
%     while cauchybw_current > bwatthr
%         alpha_left = alpha_current;
%         alpha_current = 10*alpha_current;
%         cauchybw_current = cauchybwatthr(alpha_current);
%     end
% elseif cauchybw_current < bwatthr
%     while cauchybw_current < bwatthr
%         alpha_current = 0.1*alpha_current;
%         alpha_left = alpha_current;
%         cauchybw_current = cauchybwatthr(alpha_current);
%     end
% else 
%     alpha = alpha_current;
%     return
% end
% 
% for kk = 1:steps
%    exponent = 2^(-kk); 
%    alpha_current = alpha_left*10^exponent;
%    cauchybw_current = cauchybwatthr(alpha_current);
%    if cauchybw_current > bwatthr
%        alpha_left = alpha_current;
%    end
% end
% 
% alpha = alpha_current;
% 
% end

function w = octave_lambertw(b,z)
% Copyright (C) 1998 by Nicol N. Schraudolph <schraudo@inf.ethz.ch>
%
% @deftypefn {Function File} {@var{x} = } lambertw (@var{z})
% @deftypefnx {Function File} {@var{x} = } lambertw (@var{n}, @var{z})
% Compute the Lambert W function of @var{z}.
%
% This function satisfies W(z).*exp(W(z)) = z, and can thus be used to express%
% solutions of transcendental equations involving exponentials or logarithms.%%
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
    % some error checking
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
%%
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

end 
%%error('PRECISION:iteration limit reached, result of lambertw may be inaccurate');
