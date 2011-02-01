function b=gammatonefir(fc,fs,varargin);
%GAMMATONEFIR  Gammatone filter coefficients
%   Usage: b = gammatonefir(fc,fs,n,betamul);
%          b = gammatonefir(fc,fs,n);
%          b = gammatonefir(fc,fs);
%
%   Input parameters:
%      fc    -  center frequency in Hz.
%      fs    -  sampling rate in Hz.
%      n     -  filter order.
%      beta  -  bandwidth of the filter.
%
%   Output parameters:
%      b     -  FIR filters as columns
%
%   GAMMATONEFIR(fc,fs,n,betamul) computes the filter coefficients of a digital
%   FIR gammatone filter of length n with center frequency fc, 4th order
%   rising slope, sampling rate fs and bandwith determined by betamul. The
%   bandwidth _beta of each filter is determined as betamul times AUDFILTBW
%   of the center frequency of corresponding filter.
%
%   GAMMATONEFIR(fc,fs,n) will do the same but choose a filter bandwidth
%   according to Glasberg and Moore (1990).  betamul is choosen to be 1.0183.
%
%   GAMMATONEFIR(fc,fs) will do as above and choose a sufficiently long
%   filter to accurately represent the lowest subband channel.
%
%   If fc is a vector, each entry of fc is considered as one center
%   frequency, and the corresponding coefficients are returned as column
%   vectors in the output.
%
%   The inpulse response of the gammatone filter is given by
%
%M    g(t) = a*t^(n-1)*cos(2*pi*fc*t)*exp(-2*pi*beta*t)
%F  \[g(t) = at^{n-1}cos(2\pi\cdot fc\cdot t)e^{-2\pi \beta \cdot t}\]
%
%   The gammatone filters as implemented by this function generate
%   complex valued output, because the filters are modulated by the
%   exponential function. Using REAL on the output will give the
%   coefficients of the corresponding cosine modulated filters.
%
%   To create the filter coefficients of a 1-erb spaced filter bank using
%   gammatone filters use the following construction
%
%C    [b,a] = gammatonefir(fs,erbspacebw(flow,fhigh));
%  
%R  aertsen1980strI glasberg1990daf
  
%   AUTHOR : Peter L. Soendergaard

% ------ Checking of input parameters ---------

if nargin<2
  error('Too few input arguments.');
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

if ~isnumeric(fc) || ~isvector(fc) || any(fc<0) || any(fc>fs/2)
  error(['%s: fc must be a vector of positive values that are less than half ' ...
         'the sampling rate.'],upper(mfilename));
end;

definput.flags.real={'real','complex'};
definput.flags.norm={'equalfreq','equalenergy'};
definput.keyvals.n=[];

definput.keyvals.betamul=1.0183;

[flags,keyvals,n,betamul]  = ltfatarghelper({'n','betamul'},definput,varargin);

nchannels = length(fc);

% ourbeta is used in order not to mask the beta function.

ourbeta = betamul*audfiltbw(fc);

if isempty(n)
  % Calculate a good value for n
  % FIXME actually do this
  n=2000;
end;

b={};

for ii = 1:nchannels

  delay = 3/(2*pi*ourbeta(ii));
  
  scalconst = 2*(2*pi*ourbeta(ii))^4/factorial(4-1)/fs;
  
  nfirst = ceil(fs*delay);
  nlast = n-nfirst;

  t=[(0:nlast-1)/fs+delay,...
     (0:nfirst-1)/fs-nfirst/fs+delay].';  

  % g(t) = a*t^(n-1)*cos(2*pi*fc*t)*exp(-2*pi*beta*t)
  if flags.do_real
    b{ii} = scalconst*t.^(4-1).*cos(2*pi*fc(ii)*t).*exp(-2*pi* ...
                                                      ourbeta(ii)*t);
  else
    b{ii} = scalconst*t.^(4-1).*exp(2*pi*i*fc(ii)*t).*exp(-2*pi* ...
                                                      ourbeta(ii)*t);
  end;    
    
  if flags.do_equalenergy
    b{ii}=b{ii}/rms(b{ii});
  end;
  
  
end;
