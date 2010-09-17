function xo=uquant(xi,nbits,xmax,varargin);
%UQUANT  Simulate uniform quantization.
%   Usage:  x=uquant(x,nbits,xmax);
%           x=uquant(x,nbits,xmax,...);
%
%   UQUANT(x,nbits,xmax) simulates the effect of uniform quantization of x using
%   nbits bit. The output is simply x rounded to $2^{nbits}$ different values.
%   The xmax parameters specify the maximal value that should be quantifiable.
%
%   UQUANT takes the following flags at the end of the input arguments.
%
%-   's' - Signed quantization. This assumes that the signal
%            has a both positive and negative part. Usefull for sound
%            signals. This is the default
%
%-   'u' - Unsigned quantization. Assumes the signal is positive.
%          Negative values are silently rounded to zero.
%          Usefull for images.
%
%   If this function is applied to a complex signal, it will simply be
%   applied to the real and imaginary part separately.
%

%   AUTHOR : Peter Soendergaard and Bruno Torresani.  
%   TESTING: OK
%   REFERENCE: OK

if nargin<4
  error('Too few input parameters.');
end;

% Define initial value for flags and key/value pairs.
definput.flags.sign={'s','u'};

[flags,keyvals]=ltfatarghelper({},definput,varargin);

% ------ handle complex values ------------------
if ~isreal(xi)
  xo = uquant(real(xi),nbits,xmax,varargin{:}) + ...
	i*uquant(imag(xi),nbits,xmax,varargin{:});
  return
end;

if nbits<1
  error('Must specify at least 2 bits.');
end;

% Calculate number of buckets.
nbuck=2^nbits;    

if xmax<max(abs(xi(:)))
  error('Signal contains values higher than xmax.');
end;

if flags.do_s
  
  % ------------ unsigned case -----------------
  
  bucksize=xmax/(nbuck/2-1);
  
  xo=round(xi/bucksize)*bucksize;        
  
else

  % ------------- signed case------------
  
  bucksize=xmax/(nbuck-.51);
  
  % Thresh all negative values to zero.
  xi=xi.*(xi>0);
  
  xo=round(xi/bucksize)*bucksize;        
  
end;

