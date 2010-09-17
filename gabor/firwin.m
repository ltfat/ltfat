function g=firwin(name,M,varargin);
%FIRWIN  FIR window
%   Usage:  g=firwin(name,M);
%           g=firwin(name,M,...);
%
%   FIRWIN(name,M) will return a FIR window of length M of type
%   name.
%
%   The windows are normalized, such that if used for Gabor systems
%   with parameters _a=M/2 and M, they will generate Gabor frames
%   with framebounds close to 1.
%
%   All windows are symmetric and can be used for Wilson bases, except
%   when noted otherwise.
%
%   If a window g forms a "partition of unity" (PU) it means specifically
%   that
%
%C         g+fftshift(g)==ones(L,1);
%
%   A perfect PU can only be formed if the window length is even, but
%   some windows may work for odd lengths anyway.
% 
%   FIRWIN understands the following flags at the end of the list of input
%   parameters:
%
%-     'wp'      - Output is whole point even. This is the default.
%
%-     'hp'      - Output is half point even, as most Matlab filter
%                  routines.
%
%-     'delay',d - Delay the output by d samples. Default is zero delay.
%
%-     'causual' - Delay the window enough to make it causal.
%
%   The windows available are:
%
%-      hanning    - Hanning window. Forms a PU.
%
%-      sqrthann   - Square root of a Hanning window. Normalized so it
%                    generates a normalized tight Gabor system with parameters
%                    a=M/2 and M or an orthonormal Wilson/WMDCT basis with
%                    M channels.
%
%-      square     - (Almost) rectangular window. Forms a PU. Alias: 'rect'
%
%-      sqrtsquare - Square root of the square window. As sqrthann.
%
%-      tria - (Almost) triangular window. Forms a PU. Alias: 'bartlett'
%
%-      sqrttria   - Square root of the triangular window. As sqrthann.
%
%-      hamming    - Hamming window. Forms a PU that sums to 1.08 instead
%                    of 1.0 as usual. This window should not be used for
%                    a Wilson basis, as a reconstruction window cannot be
%                    found by WILDUAL.
%
%-      sqrtham    - Square root of a Hamming window. As sqrthan.
%
%-      blackman   - Blackman window
%
%-      ogg        - Iterated sine window used in the ogg sound
%                    Generates an ortonormal Wilson/WMDCT basis.
% 
%   See also:  pgauss, pbspline, firkaiser
%
%R  opsc89

%   AUTHOR : Peter Soendergaard.
%   REFERENCE: NA
  
if nargin<2
  error('Too few input parameters.');
end;

name=lower(name);

if rem(M,2)==1
  % Some windows are not defined for odd lengths.
  
  switch name
 case {'square','rect','sqrtsquare','sqrtrect','tria','sqrttria','sqrtblack'}
    error('The length of the choosen window must be even.');
  end;
end;

% Define initial value for flags and key/value pairs.
definput.flags.centering={'wp','hp'};
definput.flags.delay={'nodelay','delay','causal'};

definput.keyvals.delay=0;

[flags,keyvals]=ltfatarghelper({},definput,varargin);

if flags.do_wp
  cent=0;
else
  cent=0.5;
end;

% This is the normally used sampling vector.
x=((0:M-1)'+cent)/M;
switch name    
 case {'hanning','hann'}
  g=(0.5+0.5*cos(2*pi*x));
  
 case {'sqrthann'}
  g=sqrt(firwin('hanning',M,varargin{:}))/sqrt(M/2);
  
 case {'hamming'}
  g=0.54-0.46*cos(2*pi*(x-.5));
  
 case {'sqrtham'}
  g=sqrt(2)*sqrt(firwin('hamming',M,varargin{:})/1.08)/sqrt(M);
  
 case {'square','rect'} 
  g=middlepad(ones(M/2,1),M,flags.centering);
  
 case {'sqrtsquare','sqrtrect'}
  g=sqrt(middlepad(ones(M/2,1),M,flags.centering))/sqrt(M/2);

 case {'tria','triangular','bartlett'}
  g=pbspline(M,1,M/2,flags.centering)*sqrt(M/2);

 case {'sqrttria'}
  g=sqrt(pbspline(M,1,M/2,flags.centering)/sqrt(M/2));
  
 case {'blackman'}
  g=0.42-0.5*cos(2*pi*(x-.5))+0.08*cos(4*pi*(x-.5));

 case {'ogg'}
  g=sin(pi/2*sin(pi*(x-.5)).^2)/sqrt(M/2);
  
 case {'testfun'}
  g=fftshift(exp(-1./(1/10*x.*(1-x)))/exp(-4));  
 otherwise
  error('Unknown window: %s.',name);
end;
