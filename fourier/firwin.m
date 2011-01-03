function g=firwin(name,M,varargin);
%FIRWIN  FIR window
%   Usage:  g=firwin(name,M);
%           g=firwin(name,M,...);
%
%   FIRWIN(name,M) will return a FIR window of length M of type name.
%
%   All windows are symmetric, they are zero delay and zero phase
%   filters. They can be used for Wilson bases, except when noted otherwise.
%
%   If a window g forms a "partition of unity" (PU) it means specifically
%   that
%
%C         g+fftshift(g)==ones(L,1);
%
%   A perfect PU can only be formed if the window length is even, but
%   some windows may work for odd lengths anyway.
%
%   If a window is the square root of a window that forms a PU, the window
%   will generate a tight Gabor frame / orthonormal Wilson/WMDCT basis if
%   the number of channels is less than M.
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
%-      hann       - Hann window. Forms a PU.
%
%-      sine       - Sine window. This is the square root of the Hanning
%                    window. Aliases: 'cosine', 'sqrthann'
%
%-      square     - (Almost) rectangular window. Forms a PU. Alias: 'rect'
%
%-      sqrtsquare - Square root of the square window.
%
%-      tria -       (Almost) triangular window. Forms a PU. Alias: 'bartlett'
%
%-      sqrttria   - Square root of the triangular window.
%
%-      hamming    - Hamming window. Forms a PU that sums to 1.08 instead
%                    of 1.0 as usual. This window should not be used for
%                    a Wilson basis, as a reconstruction window cannot be
%                    found by WILDUAL.
%
%-      sqrtham    - Square root of a Hamming window.
%
%-      blackman   - Blackman window
%
%-      nuttall    - Nuttall window
%
%-      ogg        - Iterated sine window used in the ogg sound
%                    Generates an ortonormal Wilson/WMDCT basis.
% 
%   See also:  pgauss, pbspline, firkaiser
%
%R  opsc89 harris1978 nuttall1981
  
  
%     Normalized so it
%                    generates a normalized tight Gabor system with parameters
%                    a=M/2 and M or an orthonormal Wilson/WMDCT basis with
%                    M channels.

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

definput.keyvals.taper=1;
definput.keyvals.delay=0;

[flags,keyvals]=ltfatarghelper({},definput,varargin);

if flags.do_wp
  cent=0;
else
  cent=0.5;
end;

% Deal with tapering
if keyvals.taper<1
  if keyvals.taper==0
    % Window is only tapering, do this and bail out, because subsequent
    % code may fail.
    g=ones(M,1);
    return;
  end;

  % Modify M to fit with tapering
  Morig=M;
  M=round(M*keyvals.taper);
  Mtaper=Morig-M;

  p1=round(Mtaper/2);
  p2=Mtaper-p1;

  % Switch centering if necessary
  if cent==0 && p1~=p2
    cent=0.5;
  end;
  
  if cent==0.5 && p1~=p2
    cent=1;
  end;
    
end;

% This is the normally used sampling vector.
x=((0:M-1)'+cent)/M;
switch name    
 case {'hanning','hann'}
  g=(0.5+0.5*cos(2*pi*x));
  
 case {'sine','cosine','sqrthann'}
  g=sqrt(firwin('hanning',M,varargin{:}));
  
 case {'hamming'}
  g=0.54-0.46*cos(2*pi*(x-.5));
  
 case {'sqrtham','sqrthamming'}
  g=sqrt(firwin('hamming',M,varargin{:})/1.08);
  
 case {'square','rect'} 
  g=middlepad(ones(M/2,1),M,flags.centering);
  
 case {'sqrtsquare','sqrtrect'}
  g=sqrt(middlepad(ones(M/2,1),M,flags.centering));

 case {'tria','triangular','bartlett'}
  g=pbspline(M,1,M/2,flags.centering);

 case {'sqrttria'}
  g=sqrt(pbspline(M,1,M/2,flags.centering));
  
 case {'blackman'}
  g=0.42-0.5*cos(2*pi*(x-.5))+0.08*cos(4*pi*(x-.5));

 case {'nuttall'}
  g = 0.355768 - 0.487396*cos(2*pi*(x-.5)) + 0.144232*cos(4*pi*(x-.5)) -0.012604*cos(6*pi*(x-.5));
  
 case {'ogg'}
  g=sin(pi/2*sin(pi*(x-.5)).^2);
  
 otherwise
  error('Unknown window: %s.',name);
end;

if keyvals.taper<1
  % Perform the actual tapering.
  g=[ones(p1,1);g;ones(p2,1)];  
end;

