function [g,info]=firwin(name,M,varargin);
%FIRWIN  FIR window
%   Usage:  g=firwin(name,M);
%           g=firwin(name,M,...);
%
%   FIRWIN(name,M) will return a FIR window of length M of type name.
%
%   All windows are symmetric, they are zero delay and zero phase
%   filters. They can be used for the Wilson and WMDCT transform, except when noted otherwise.
%
%   In the following PSL means "Peak Sidelobe level", and the main lobe
%   width is measured in normalized frequencies.
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
%   The windows available are:
%
%-      hann       - von Hann window. Forms a PU. The Hann window has a
%                    mainlobe with of 8/M, a PSL of -31.5 dB and decay rate
%                    of 18 dB/Octave.
%
%-      sine       - Sine window. This is the square root of the Hanning
%                    window. The sine window has a mainlobe width of 8/M,
%                    a  PSL of -22.3 dB and decay rate of 12 dB/Octave.
%                    Aliases: 'cosine', 'sqrthann'
%
%-      rect       - (Almost) rectangular window. The rectangular window has a
%                    mainlobe width of 4/M, a  PSL of -13.3 dB and decay
%                    rate of 6 dB/Octave. Forms a PU. Alias: 'square'
%
%-      sqrtrect   - Square root of the rectangular window.
%
%-      tria       - (Almost) triangular window. Forms a PU. Alias: 'bartlett'
%
%-      sqrttria   - Square root of the triangular window.
%
%-      hamming    - Hamming window. Forms a PU that sums to 1.08 instead
%                    of 1.0 as usual. The Hamming window has a
%                    mainlobe width of 8/M, a  PSL of -42.7 dB and decay
%                    rate of 6 dB/Octave. This window should not be used for
%                    a Wilson basis, as a reconstruction window cannot be
%                    found by WILDUAL.
%
%-      blackman   - Blackman window. The Hamming window has a
%                    mainlobe width of 12/M, a PSL of -58.1 dB and decay
%                    rate of 18 dB/Octave.
%
%-      blackman2  - Alternate Blackman window. This window has a
%                    mainlobe width of 12/M, a PSL of -68.24 dB and decay
%                    rate of 6 dB/Octave.
%
%-      nuttall    - Nuttall window. The Nuttall window has a
%                    mainlobe width of 16/M, a PSL of -93.32 dB and decay
%                    rate of 18 dB/Octave.
%
%-      itersine -   Iterated sine window. Generates an orthonormal
%                    Wilson/WMDCT basis. This window is described in 
%                    Wesfreid and Wickerhauser (1993) and is used in  the
%                    ogg sound codec. Alias: 'ogg'
%
%   FIRWIN understands the following flags at the end of the list of input
%   parameters:
%
%-     'wp'      - Output is whole point even. This is the default.
%
%-     'hp'      - Output is half point even, as most Matlab filter
%                  routines.
%
%-     'taper',t - Extend the window by a flat section in the middle. The
%                  argument t is the ratio of the rising and falling
%                  parts as compared to the total length of the
%                  window. The default value of 1 means no
%                  tapering. Accepted values lie in the range from 0 to 1.
%
%   Additionally, FIRWIN accepts flags to normalize the output. Please see
%   the help of NORMALIZE. Default is to use 'peak' normalization, which is
%   useful for using the output for windowing in the time-domain. For
%   filtering in the time-domain, a normalization of '1' or 'area' is
%   preferable.
%
%   See also:  pgauss, pbspline, firkaiser, normalize
%
%R  opsc89 harris1978 nuttall1981 wesfreid1993
 
%   AUTHOR : Peter Soendergaard.
%   REFERENCE: NA
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~ischar(name)
  error('%s: First input argument must the name of a window.',upper(mfilename));
end;
  

% Always set to this
info.isfir=1;

% Default values, may be overwritten later in the code
info.ispu=0;
info.issqpu=0;

name=lower(name);

% Define initial value for flags and key/value pairs.
definput.import={'normalize'};
definput.importdefaults={'null'};
definput.flags.centering={'wp','hp'};
%definput.flags.delay={'nodelay','delay','causal'};

definput.keyvals.taper=1;
%definput.keyvals.delay=0;

[flags,keyvals]=ltfatarghelper({},definput,varargin);

if flags.do_wp
  cent=0;
else
  cent=0.5;
end;

if M==0
  g=[];
  return;
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
x   = ((0:M-1)'+cent)/M;

% Use this if window length must be odd or even
Modd   = M-mod(M+1,2);
Meven  = M-mod(M,2);
xodd   = ((0:Modd-1)'+cent)/Modd;
xeven  = ((0:Meven-1)'+cent)/Meven;
do_sqrt=0;
switch name    
 case {'hanning','hann'}
  g=(0.5+0.5*cos(2*pi*x));
  info.ispu=1;
  
 case {'sine','cosine','sqrthann'}
  g=firwin('hanning',M,varargin{:});
  info.issqpu=1;
  do_sqrt=1;
  
 case {'hamming'}  
  if cent==0
    g=0.54-0.46*cos(2*pi*(xodd-.5));
  else
    g=0.54-0.46*cos(2*pi*(xeven-.5));
  end;
  
 case {'square','rect'} 
  if cent==0
    g=ones(Modd,1);
  else
    g=ones(Meven,1);
  end;

 case {'halfsquare','halfrect'} 
  g=ones(Meven/2,1);
  info.ispu=1;
  
 case {'sqrthalfsquare','sqrthalfrect'}
  g=ones(Meven/2,1);
  info.issqpu=1;
  do_sqrt=1;
  
 case {'tria','triangular','bartlett'}
  if flags.do_wp
    gw=linspace(1,0,Meven/2+1).';
    g=[gw;flipud(gw(2:end-1))];
  %g=pbspline(M,1,Meven/2,flags.centering);
  else
    gw=((0:Meven/2-1)/(Meven/2)+0.5).';
    g=[gw;flipud(gw)];
  end;
  info.ispu=1;
  
 case {'sqrttria'}
  g=firwin('tria',M,flags.centering);
  info.issqpu=1;
  do_sqrt=1;
  %if flags.do_hp && rem(M,2)==1
  %  % Remove small error in this case
  %  g((M+1)/2)=0;
  %end;
  
 case {'blackman'}
  g=0.42-0.5*cos(2*pi*(x-.5))+0.08*cos(4*pi*(x-.5));

 case {'blackman2'}
  g=7938/18608-9240/18608*cos(2*pi*(x-.5))+1430/18608*cos(4*pi*(x-.5));

 case {'nuttall'}
  g = 0.355768 - 0.487396*cos(2*pi*(x-.5)) + 0.144232*cos(4*pi*(x-.5)) -0.012604*cos(6*pi*(x-.5));
  
 case {'itersine','ogg'}
  g=sin(pi/2*sin(pi*(x-.5)).^2);
  info.issqpu=1;
  
 otherwise
  error('Unknown window: %s.',name);
end;

% Add zeros if needed.
if length(g)<M
  g=middlepad(g,M,flags.centering);
end;

% Do sqrt if needed. This must be done /after/ middlepad, so do it at
% this point.
if do_sqrt
  g=sqrt(g);
end;

if keyvals.taper<1
  % Perform the actual tapering.
  g=[ones(p1,1);g;ones(p2,1)];  
end;

g=normalize(g,flags.norm);
