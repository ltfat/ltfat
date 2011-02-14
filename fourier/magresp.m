function magresp(g,varargin);
%MAGRESP   Magnitude response plot of window
%   Usage:   magresp(g,...);
%            magresp(g,fs,...);
%            magresp(g,fs,dynrange,....);
%
%   MAGRESP(g) will display the magnitude response of the window on a log
%   scale (dB);
%
%   MAGRESP(g,fs) does the same for windows that are intended to be used
%   with signals with sampling rate fs. The x-axis will display Hz.
%
%   MAGRESP(g,fs,dynrange) will limit the dynamic range (see below).
%   
%   MAGRESP takes the following parameters at the end of the line of
%   input arguments.
%
%-     'dynrange',r - Limit the dynamic range of the plot to r dB.
%
%-     'fir'  - Indicate that the input is an FIR window. MAGRESP will
%               zero-extend the window to display a smooth magnitude
%               response.
%
%-     'L',L  - Zero-extend the window to length L.
%
%-     'posfreq' - Show only positive frequencies.
%
%-     'nf'      - Show negative frequencies
%
%-     'autoposfreq' - Show positive frequencies for real-valued signals,
%               otherwise show also the negative frequencies. This is the default.
%
%   In addition to these flags, it is possible to speficy any of the
%   normalization flags from NORMALIZE to normalize the input before
%   calculation of the magnitude response. Specifying '1' or 'area' will
%   display a magnitude response which peaks at 0 dB.
%
%   The following will display the magnitude response of a Hann window
%   of length 20:
%
%C     magresp({'hann',20});
%
%   The following will display the magnitude response of a Gaussian window of length 100
%
%C     magresp('gauss','L',100);
%
%   Demos: demo_gabfir     

%   AUTHOR : Peter Soendergaard.
%   TESTING: NA
%   REFERENCE: NA

if nargin<1
  error('Too few input arguments.');
end;

L=[];
fs=[];
donf=0;

% Define initial value for flags and key/value pairs.

definput.flags.posfreq={'autoposfreq','posfreq','nf'};

definput.import={'normalize'};
definput.importdefaults={'null'};
definput.keyvals.fs=[];
definput.keyvals.opts={};
definput.keyvals.L=[];
definput.flags.wintype={'notype','fir','long'};
definput.keyvals.dynrange=[];

[flags,keyvals,fs,L]=ltfatarghelper({'fs','L'},definput,varargin);

[g,info] = comp_fourierwindow(g,L,'MAGRESP');

do_real=flags.do_posfreq;
if flags.do_autoposfreq
  do_real=isreal(g);
end;

if flags.do_fir
  info.isfir=1;
end;

if isempty(L) 
  if info.isfir
    % Choose a strange length, such that we don't accidentically hit all
    % the zeros in the response.
    L=length(g)*13+47;
  else
    L=length(g);
  end;
end;

g=fir2long(g,L);

g=normalize(g,flags.norm);
if do_real

  % Compute spectrum and normalize
  FF=abs(fftreal(real(g)));
    
  % Convert to Db. Add eps to avoid log of zero.
  FF=20*log10(FF+realmin);

  xmin=0;

else

  % Compute spectrum and normalize. fftshift to center correctly for plotting.
  FF=fftshift(abs(fft(g)));
  
  % Convert to Db. Add eps to avoid log of zero.
  FF=20*log10(FF+realmin);

  xmin=-1;
end;

ymax=max(FF);
if ~isempty(keyvals.dynrange)
  ymin=ymax-keyvals.dynrange;
else
  ymin=min(FF);
end;

Lplot=length(FF);

% Only plot positive frequencies for real-valued signals.
if isempty(fs)
  xrange=linspace(xmin,1,Lplot).';
  axisvec=[xmin 1 ymin ymax];
else
  xrange=linspace(xmin*floor(fs/2),floor(fs/2),Lplot).';
  axisvec=[xmin fs/2 ymin ymax];
end;

plot(xrange,FF,keyvals.opts{:});
axis(axisvec);
ylabel('Magnitude response (Db)');

if isempty(fs)
  xlabel('Frequency (normalized) ');
else
  xlabel('Frequency (Hz)');
end;

legend('off');

