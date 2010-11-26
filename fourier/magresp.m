function magresp(g,varargin);
%MAGRESP   Magnitude response plot of window
%   Usage:   magresp(g,...);
%            magresp(g,fs,...);
%            magresp(g,fs,L,....);
%
%   MAGRESP(g) will display the magnitude response of the window g.
%   This is the DFT of g shown on a log scale normalized such that
%   the peak is 0 db.
%
%   MAGRESP(g,fs) does the same, but extends the window to length L.
%   Always use this mode for FIR windows, and select an L somewhat
%   longer than the window to make an accurate plot.
%
%   MAGRESP(g,fs,L) MAGRESP(g,[],L) will do the same for a window
%   intended to be used with signals with sampling rate sr. The x-axis
%   will display Hz.
%
%   If the input window is real, only the positive frequencies will be
%   shown. Adding the option 'nf' as the last parameter will show the
%   negative frequencies anyway.
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
%   Demos: demo_gaborfir     

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

definput.flags.dynrange={'nodynrange','dynrange'};

if isreal(g)
  definput.flags.posfreq={'posfreq','nf'};
else
  definput.flags.posfreq={'nf','posfreq'};
end;

definput.keyvals.fs=[];
definput.keyvals.opts={};
definput.keyvals.L=[];
definput.flags.wintype={'notype','fir','long'};
definput.keyvals.dynrange=100;

[flags,keyvals,fs,L]=ltfatarghelper({'fs','L'},definput,varargin);

[g,info] = comp_fourierwindow(g,L,'MAGRESP');

if flags.do_fir
  info.isfir=1;
end;

if isempty(L) 
  if info.isfir
    L=length(g)*13+47;
  else
    L=length(g);
  end;
end;

g=fir2long(g,L);

% Perform unitaty DFT.
FF=abs(dft(g));

% Scale
FF=FF/max(max(abs(FF)));

% Convert to Db. Add eps to avoid log of zero.
FF=20*log10(FF+realmin);

ymin=max(-keyvals.dynrange,min(min(FF)));

donf=0;
if ~isreal(g)
  donf=1;
end;

if donf
  plotff=fftshift(FF);
  if isempty(fs)
    xrange=-floor(L/2):ceil(L/2)-1;
    axisvec=[-L/2 L/2 ymin 0];
  else
    xrange=linspace(-floor(fs/2),ceil(fs/2)-1,L).';
    axisvec=[-fs/2 fs/2 ymin 0];
  end;
else
  % Only plot positive frequencies for real-valued signals.
  if isempty(fs)
    xrange=linspace(0,1,floor(L/2)+1);
    axisvec=[0 1 ymin 0];
  else
    xrange=linspace(0,floor(fs/2),L/2+1).';
    axisvec=[0 fs/2 ymin 0];
  end;
  plotff=FF(1:floor(L/2)+1);
end;

plot(xrange,plotff,keyvals.opts{:});
axis(axisvec);
ylabel('Magnitude response (Db)');

if isempty(fs)
  xlabel('Frequency (normalized) ');
else
  xlabel('Frequency (Hz)');
end;

legend('off');

