function []=phaseplot(f,varargin)
%PHASEPLOT  Phase plot
%   Usage: phaseplot(f,op1,op2, ... );
%          phaseplot(f,fs,op1,op2, ... );
%
%   PHASEPLOT(f) plots the phase of f using a DGT.
%
%   PHASEPLOT(f,fs) does the same for a signal with sampling rate fs
%   (sampled with fs samples per second);
%
%   PHASEPLOT should only be used for short signals (shorter than the
%   resolution of the screen), as there will otherwise be some visual
%   aliasing, such that very fast changing areas will look very smooth.
%   PHASEPLOT always calculates the phase of the full time/frequency plane
%   (as opposed to SGRAM), and you therefore risk running out of memory
%   for long signals.
%
%   Additional arguments can be supplied like this:
%   PHASEPLOT(f,'nf','tfr',2,'phl'). The arguments must be character strings
%   possibly followed by an argument:
%
%-  'tfr',v   - Set the ratio of frequency resolution to time resolution.
%               A value v=1 is the default. Setting v>1 will give better
%               frequency resolution at the expense of a worse time
%               resolution. A value of 0<v<1 will do the opposite.
%
%-  'p'       - Use a periodic bounday condition. This is the default.
%
%-  'e'       - Use an even boundary condition. The produces the fewest
%               visible artifacts at the boundary. Note: Currently this
%               doubles the computational time.
%
%-  'nf'      - Display negative frequencies, with the zero-frequency
%               centered in the middle. For real signals, this will just
%               mirror the upper half plane. This is standard for complex
%               signals.
%
%-  'tc'      - Time centering. Move the beginning of the signal to the
%               middle of the plot. This is usefull for visualizing the
%               window functions of the toolbox.
%
%-  'thr'     - Use the amplitude values to threshold the phase values.
%               For small amplitude values the phase values can be meaningless.
%               This is because an error, that is small in amplitude can still
%               lead to arbitrary phase values in the vicinity of zero.
%               Setting this flag will set the phase is to a constant value
%               (0) for all pairs (m,n), where the amplitude is below a 
%               certain relative value, set to 0.001 of the maximum amplitude.
%
%-  'nophl'   - By default the phase is displayed as it would be if calculated
%               by a traditional filter bank implementation (a time
%               invariant Gabor system). This is called a phase-locked
%               dgt. This switch will instead display the phase as it is
%               calculated by the DGT (a frequency invariant Gabor system).
%
%-  'fmax',y  - Display y as the highest frequency.
% 
%   See also: phaselock
%
%   Demos: demo_phaseplot
%
%R  Carmona98practical

%   AUTHOR: Peter Soendergaard
%   REFERENCE: NA
%   TESTING: NA

if nargin<1
  error('Too few input arguments.');
end;

if sum(size(f)>1)>1
  error('Input must be a vector.');
end;

% Standard values controlled by optional arguments
if isreal(f)
  donf=0;
else
  donf=1;
end;

tfr_mul=1;
dotc=0;
domask=0;
dophaslock=1;
dofs=0;
doeven=0;
dosymphl=0;
dofmax=0;

% Resampling rate: Used when fmax is issued.
resamp=1;

% Parse optional arguments
if ~isempty(varargin)

  startii=1;

  if isnumeric(varargin{1})
    dofs=1;
    fs=varargin{1};
    startii=2;
  end;

  ii=startii-1;
  while ii<length(varargin)

    ii=ii+1;

    optarg=varargin{ii};

    if ischar(optarg)
      switch lower(optarg)
	case 'nf'
	  donf=1;
	case 'tc'
	  dotc=1;
	  f=fftshift(f);
	case 'p'
	case 'e'
	  doeven=1;
    case 'thr'
      domask=1;
      % mask_val=varargin{ii+1}; Old version, not fitting description
      mask_val=0.001;
      ii=ii+1;
    case 'fmax'
      dofmax=1;
      fmax=varargin{ii+1};
      ii=ii+1;
    case 'nophl'
      dophaslock=0;          
    case 'phl'
      dophaslock=1;
    case 'symphl'
      dosymphl=1;
    otherwise
	  error([optarg,' : Unknown optional argument 1']);
      end;
    end;

    if iscell(optarg)
      if isempty(optarg) || ~ischar(optarg{1})
	error('First element of optinal argument cell array must be a character string.');
      end;

      switch lower(optarg{1})
	case 'tfr'
	  tfr_mul=optarg{2};
	otherwise
	  error([optarg{1},' : Unknown optional argument 2']);
      end;

    end;
  end;
  
end;

% Downsample
if dofmax
  if dofs
    resamp=fmax*2/fs;
  else
    resamp=fmax*2/length(f);
  end;
  
  f=fftresample(f,round(length(f)*resamp));
end;


% Always do the full STFT
L=length(f);
a=1;
b=1;
M=L;
N=L;

if doeven
  if mod(a,2)==0
    g=pgauss(2*L,.5*tfr_mul,.5);
  else
    g=pgauss(2*L,.5*tfr_mul,0);
  end;
else
  g=pgauss(L,tfr_mul);
end;

if doeven
  f=[f;flipud(f)];
  eL=length(f);
  f=f.*(exp(-2*pi*i*(eL-1)).');
  g=circshift(g,floor(a/2));
  coef=dgt(f,g,a,M);
  coef=coef(:,1:size(coef,2)/2);

  % Correct the phase.
  if dophaslock
    % Correct the phase.
    TimeInd = (0:(N-1))*a+.5;
    FreqInd = (0:(M-1)+.5)/M;
    
    phase = FreqInd'*TimeInd;
    phase = exp(-2*i*pi*phase);
    coef=coef.*phase;

    %for n=0:size(coef,2)-1
    %  for m=0:M-1
%	coef(m+1,n+1)=coef(m+1,n+1)*exp(-2*pi*i*(.5-n*a)*(m+.5)/M);
%      end;
 %   end;
  else
    phasecor=exp(-2*pi*i*.5*((0:M-1).'+.5)/M);
    coef=coef.*repmat(phasecor,1,size(coef,2));
  end;

else
  coef=dgt(f,g,a,M);

  if dophaslock
    % use phase locking
    coef = phaselock(coef,a);
  end
  
  if dosymphl
    TimeInd = (0:(N-1))/N;
    FreqInd = (0:(M-1))*b;
    
    phase = FreqInd'*TimeInd;
    phase = exp(i*pi*phase);
    coef=coef.*phase;
  end;

end;

if domask
  % keep only the largest coefficients.
  maxc=max(abs(coef(:)));
  mask=abs(coef)<maxc*mask_val;
  coef(mask)=0;
end

coef = angle(coef);

if donf
  % Calculate negative frequencies, use DGT
  % Display zero frequency in middle of plot.
 
  % Move zero frequency to the center.
  coef=fftshift(coef,1);

  if dofmax
    yr=-fmax:fmax/M:fmax;
  else
    if dofs
      yr=-fs/2:fs/M:fs/2;
    else
      yr=-L/2:b:L/2;
    end;
  end;

else
  % Dont calculate negative frequencies, try to use DGT of twice the size.
  coef=coef(1:ceil((M+1)/2),:);
  if dofmax
    yr=0:fmax/M:fmax;
  else
    if dofs
      yr=0:fs/M:fs/2;
    else
      yr=0:b:L/2;
    end;
  end;
end;

if dotc
  xr=-floor(N/2)*a:a:floor((N-1)/2)*a;
else
  xr=0:a:N*a-1;
end;

if dofs
  % Scale x-axis by sampling rate.
  xr=xr/fs;  
end;

imagesc(xr,yr,coef);
axis('xy');
if dofs
  xlabel('Time (s)')
  ylabel('Frequency (Hz)')
else
  xlabel('Time')
  ylabel('Frequency')
end;

