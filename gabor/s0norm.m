function y = s0norm(f,varargin)
%S0NORM S0-norm of signal
%   Usage: y = s0norm(f);
%          y = s0norm(f,...);
%
%   `s0norm(f)` computes the $S_0$-norm of a vector.
%
%   If the input is a matrix or ND-array, the $S_0$-norm is computed along
%   the first (non-singleton) dimension, and a vector of values is returned.
%
%   **WARNING**: The $S_0$-norm is computed by computing a full Short-time
%   Fourier transform of a signal, which can be quite time-consuming. Use
%   this function with care for long signals.
%
%   `s0norm` takes the following flags at the end of the line of input
%   parameters:
%
%     'dim',d   Work along specified dimension. The default value of []
%               means to work along the first non-singleton one.
%
%     'rel'     Return the result relative to the $l^2$ norm (the energy) of the
%               signal.

%   AUTHOR : Peter L. Søndergaard
  
%% ------ Checking of input parameters ---------

if ~isnumeric(f) 
  error('%s: Input must be numerical.',upper(mfilename));
end;

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.flags.rel={'norel','rel'};
definput.keyvals.dim=[];
[flags,kv]=ltfatarghelper({},definput,varargin);

%% ------ Computation --------------------------
 
[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],kv.dim, ...
                                                  upper(mfilename));
permutedsize(1)=1;
y=zeros(permutedsize,assert_classname(f));

g=pgauss(L);

for ii=1:W  
  % Compute the STFT by the simple algorithm and sum each column of the
  % STFT as they are computed, to avoid L^2 memory usage.
  for jj=0:L-1
    y(1,ii)=y(1,ii)+sum(abs(fft(f(:,ii).*circshift(g,jj))));
  end;
  
  if flags.do_rel
    y(1,ii)=y(1,ii)/norm(f(:,ii));
  end;
end;
  
y=y/L;

y=assert_sigreshape_post(y,dim,permutedsize,order);

