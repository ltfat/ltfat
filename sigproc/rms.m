function y = rms(f,varargin)
%RMS RMS value of signal
%   Usage: y = rms(f);
%          y = rms(f,...);
%
%   `RMS(f)` computes the RMS (Root Mean Square) value of a finite sampled
%   signal sampled at a uniform sampling rate. This is a vector norm
%   equal to the $l^2$ averaged by the length of the signal.
%
%   If the input is a matrix or ND-array, the RMS is computed along the
%   first (non-singleton) dimension, and a vector of values is returned.
%
%   The RMS value of a signal *x* of length *N* is computed by
%
%   ..                       N
%      rms(f) = 1/sqrt(N) ( sum |f(n)|^2 )^(1/2)
%                           n=1
%
%   .. math:: rms(f) = \frac{1}{\sqrt N} \left( \sum_{n=1}^N |f(n)|^2
%      \right)^{\frac{1}{2}}
%
%   `RMS` takes the following flags at the end of the line of input
%   parameters:
%
%     'ac'       Consider only the AC component of the signal (i.e. the mean is
%                removed).
%
%     'dim',d    Work along specified dimension. The default value of `[]`
%                means to work along the first non-singleton one.
%

%   AUTHOR : Peter L. Søndergaard
  
%% ------ Checking of input parameters ---------

if ~isnumeric(f) 
  error('%s: Input must be numerical.',upper(mfilename));
end;

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.dim=[];
definput.flags.mean={'noac','ac'};
[flags,kv]=ltfatarghelper({'dim'},definput,varargin);

%% ------ Computation --------------------------

% It is better to use 'norm' instead of explicitly summing the squares, as
% norm (hopefully) attempts to avoid numerical overflow.
 
[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],kv.dim, ...
                                                  upper(mfilename));
permutedsize(1)=1;
y=zeros(permutedsize);
if flags.do_ac

  for ii=1:W        
    y(1,ii) = norm(f(:,ii)-mean(f(:,ii)))/sqrt(L);
   end;

else

  for ii=1:W
    y(1,ii)=norm(f(:,ii))/sqrt(L);
  end;

end;
  
y=assert_sigreshape_post(y,kv.dim,permutedsize,order);
