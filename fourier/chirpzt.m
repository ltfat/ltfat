function c = chirpzt(f,K,deltao,o,dim)
%CHIRPZT Chirped Z-transform
%   Usage:  c = chirpzt(f,o0,deltao,K,dim)
%
%   Input parameters:
%         f      : Input data.
%         K      : Number of values.
%         deltao : Angle increment.
%         o      : Starting angle. 
%
%   Output parameters:
%         c      : Coefficient vector.
%
%   `c = chirpzt(f,o0,deltao,K)` computes *K* samples of the discrete-time 
%   fourier transform DTFT *c* of *f* at values $c(k)=F(o0+k*deltao)$ 
%   for $k=0,\dots,K-1$. The values are computed along dimension *dim*.

%% Check the input arguments
if nargin < 1
    error('%s: Not enough input arguments.',upper(mfilename))
end

if isempty(f)
    error('%s: X must be a nonempty vector or a matrix.',upper(mfilename))
end

if nargin<5
  dim=[];  
end;

[f,~,Ls,~,dim,permutedsize,order]=assert_sigreshape_pre(f,[],dim,'CHIRPZT');

if nargin > 1  && ~isempty(K)
   if ~isreal(K) || ~isscalar(K)
      error('%s: o0 must be a real scalar.',upper(mfilename))
   end
else
   K = size(f,1);
end

if nargin > 2  && ~isempty(deltao)
   if ~isreal(K) || ~isscalar(K)
      error('%s: deltao must be a real scalar.',upper(mfilename))
   end
else
   deltao = 2*pi/K;
end

if nargin > 3  && ~isempty(o)
   if ~isreal(K) || ~isscalar(K)
      error('%s: o must be a real scalar.',upper(mfilename))
   end
else
   o = 0;
end

c = comp_chirpzt(f,K,deltao,o);


permutedsize(1)=K;
c=assert_sigreshape_post(c,dim,permutedsize,order);

