function [w,s,xvals] = wavfun(g,N,varargin)
% WAVFUN  Wavelet Function
%    Usage: [w,s,xvals] = wavfun(g,N) 
%
%   Input parameters:
%         g     : Wavelet filterbank
%         N     : Number of iterations
%   Output parameters:
%         w     : Approximation of wavelet function(s)
%         s     : Approximation of the scaling function
%         xvals : Correct x-axis values
%
%   Iteratively generate a discrete approximation of wavelet and scaling
%   functions. The algorithm is equal to the DWT reconstruction of a single
%   coefficient at level $N$ set to 1.

definput.keyvals.a = [];
[flags,kv,a]=ltfatarghelper({'a'},definput,varargin);
if(isstruct(g))
    if(isempty(a))
       a = g.a; 
    end
    g=g.g;
end
gLen = length(g);

if(isempty(a))
       a = gLen*ones(gLen,1);
else
    if(length(a)==1)
        a = a*ones(gLen,1); 
    else
       if(length(a)~=gLen)
            error('%s: Number of the subsampling factors is not equal to the number of filters.',upper(mfilename));  
       end
    end
end



lo = g{1}.h(:);
s = lo;
wtemp = cell(gLen-1,1);
for ii=1:gLen-1
    wtemp{ii} = g{ii+1}.h(:);
end

for n=1:N
    for ii=1:gLen-1 
       wtemp{ii} = convolve(ups(wtemp{ii},a(1),1),lo);
    end
    s = convolve(ups(s,a(1),1),lo);
end

w = zeros(length(wtemp{1}),gLen-1);
for ii=1:gLen-1 
   w(:,ii) = wtemp{ii};
end


if(nargout>2)
    xvals=linspace(0,length(lo)-1/length(s),length(s));
end
