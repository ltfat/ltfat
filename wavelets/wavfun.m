function [w,s,xvals] = wavfun(g,N)
% WAVFUN  Wavelet Function
%    Usage: [w,s,xvals] = wavfun(g,N) 
%
%   Input parameters:
%         g     : Wavelet filterbank
%         N     : Number of iterations
%   Output parameters:
%         w     : Approximation of wavelet function(s)
%         s     : Approximation of the scaling function
%         xvals : Correct 
%
%   Iteratively generate a discrete approximation of wavelet and scaling functions.  
%


if(isstruct(g))
    g=g.g;
end

gLen = length(g);

lo = g{1}(:);
s = lo;
wtemp = cell(gLen-1,1);
for ii=1:gLen-1
    wtemp{ii} = g{ii+1}(:);
end

for n=1:N
    for ii=1:gLen-1 
       wtemp{ii} = convolve(ups(wtemp{ii},2,1),lo);
    end
    s = convolve(ups(s,2,1),lo);
end

w = zeros(length(wtemp{1}),gLen-1);
for ii=1:gLen-1 
   w(:,ii) = wtemp{ii};
end


if(nargout>2)
    xvals=linspace(0,length(lo),length(s));
end
