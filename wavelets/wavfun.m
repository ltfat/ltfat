function [w,s,xvals] = wavfun(g,N)

lo = g{1}(:);
hi = g{2}(:);
s = lo;
w = hi;

for n=1:N
    w = conv(ups(w,2,1),lo);
    s = conv(ups(s,2,1),lo);
end

if(nargout>2)
    xvals=linspace(0,length(lo),length(s));
end
