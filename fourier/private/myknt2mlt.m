function [m,t] = myknt2mlt(t)

% this is a little exercise how to compute multiplicities in an ordered
% vector of knots
% the copyright is with Matlab's spline toolbox


[r,c] = size(t);

if r*c<2 m=0; 
    return, 
end

difft = diff(t);
if any(difft<0) 
    t = sort(t); 
    difft = diff(t); 
end

index = zeros([r,c]);
index(2:end) = difft==0;
m = cumsum(index);
zz = find(diff(index)<0);

if isempty(zz) 
    return, 
end

z = zeros([r,c]);
pt = m(zz);
pt(2:end) = diff(pt);
z(zz+1) = pt;
m = m - cumsum(z);