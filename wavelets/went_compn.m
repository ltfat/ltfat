function E = went_compn(x,p,f)
% 0<p<=2
assert(0<p&&p<=2,'Parameter p have to be in the range ]0,2]');
assert(0<f&&f<1,'Parameter f  have to be in the range ]0,1[');

%u = x./norm(x);
v = sort(abs(x(:)),'descend');

wnom = cumsum(v.^p);
wden = sum(v.^p);

w = wnom./wden;

[~,E]=min(abs(w-f));