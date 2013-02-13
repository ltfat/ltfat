function E = went_compa(x,p)
% 0<p<=2
assert(0<p&&p<=2,'Parameter p have to be in the range ]0,2]');
N = numel(x);
%u = x./norm(x);
v = sort(abs(x(:)),'descend');

wnom = cumsum(v.^p);
wden = sum(v.^p);

w = wnom./wden;

E=N-sum(w);